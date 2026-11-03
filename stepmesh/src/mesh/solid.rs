//! One B-rep solid as a closed, oriented triangle mesh.

use std::collections::{HashMap, HashSet};

use super::edges::discretize_edge;
use super::face::{FaceMesh, VRef, clip_ears, tessellate};
use crate::fxhash::FxMap;
use crate::geom::Surface;
use crate::geom::bounds::has_spline;
use crate::model::{Brep, Face, Loop, P3, cross, dot, norm, sub};
use crate::{Error, Options};

type Polylines = HashMap<u64, Vec<P3>>;

/// Times the samples of an edge are made finer (each halving the chord tolerance) when its
/// polyline crosses another in a face's parameter space.
const MAX_EDGE_REFINEMENTS: u32 = 8;

/// Vertex identity across the faces of a solid: edge end points are the B-rep's vertices.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Key {
    Vertex(u64),
    EdgeSample(u64, usize),
    Inner(usize, usize),
}

/// Mesh in the solid's own coordinates.
pub struct LocalMesh {
    pub positions: Vec<P3>,
    pub triangles: Vec<[u32; 3]>,
}

/// Meshes every face of `b`; `solid` is the index reported in errors. `budget` is the number of
/// triangles still allowed.
pub fn mesh_brep(
    b: &Brep,
    o: &Options,
    solid: usize,
    budget: Option<usize>,
) -> Result<LocalMesh, Error> {
    let (polys, meshes) = mesh_faces(b, o)?;
    let joined = joined_vertices(b);
    let vertex = |v: u64| joined.get(&v).copied().unwrap_or(v);
    let mut index: FxMap<Key, u32> = FxMap::default();
    let mut positions: Vec<P3> = Vec::new();
    let mut triangles = Vec::new();
    let mut face_no = 0;
    let mut meshes = meshes.into_iter();
    for shell in &b.shells {
        let mut shell_tris = Vec::new();
        for fm in meshes.by_ref().take(shell.faces.len()) {
            let mut id = |r: VRef| -> u32 {
                let (key, p) = match r {
                    VRef::Edge { edge, i } => {
                        let e = &b.edges[&edge];
                        let last = polys[&edge].len() - 1;
                        let key = match i {
                            0 => Key::Vertex(vertex(e.start_vertex)),
                            i if i == last => Key::Vertex(vertex(e.end_vertex)),
                            i => Key::EdgeSample(edge, i),
                        };
                        (key, polys[&edge][i])
                    }
                    VRef::Inner(i) => (Key::Inner(face_no, i), fm.inner[i]),
                };
                *index.entry(key).or_insert_with(|| {
                    positions.push(p);
                    (positions.len() - 1) as u32
                })
            };
            for t in &fm.tris {
                let [a, b, c] = t.map(&mut id);
                if a != b && b != c && a != c {
                    shell_tris.push(if shell.reversed { [a, c, b] } else { [a, b, c] });
                }
            }
            face_no += 1;
            if let Some(limit) = budget {
                if triangles.len() + shell_tris.len() > limit {
                    return Err(Error::TriangleLimit {
                        count: triangles.len() + shell_tris.len(),
                        limit,
                    });
                }
            }
        }
        let v = super::check::signed_volume(&positions, &shell_tris);
        let wrong_sign = if shell.void { v >= 0.0 } else { v <= 0.0 };
        if !super::check::closed_oriented(&shell_tris) || wrong_sign {
            return Err(Error::NotWatertight { solid });
        }
        triangles.extend(shell_tris);
    }
    Ok(LocalMesh {
        positions,
        triangles,
    })
}

/// Edges of spline faces. Such a face's certified bound includes its boundary chords, which the
/// face cannot refine: their edges leave it half the deflection.
fn spline_edges(b: &Brep) -> HashSet<u64> {
    b.shells
        .iter()
        .flat_map(|s| &s.faces)
        .filter(|f| has_spline(&f.surface))
        .flat_map(|f| &f.loops)
        .filter_map(|(l, _)| match l {
            Loop::Edges(list) => Some(list.iter().map(|&(e, _)| e)),
            _ => None,
        })
        .flatten()
        .collect()
}

/// Sampling options of an edge, `level` times made finer.
fn edge_options(o: &Options, on_spline: bool, level: u32) -> Options {
    let share = if on_spline { 2. } else { 1. };
    Options {
        linear_deflection: o.linear_deflection / (share * f64::from(1u32 << level)),
        ..o.clone()
    }
}

/// The polyline of every edge of `b`, sampled as its faces need.
pub fn edge_polylines(b: &Brep, o: &Options) -> Result<Polylines, Error> {
    let spline = spline_edges(b);
    b.edges
        .iter()
        .map(|(&id, e)| {
            Ok((
                id,
                discretize_edge(e, &edge_options(o, spline.contains(&id), 0))?,
            ))
        })
        .collect()
}

/// Edge polylines and the mesh of every face, in shell order. Edges sampled to the tolerances
/// can still cross in a face's parameter space where two curves run closer than the tolerance
/// (a wall thinner than it, a face narrower than it): such edges are sampled again, finer, and
/// the faces using them meshed again.
fn mesh_faces(b: &Brep, o: &Options) -> Result<(Polylines, Vec<FaceMesh>), Error> {
    let faces: Vec<&Face> = b.shells.iter().flat_map(|s| &s.faces).collect();
    let spline = spline_edges(b);
    let mut polys = edge_polylines(b, o)?;
    let mut meshes: Vec<Option<FaceMesh>> = vec![None; faces.len()];
    let mut refinements: HashMap<u64, u32> = HashMap::new();
    // Faces before `next` are meshed.
    let mut next = 0;
    while let Some(k) = (next..faces.len()).find(|&k| meshes[k].is_none()) {
        next = k;
        let (error, crossing) = match tessellate(faces[k], &polys, o) {
            Ok(fm) => {
                meshes[k] = Some(fm);
                continue;
            }
            Err(failure) => failure,
        };
        if crossing.is_empty() {
            return Err(error);
        }
        // Curves that cross within the tolerance stay crossed however fine their samples; on a
        // plane the loop is cut into ears instead.
        let spent = crossing.iter().any(|id| {
            refinements
                .get(id)
                .is_some_and(|&n| n >= MAX_EDGE_REFINEMENTS)
        });
        if spent {
            if !matches!(faces[k].surface, Surface::Plane(_)) {
                return Err(error);
            }
            meshes[k] = Some(clip_ears(faces[k], &polys, o)?);
            continue;
        }
        for &id in &crossing {
            let level = refinements.entry(id).or_default();
            *level += 1;
            polys.insert(
                id,
                discretize_edge(
                    &b.edges[&id],
                    &edge_options(o, spline.contains(&id), *level),
                )?,
            );
        }
        for (j, f) in faces.iter().enumerate() {
            let uses = f.loops.iter().any(|(l, _)| match l {
                Loop::Edges(list) => list.iter().any(|(e, _)| crossing.contains(e)),
                _ => false,
            });
            if uses {
                meshes[j] = None;
                next = next.min(j);
            }
        }
    }
    Ok((polys, meshes.into_iter().flatten().collect()))
}

/// Vertices joined by a zero-length edge (a sliver face's), each mapped to the smallest id of
/// its group: they are one point, so one mesh vertex.
fn joined_vertices(b: &Brep) -> HashMap<u64, u64> {
    let mut parent: HashMap<u64, u64> = HashMap::new();
    let root = |parent: &HashMap<u64, u64>, mut v: u64| {
        while let Some(&p) = parent.get(&v) {
            v = p;
        }
        v
    };
    for e in b.edges.values() {
        let zero = norm(sub(e.start, e.end)) <= 1e-9 * (1.0 + norm(e.start));
        if zero && e.start_vertex != e.end_vertex {
            let (x, y) = (root(&parent, e.start_vertex), root(&parent, e.end_vertex));
            if x != y {
                parent.insert(x.max(y), x.min(y));
            }
        }
    }
    let ids: Vec<u64> = parent.keys().copied().collect();
    ids.into_iter().map(|v| (v, root(&parent, v))).collect()
}

pub(crate) fn triangle_area(a: P3, b: P3, c: P3) -> f64 {
    let n = cross(
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    0.5 * dot(n, n).sqrt()
}
