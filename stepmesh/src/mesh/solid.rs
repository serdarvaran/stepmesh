//! One B-rep solid as a closed, oriented triangle mesh.

use std::collections::HashMap;

use super::edges::discretize_edge;
use super::face::{VRef, tessellate_face};
use crate::model::{Brep, P3, cross, dot};
use crate::{Error, Options};

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
    let mut polys = HashMap::with_capacity(b.edges.len());
    for (&id, e) in &b.edges {
        polys.insert(id, discretize_edge(e, o)?);
    }
    let mut index: HashMap<Key, u32> = HashMap::new();
    let mut positions: Vec<P3> = Vec::new();
    let mut triangles = Vec::new();
    let mut face_no = 0;
    for shell in &b.shells {
        let mut shell_tris = Vec::new();
        for f in &shell.faces {
            let fm = tessellate_face(f, &polys, o)?;
            let mut id = |r: VRef| -> u32 {
                let (key, p) = match r {
                    VRef::Edge { edge, i } => {
                        let e = &b.edges[&edge];
                        let last = polys[&edge].len() - 1;
                        let key = match i {
                            0 => Key::Vertex(e.start_vertex),
                            i if i == last => Key::Vertex(e.end_vertex),
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

pub(crate) fn triangle_area(a: P3, b: P3, c: P3) -> f64 {
    let n = cross(
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    0.5 * dot(n, n).sqrt()
}
