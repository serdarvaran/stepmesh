//! Triangulation of one face in its (u, v) parameter space.

use std::collections::{HashMap, HashSet, VecDeque};

use spade::handles::FixedVertexHandle;
use spade::{ConstrainedDelaunayTriangulation, Point2, PositionInTriangulation, Triangulation};

use super::edges::sample_curve;
use crate::geom::Surface;
use crate::model::{Face, Loop, P3, add, dot, mul, norm, sub};
use crate::{Error, Options};

/// Refinement passes over the face; each inserts the centroid or edge midpoints of every
/// triangle that is still too far from the surface.
const MAX_REFINE_ROUNDS: usize = 40;
/// A boundary point may lie this far from the face's surface (file tolerance).
const ON_SURFACE: f64 = 0.05;

/// A mesh vertex: a sample of an edge polyline (shared with the neighbouring face) or a point
/// of this face only.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VRef {
    Edge { edge: u64, i: usize },
    Inner(usize),
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct FaceMesh {
    pub inner: Vec<P3>,
    /// Counter-clockwise seen from the side the face's normal points to.
    pub tris: Vec<[VRef; 3]>,
}

#[derive(Debug, Clone, Copy)]
struct BPoint {
    uv: (f64, f64),
    p: P3,
    r: VRef,
    /// The surface collapses here (pole, apex): copies of this point at other parameters are
    /// the same vertex.
    pole: bool,
}

struct Ctx<'a> {
    f: &'a Face,
    s: &'a Surface,
    o: &'a Options,
    inner: Vec<P3>,
}

pub fn tessellate_face(
    f: &Face,
    edges: &HashMap<u64, Vec<P3>>,
    o: &Options,
) -> Result<FaceMesh, Error> {
    let mut c = Ctx {
        f,
        s: &f.surface,
        o,
        inner: Vec::new(),
    };
    let fail = Error::Projection { face: f.id };
    let mut loops = Vec::new();
    for (l, orientation) in &f.loops {
        if let Loop::Edges(list) = l {
            let pts = loop_points(list, *orientation, edges).ok_or_else(|| fail.clone())?;
            loops.push(c.map_loop(&pts)?);
        }
    }
    if loops.is_empty() {
        loops = c.full_domain()?;
    } else if loops.iter().any(|l| c.wraps(l)) {
        // A loop around a periodic surface without a seam edge: not handled yet.
        return Err(fail);
    }
    c.align(&mut loops);
    c.triangulate(&loops)
}

/// Points of a bound in traversal order, each edge's last point left to the next edge.
fn loop_points(
    list: &[(u64, bool)],
    orientation: bool,
    edges: &HashMap<u64, Vec<P3>>,
) -> Option<Vec<(VRef, P3)>> {
    let order: Vec<(u64, bool)> = if orientation {
        list.to_vec()
    } else {
        list.iter().rev().map(|&(e, o)| (e, !o)).collect()
    };
    let mut out = Vec::new();
    for (e, forward) in order {
        let poly = edges.get(&e)?;
        let n = poly.len();
        let idx: Vec<usize> = if forward {
            (0..n - 1).collect()
        } else {
            (1..n).rev().collect()
        };
        out.extend(
            idx.into_iter()
                .map(|i| (VRef::Edge { edge: e, i }, poly[i])),
        );
    }
    (!out.is_empty()).then_some(out)
}

impl Ctx<'_> {
    fn fail(&self) -> Error {
        Error::Projection { face: self.f.id }
    }

    /// Whether the surface collapses along u (or v) at this parameter point.
    fn singular(&self, uv: (f64, f64)) -> (bool, bool) {
        let (p, su, sv) = self.s.derivs(uv.0, uv.1);
        let scale = 1e-9 * (1.0 + norm(p)).max(norm(su).max(norm(sv)));
        (norm(su) <= scale, norm(sv) <= scale)
    }

    /// Parameters of a loop, continuous across periodic seams. A point where the surface
    /// collapses (pole, apex) is spread over the parameter range between its neighbours, in
    /// steps of at most the angular deflection, so the triangles around it form a fan.
    fn map_loop(&self, pts: &[(VRef, P3)]) -> Result<Vec<BPoint>, Error> {
        let mut uv: Vec<Option<(f64, f64)>> = vec![None; pts.len()];
        let mut sing = vec![(false, false); pts.len()];
        let mut hint = None;
        for (i, &(_, p)) in pts.iter().enumerate() {
            let q = self.s.project(p, hint).ok_or_else(|| self.fail())?;
            if norm(sub(self.s.point(q.0, q.1), p)) > ON_SURFACE {
                return Err(self.fail());
            }
            sing[i] = self.singular(q);
            uv[i] = Some(q);
            if !sing[i].0 && !sing[i].1 {
                hint = Some(q);
            }
        }
        let n = pts.len();
        let regular = |i: usize| !sing[i].0 && !sing[i].1;
        if !(0..n).any(regular) {
            return Err(self.fail());
        }
        let neighbour = |i: usize, step: isize| {
            let mut j = i;
            loop {
                j = (j as isize + step).rem_euclid(n as isize) as usize;
                if regular(j) {
                    return uv[j].unwrap_or((0.0, 0.0));
                }
            }
        };
        let mut out = Vec::with_capacity(n + 2);
        for i in 0..n {
            let (r, p) = pts[i];
            let q = uv[i].unwrap_or((0.0, 0.0));
            if regular(i) {
                out.push(BPoint {
                    uv: q,
                    p,
                    r,
                    pole: false,
                });
                continue;
            }
            let (a, b) = (neighbour(i, -1), neighbour(i, 1));
            let (from, to) = if sing[i].0 { (a.0, b.0) } else { (a.1, b.1) };
            let steps = ((to - from).abs() / self.o.angular_deflection.min(0.5))
                .ceil()
                .max(1.0) as usize;
            for k in 0..=steps {
                let x = from + (to - from) * k as f64 / steps as f64;
                let uv = if sing[i].0 { (x, q.1) } else { (q.0, x) };
                out.push(BPoint {
                    uv,
                    p,
                    r,
                    pole: true,
                });
            }
        }
        Ok(out)
    }

    /// Whether a mapped loop ends a period away from where it starts.
    fn wraps(&self, l: &[BPoint]) -> bool {
        let (Some(first), Some(last)) = (l.first(), l.last()) else {
            return false;
        };
        let (pu, pv) = self.s.period();
        let gap = |a: f64, b: f64, p: Option<f64>| p.is_some_and(|p| (a - b).abs() > 0.5 * p);
        gap(first.uv.0, last.uv.0, pu) || gap(first.uv.1, last.uv.1, pv)
    }

    /// Shifts every loop by whole periods next to the one with the widest parameter range.
    fn align(&self, loops: &mut [Vec<BPoint>]) {
        let (pu, pv) = self.s.period();
        let centre = |l: &[BPoint]| {
            let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
            for b in l {
                lo = (lo.0.min(b.uv.0), lo.1.min(b.uv.1));
                hi = (hi.0.max(b.uv.0), hi.1.max(b.uv.1));
            }
            (
                ((lo.0 + hi.0) / 2.0, (lo.1 + hi.1) / 2.0),
                (hi.0 - lo.0) + (hi.1 - lo.1),
            )
        };
        let Some(reference) =
            (0..loops.len()).max_by(|&a, &b| centre(&loops[a]).1.total_cmp(&centre(&loops[b]).1))
        else {
            return;
        };
        let (rc, _) = centre(&loops[reference]);
        for l in loops.iter_mut() {
            let (c, _) = centre(l);
            let du = pu.map_or(0.0, |p| ((rc.0 - c.0) / p).round() * p);
            let dv = pv.map_or(0.0, |p| ((rc.1 - c.1) / p).round() * p);
            for b in l.iter_mut() {
                b.uv = (b.uv.0 + du, b.uv.1 + dv);
            }
        }
    }

    fn new_inner(&mut self, p: P3) -> VRef {
        self.inner.push(p);
        VRef::Inner(self.inner.len() - 1)
    }

    /// Boundary of a face that covers its whole closed surface (a sphere bounded only by a
    /// vertex loop): the parameter rectangle, with opposite periodic sides sharing their points
    /// and collapsed sides reduced to one point.
    fn full_domain(&mut self) -> Result<Vec<Vec<BPoint>>, Error> {
        let (pu, pv) = self.s.period();
        let [(_, _), (vb0, vb1)] = self.s.param_bounds();
        let Some(pu) = pu else {
            return Err(self.fail());
        };
        let (u0, u1) = (0.0, pu);
        let (v0, v1) = match pv {
            Some(p) => (0.0, p),
            None if vb0.is_finite() && vb1.is_finite() => (vb0, vb1),
            None => return Err(self.fail()),
        };
        let s = self.s.clone();
        let o = self.o.clone();
        // Side samples (parameters) along u at v = v0 and along v at u = u0.
        let bottom_collapsed = self.singular((u0, v0)).0;
        let top_collapsed = self.singular((u0, v1)).0;
        let us = if bottom_collapsed {
            let n = (pu / o.angular_deflection.min(0.5)).ceil() as usize;
            (0..=n)
                .map(|i| u0 + (u1 - u0) * i as f64 / n as f64)
                .collect()
        } else {
            sample_curve(
                &|u| {
                    let (p, su, _) = s.derivs(u, v0);
                    (p, su)
                },
                u0,
                u1,
                &o,
            )
        };
        let vs = sample_curve(
            &|v| {
                let (p, _, sv) = s.derivs(u0, v);
                (p, sv)
            },
            v0,
            v1,
            &o,
        );
        // Left side (u = u0) shared with the right side (u = u1).
        let side_v: Vec<VRef> = vs.iter().map(|&v| self.new_inner(s.point(u0, v))).collect();
        let last = us.len() - 1;
        let mut row = |v: f64, collapsed: bool, end: VRef| -> Vec<VRef> {
            let mut r: Vec<VRef> = if collapsed {
                vec![end; us.len()]
            } else {
                us.iter().map(|&u| self.new_inner(s.point(u, v))).collect()
            };
            r[0] = end;
            r[last] = end;
            r
        };
        let bottom = row(v0, bottom_collapsed, side_v[0]);
        let top = if pv.is_some() {
            bottom.clone()
        } else {
            row(v1, top_collapsed, side_v[vs.len() - 1])
        };
        let pos = |r: VRef, this: &Self| match r {
            VRef::Inner(i) => this.inner[i],
            VRef::Edge { .. } => [0.0; 3],
        };
        let mut l = Vec::new();
        let mut push = |r: VRef, uv: (f64, f64), pole: bool, this: &Self| {
            l.push(BPoint {
                uv,
                p: pos(r, this),
                r,
                pole,
            })
        };
        for (i, &u) in us.iter().enumerate().take(last) {
            push(bottom[i], (u, v0), bottom_collapsed, self);
        }
        for (j, &v) in vs.iter().enumerate().take(vs.len() - 1) {
            push(side_v[j], (u1, v), j == 0 && bottom_collapsed, self);
        }
        for i in (1..=last).rev() {
            push(top[i], (us[i], v1), top_collapsed, self);
        }
        for j in (1..vs.len()).rev() {
            push(
                side_v[j],
                (u0, vs[j]),
                j == vs.len() - 1 && top_collapsed,
                self,
            );
        }
        Ok(vec![l])
    }

    fn triangulate(&mut self, loops: &[Vec<BPoint>]) -> Result<FaceMesh, Error> {
        let (su, sv) = self.metric(loops);
        let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
            ConstrainedDelaunayTriangulation::new();
        // Per vertex handle: its reference, parameters and position.
        let mut info: Vec<BPoint> = Vec::new();
        let face = self.f.id;
        let insert = |cdt: &mut ConstrainedDelaunayTriangulation<Point2<f64>>,
                      info: &mut Vec<BPoint>,
                      b: BPoint|
         -> Result<FixedVertexHandle, Error> {
            let h = cdt
                .insert(Point2::new(b.uv.0 * su, b.uv.1 * sv))
                .map_err(|_| Error::Projection { face })?;
            if h.index() == info.len() {
                info.push(b);
            } else if info[h.index()].r != b.r && norm(sub(info[h.index()].p, b.p)) > ON_SURFACE {
                return Err(Error::Projection { face });
            }
            Ok(h)
        };
        for l in loops {
            let handles = l
                .iter()
                .map(|&b| insert(&mut cdt, &mut info, b))
                .collect::<Result<Vec<_>, _>>()?;
            for i in 0..handles.len() {
                let (a, b) = (handles[i], handles[(i + 1) % handles.len()]);
                if a != b && !cdt.can_add_constraint(a, b) && !cdt.exists_constraint(a, b) {
                    return Err(self.fail());
                }
                if a != b && !cdt.exists_constraint(a, b) {
                    cdt.add_constraint(a, b);
                }
            }
        }
        let planar = matches!(self.s, Surface::Plane(_));
        for _ in 0..if planar { 0 } else { MAX_REFINE_ROUNDS } {
            let inside = classify(&cdt);
            let candidates = self.refinement_points(&cdt, &inside, &info);
            if candidates.is_empty() {
                break;
            }
            for uv in candidates {
                // A point on a boundary edge would split it in this face only: a T-junction.
                match cdt.locate(Point2::new(uv.0 * su, uv.1 * sv)) {
                    PositionInTriangulation::OnVertex(_) => continue,
                    PositionInTriangulation::OnEdge(e)
                        if cdt.is_constraint_edge(e.as_undirected()) =>
                    {
                        continue;
                    }
                    _ => {}
                }
                let p = self.s.point(uv.0, uv.1);
                let r = self.new_inner(p);
                insert(
                    &mut cdt,
                    &mut info,
                    BPoint {
                        uv,
                        p,
                        r,
                        pole: false,
                    },
                )?;
            }
        }
        let inside = classify(&cdt);
        let flip = !self.f.same_sense;
        let mut tris = Vec::new();
        for face in cdt.inner_faces() {
            if !inside.contains(&face.fix().index()) {
                continue;
            }
            let [a, b, c] = face.vertices().map(|v| info[v.fix().index()].r);
            if a == b || b == c || a == c {
                continue; // collapsed onto a pole or apex
            }
            tris.push(if flip { [a, c, b] } else { [a, b, c] });
        }
        Ok(FaceMesh {
            inner: std::mem::take(&mut self.inner),
            tris,
        })
    }

    /// Scale of u and v so that parameter distances approximate surface distances.
    fn metric(&self, loops: &[Vec<BPoint>]) -> (f64, f64) {
        let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
        for b in loops.iter().flatten() {
            lo = (lo.0.min(b.uv.0), lo.1.min(b.uv.1));
            hi = (hi.0.max(b.uv.0), hi.1.max(b.uv.1));
        }
        let (mut su, mut sv, mut n) = (0.0, 0.0, 0.0);
        for i in 0..5 {
            for j in 0..5 {
                let (u, v) = (
                    lo.0 + (hi.0 - lo.0) * i as f64 / 4.0,
                    lo.1 + (hi.1 - lo.1) * j as f64 / 4.0,
                );
                let (_, du, dv) = self.s.derivs(u, v);
                su += norm(du);
                sv += norm(dv);
                n += 1.0;
            }
        }
        let fix = |x: f64| if x > 0.0 && x.is_finite() { x / n } else { 1.0 };
        (fix(su), fix(sv))
    }

    /// Centroids and edge midpoints of inside triangles that are too far from the surface or
    /// that tilt too far from the surface.
    fn refinement_points(
        &self,
        cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>,
        inside: &HashSet<usize>,
        info: &[BPoint],
    ) -> Vec<(f64, f64)> {
        let (d, a) = (self.o.linear_deflection, self.o.angular_deflection);
        let mut out = Vec::new();
        let mut seen_edges = HashSet::new();
        let cos_a = a.cos();
        for face in cdt.inner_faces() {
            if !inside.contains(&face.fix().index()) {
                continue;
            }
            // A triangle with two copies of a pole is a sliver of the fan: nothing to refine. Two
            // copies of a seam point, though, mean the triangle spans the whole face.
            let v = face.vertices().map(|h| info[h.fix().index()]);
            if (0..3).any(|i| v[i].pole && v[i].r == v[(i + 1) % 3].r) {
                continue;
            }
            let uv_c = (
                (v[0].uv.0 + v[1].uv.0 + v[2].uv.0) / 3.0,
                (v[0].uv.1 + v[1].uv.1 + v[2].uv.1) / 3.0,
            );
            let p_c = mul(add(add(v[0].p, v[1].p), v[2].p), 1.0 / 3.0);
            let off = norm(sub(self.s.point(uv_c.0, uv_c.1), p_c));
            // The surface turns too much between the corners (pole copies excluded: there the
            // normal depends on the side the pole is approached from).
            let normals: Vec<P3> = v
                .iter()
                .filter(|x| !x.pole)
                .map(|x| self.s.normal(x.uv.0, x.uv.1))
                .collect();
            let turns = (0..normals.len())
                .any(|i| (i + 1..normals.len()).any(|j| dot(normals[i], normals[j]) < cos_a));
            if off > d || turns {
                out.push(uv_c);
            }
            for e in face.adjacent_edges() {
                let undirected = e.fix().as_undirected();
                if cdt.is_constraint_edge(undirected) || !seen_edges.insert(undirected.index()) {
                    continue;
                }
                let [x, y] = e.vertices().map(|h| info[h.fix().index()]);
                let uv_m = ((x.uv.0 + y.uv.0) / 2.0, (x.uv.1 + y.uv.1) / 2.0);
                if norm(sub(self.s.point(uv_m.0, uv_m.1), mul(add(x.p, y.p), 0.5))) > d {
                    out.push(uv_m);
                }
            }
        }
        out
    }
}

/// Inner faces inside the boundary: an odd number of constraint edges away from the outside.
fn classify(cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>) -> HashSet<usize> {
    let mut parity: HashMap<usize, bool> = HashMap::new();
    let mut queue = VecDeque::new();
    for face in cdt.inner_faces() {
        for e in face.adjacent_edges() {
            if e.rev().face().is_outer() {
                let odd = cdt.is_constraint_edge(e.fix().as_undirected());
                if let std::collections::hash_map::Entry::Vacant(slot) =
                    parity.entry(face.fix().index())
                {
                    slot.insert(odd);
                    queue.push_back(face.fix());
                }
            }
        }
    }
    while let Some(f) = queue.pop_front() {
        let odd = parity[&f.index()];
        for e in cdt.face(f).adjacent_edges() {
            let Some(g) = e.rev().face().as_inner() else {
                continue;
            };
            if parity.contains_key(&g.fix().index()) {
                continue;
            }
            parity.insert(
                g.fix().index(),
                odd != cdt.is_constraint_edge(e.fix().as_undirected()),
            );
            queue.push_back(g.fix());
        }
    }
    parity
        .into_iter()
        .filter(|&(_, odd)| odd)
        .map(|(i, _)| i)
        .collect()
}
