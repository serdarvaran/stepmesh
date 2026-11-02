//! Triangulation of one face in its (u, v) parameter space.

use std::collections::{HashMap, VecDeque};

use spade::handles::{FixedFaceHandle, FixedVertexHandle, InnerTag};
use spade::{ConstrainedDelaunayTriangulation, Point2, PositionInTriangulation, Triangulation};

use super::edges::sample_curve;
use super::param::Param;
use super::seam::{free_param, nearest_along_u, restart_at, seam_pair, seam_route, seam_start};
use crate::fxhash::FxSet;
use crate::geom::Surface;
use crate::model::{Face, Loop, P3, add, dot, mul, norm, sub};
use crate::{Error, Options};

/// Refinement passes over the face; each inserts the centroid or edge midpoints of every
/// triangle that is still too far from the surface, and the next looks only at the triangles
/// around them. Along a long thin face the bad triangles move on a few at a time, so passes are
/// many but cheap; a face still refining after this many is an error.
const MAX_REFINE_ROUNDS: usize = 10_000;
/// Twice the parameter area, as a share of the squared longest side, below which a triangle
/// counts as flat.
const FLAT: f64 = 1e-9;
/// Distance from a boundary edge, as a share of its length, within which a refinement point
/// counts as on it.
const ON_EDGE: f64 = 1e-4;
/// Triangles shorter than this many deflections are left to the distance test alone.
const CREASE_SIZE: f64 = 10.0;
/// Ratio of the parameter derivatives below which a direction counts as collapsed.
const COLLAPSED: f64 = 1e-5;
/// Change of the periodic parameter between neighbouring boundary points, as a share of the
/// period, that only passing a pole explains (edges are sampled in far smaller angular steps).
const POLE_JUMP: f64 = 0.25;
/// Most copies of one pole along its parameter row.
const MAX_POLE_COPIES: usize = 64;
/// A boundary point may lie this far from the face's surface: files carry edge curves that
/// leave their faces by their modelling tolerance (0.08 mm seen in CAx-IF test data).
const ON_SURFACE: f64 = 0.1;

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
pub(super) struct BPoint {
    pub(super) uv: (f64, f64),
    pub(super) p: P3,
    pub(super) r: VRef,
    /// The surface collapses here (pole, apex): copies of this point at other parameters are
    /// the same vertex.
    pub(super) pole: bool,
}

struct Ctx<'a> {
    f: &'a Face,
    s: Param<'a>,
    /// Whether the face's normal is that of `s` (the file's same_sense, reversed when swapped).
    same_sense: bool,
    o: &'a Options,
    inner: Vec<P3>,
    /// Edges whose polylines cross in parameter space, when they do.
    crossing: Vec<u64>,
}

pub fn tessellate_face(
    f: &Face,
    edges: &HashMap<u64, Vec<P3>>,
    o: &Options,
) -> Result<FaceMesh, Error> {
    tessellate(f, edges, o).map_err(|(e, _)| e)
}

/// As [`tessellate_face`]; when the boundary crosses itself, the error comes with the edges whose
/// polylines cross, which finer samples may part.
pub(super) fn tessellate(
    f: &Face,
    edges: &HashMap<u64, Vec<P3>>,
    o: &Options,
) -> Result<FaceMesh, (Error, Vec<u64>)> {
    let plain = |e: Error| (e, Vec::new());
    let ctx = |swapped: bool| Ctx {
        f,
        s: Param {
            surface: &f.surface,
            swapped,
        },
        same_sense: f.same_sense != swapped,
        o,
        inner: Vec::new(),
        crossing: Vec::new(),
    };
    let mut c = ctx(false);
    let mut loops = c.map_loops(edges).map_err(plain)?;
    // Loops going round along v only are closed with u and v swapped.
    let along = |c: &Ctx, l: &[BPoint]| c.advance(l);
    if loops.iter().any(|l| along(&c, l).1 != 0.0) && loops.iter().all(|l| along(&c, l).0 == 0.0) {
        c = ctx(true);
        loops = c.map_loops(edges).map_err(plain)?;
    }
    if loops.is_empty() {
        loops = c.full_domain((0.0, 0.0)).map_err(plain)?;
    } else {
        loops = c.close_wraps(loops).map_err(plain)?;
        // A face can lie outside all its loops on a closed surface (a sphere with holes): the
        // whole parameter domain is then its outer boundary. Not when a loop is given as the
        // outer bound: some exporters write such faces against their surface's normal, and the
        // face is the region inside it.
        let side = if c.same_sense { 1.0 } else { -1.0 };
        // The face lies left of its outer bound seen from its normal. A lone outer bound running
        // the other way was written against the surface's normal (vertex blends on surfaces of
        // revolution): its edges, shared with the neighbours, fix the side the mesh faces.
        if let ([l], Some(_)) = (loops.as_slice(), f.outer) {
            if side * signed_area(l) < 0.0 {
                c.same_sense = !c.same_sense;
            }
        }
        let (pu, pv) = c.s.period();
        if let (Some(pu), None) = (pu, f.outer) {
            if loops.iter().all(|l| side * signed_area(l) < 0.0) {
                // Its seams go where no loop reaches, so that the loops lie inside it.
                let start = (
                    free_param(&loops, 0, pu),
                    pv.map_or(0.0, |p| free_param(&loops, 1, p)),
                );
                loops.extend(c.full_domain(start).map_err(plain)?);
            }
        }
    }
    c.align(&mut loops);
    c.triangulate(&loops)
        .map_err(|e| (e, std::mem::take(&mut c.crossing)))
}

/// Mesh of a plane face bounded by one loop whose polylines cross (two curves closer than the
/// tolerance, a cusp): the loop cut into ears. Any triangulation of the loop that keeps its
/// orientation covers the same signed area of the plane, so the solid's volume stays exact where
/// the boundary overlaps itself.
pub(super) fn clip_ears(
    f: &Face,
    edges: &HashMap<u64, Vec<P3>>,
    o: &Options,
) -> Result<FaceMesh, Error> {
    let c = Ctx {
        f,
        s: Param {
            surface: &f.surface,
            swapped: false,
        },
        same_sense: f.same_sense,
        o,
        inner: Vec::new(),
        crossing: Vec::new(),
    };
    let crossing = || c.fail("the boundary crosses itself in parameter space");
    let loops = c.map_loops(edges)?;
    let [l] = loops.as_slice() else {
        return Err(crossing());
    };
    let ccw = signed_area(l) > 0.0;
    let turn = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    let mut ring: Vec<usize> = (0..l.len()).collect();
    let mut tris = Vec::with_capacity(l.len());
    // The search for the next ear starts where the last one was cut.
    let mut last = 0;
    while ring.len() > 3 {
        let n = ring.len();
        let corner = |k: usize| {
            let (a, b, c) = (ring[(k + n - 1) % n], ring[k], ring[(k + 1) % n]);
            (a, b, c, turn(l[a].uv, l[b].uv, l[c].uv))
        };
        let convex = |t: f64| if ccw { t > 0.0 } else { t < 0.0 };
        let empty = |a: usize, b: usize, c: usize| {
            ring.iter().all(|&j| {
                j == a || j == b || j == c || {
                    let p = l[j].uv;
                    let s = [
                        turn(l[a].uv, l[b].uv, p),
                        turn(l[b].uv, l[c].uv, p),
                        turn(l[c].uv, l[a].uv, p),
                    ];
                    !(s.iter().all(|&x| x >= 0.0) || s.iter().all(|&x| x <= 0.0))
                }
            })
        };
        // An ear empty of other boundary points, else (where the boundary overlaps itself) the
        // smallest convex corner, else the flattest one.
        let k = (0..n)
            .map(|i| (last + i) % n)
            .find(|&k| {
                let (a, b, c, t) = corner(k);
                convex(t) && empty(a, b, c)
            })
            .or_else(|| {
                (0..n)
                    .filter(|&k| convex(corner(k).3))
                    .min_by(|&x, &y| corner(x).3.abs().total_cmp(&corner(y).3.abs()))
            })
            .unwrap_or_else(|| {
                (0..n)
                    .min_by(|&x, &y| corner(x).3.abs().total_cmp(&corner(y).3.abs()))
                    .unwrap_or(0)
            });
        let (a, b, c, _) = corner(k);
        tris.push([l[a].r, l[b].r, l[c].r]);
        ring.remove(k);
        last = k.saturating_sub(1);
    }
    if let [a, b, c] = ring[..] {
        tris.push([l[a].r, l[b].r, l[c].r]);
    }
    // Oriented as `triangulate` orients: counter-clockwise in parameter space seen from the
    // surface's normal, reversed for a face against it.
    let flip = ccw != c.same_sense;
    Ok(FaceMesh {
        inner: Vec::new(),
        tris: tris
            .into_iter()
            .filter(|[a, b, c]| a != b && b != c && a != c)
            .map(|[a, b, c]| if flip { [a, c, b] } else { [a, b, c] })
            .collect(),
    })
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
    fn error(face: u64, reason: &'static str) -> Error {
        Error::Projection { face, reason }
    }

    fn fail(&self, reason: &'static str) -> Error {
        Ctx::error(self.f.id, reason)
    }

    /// Parameters of every edge loop of the face (vertex loops add nothing).
    fn map_loops(&self, edges: &HashMap<u64, Vec<P3>>) -> Result<Vec<Vec<BPoint>>, Error> {
        let f = self.f;
        let edge_loops = f
            .loops
            .iter()
            .filter(|(l, _)| matches!(l, Loop::Edges(_)))
            .count();
        let mut loops = Vec::new();
        for (index, (l, orientation)) in f.loops.iter().enumerate() {
            if let Loop::Edges(list) = l {
                let pts = loop_points(list, *orientation, edges)
                    .ok_or_else(|| self.fail("a loop edge has no polyline"))?;
                let outer = match f.outer {
                    Some(o) => Some(o == index),
                    None => (edge_loops == 1).then_some(true),
                };
                loops.push(self.map_loop(&pts, outer)?);
            }
        }
        Ok(loops)
    }

    /// Whether the surface collapses along u (or v) at this parameter point.
    /// A direction counts as collapsed when its derivative is tiny next to the other one: a
    /// boundary point a hair off a pole (file tolerance) is still the pole.
    fn singular(&self, uv: (f64, f64)) -> (bool, bool) {
        let (p, su, sv) = self.s.derivs(uv.0, uv.1);
        let (a, b) = (norm(su), norm(sv));
        let floor = 1e-9 * (1.0 + norm(p));
        (a <= floor.max(COLLAPSED * b), b <= floor.max(COLLAPSED * a))
    }

    /// Parameters of a loop, continuous across periodic seams. A point where the surface
    /// collapses (pole, apex) is spread over the parameter range between its neighbours, in
    /// steps of at most the angular deflection, so the triangles around it form a fan.
    fn map_loop(&self, pts: &[(VRef, P3)], outer: Option<bool>) -> Result<Vec<BPoint>, Error> {
        let mut uv: Vec<(f64, f64)> = vec![(0.0, 0.0); pts.len()];
        let mut sing = vec![(false, false); pts.len()];
        let mut hint = None;
        for (i, &(_, p)) in pts.iter().enumerate() {
            let q = self
                .s
                .project(p, hint)
                .ok_or_else(|| self.fail("projection onto the surface did not converge"))?;
            if norm(sub(self.s.point(q.0, q.1), p)) > ON_SURFACE {
                return Err(self.fail("a boundary point lies off the surface"));
            }
            sing[i] = self.singular(q);
            uv[i] = q;
            if !sing[i].0 && !sing[i].1 {
                hint = Some(q);
            }
        }
        if !sing.iter().any(|s| !s.0 && !s.1) {
            return Err(self.fail("the loop has no regular point"));
        }
        let plain = self.spread(pts, &uv, &sing);
        // After a pole the next meridian's angle is known only up to a period; when the loop does
        // not close, or closes the wrong way round, try the other choices. A loop passing a hair
        // from a pole (a great circle through both) jumps there by half a period between two
        // regular points: the same choice.
        let Some(period) = self.s.period().0 else {
            return Ok(plain);
        };
        let mut segment = vec![0usize; pts.len()];
        let mut blocks = 0;
        for i in 0..pts.len() {
            let jump = i > 0
                && !sing[i].0
                && !sing[i - 1].0
                && (uv[i].0 - uv[i - 1].0).abs() > POLE_JUMP * period;
            if (sing[i].0 && (i == 0 || !sing[i - 1].0)) || jump {
                blocks += 1;
            }
            segment[i] = blocks;
        }
        if blocks == 0 || blocks > 4 {
            return Ok(plain);
        }
        let want = outer.map(|o| if o == self.same_sense { 1.0 } else { -1.0 });
        let fits = |l: &[BPoint]| {
            self.advance(l) == (0.0, 0.0) && want.is_none_or(|w| w * signed_area(l) > 0.0)
        };
        if fits(&plain) {
            return Ok(plain);
        }
        let mut best: Option<(usize, Vec<BPoint>)> = None;
        for combo in 0..3usize.pow(blocks as u32) {
            let ks: Vec<i32> = (0..blocks)
                .map(|b| (combo / 3usize.pow(b as u32) % 3) as i32 - 1)
                .collect();
            let shifted: Vec<(f64, f64)> = (0..pts.len())
                .map(|i| {
                    let k: i32 = ks[..segment[i]].iter().sum();
                    (uv[i].0 + k as f64 * period, uv[i].1)
                })
                .collect();
            let l = self.spread(pts, &shifted, &sing);
            let cost = ks.iter().map(|k| k.unsigned_abs() as usize).sum();
            if fits(&l) && best.as_ref().is_none_or(|(c, _)| cost < *c) {
                best = Some((cost, l));
            }
        }
        Ok(best.map_or(plain, |(_, l)| l))
    }

    /// The loop's points in parameter space; a point where the surface collapses is spread over
    /// the row between its neighbours.
    fn spread(&self, pts: &[(VRef, P3)], uv: &[(f64, f64)], sing: &[(bool, bool)]) -> Vec<BPoint> {
        let n = pts.len();
        let regular = |i: usize| !sing[i].0 && !sing[i].1;
        let neighbour = |i: usize, step: isize| {
            let mut j = i;
            loop {
                j = (j as isize + step).rem_euclid(n as isize) as usize;
                if regular(j) {
                    return uv[j];
                }
            }
        };
        let mut out = Vec::with_capacity(n + 2);
        for i in 0..n {
            let (r, p) = pts[i];
            let q = uv[i];
            if regular(i) {
                out.push(BPoint {
                    uv: q,
                    p,
                    r,
                    pole: false,
                });
                continue;
            }
            // Points a hair from a pole count as the pole too: a run of them shares the row
            // between its regular neighbours, each taking its own part in order.
            let in_run = |step: isize| {
                (1..n)
                    .map(|k| (i as isize + step * k as isize).rem_euclid(n as isize) as usize)
                    .take_while(|&j| !regular(j) && sing[j] == sing[i])
                    .count()
            };
            let (before, run) = (in_run(-1), in_run(-1) + in_run(1) + 1);
            let (a, b) = (neighbour(i, -1), neighbour(i, 1));
            let (start, end) = if sing[i].0 { (a.0, b.0) } else { (a.1, b.1) };
            let part = (end - start) / run as f64;
            let (from, to) = (
                start + part * before as f64,
                start + part * (before + 1) as f64,
            );
            let steps = self.pole_steps(to - from);
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
        out
    }

    /// How far a mapped loop advances along a periodic direction: a whole period (signed) when
    /// it goes round the surface without crossing a seam, `(0, 0)` otherwise.
    fn advance(&self, l: &[BPoint]) -> (f64, f64) {
        let (Some(first), Some(last)) = (l.first(), l.last()) else {
            return (0.0, 0.0);
        };
        let (pu, pv) = self.s.period();
        let gap = |a: f64, b: f64, p: Option<f64>| match p {
            Some(p) if (b - a).abs() > 0.5 * p => ((b - a) / p).round() * p,
            _ => 0.0,
        };
        (
            gap(first.uv.0, last.uv.0, pu),
            gap(first.uv.1, last.uv.1, pv),
        )
    }

    /// Loops that go round a periodic surface without a seam edge become one parameter polygon:
    /// a single one is closed through the pole on the side the face lies (spherical cap), two
    /// running opposite ways are joined into a band (cylinder without a seam).
    fn close_wraps(&mut self, loops: Vec<Vec<BPoint>>) -> Result<Vec<Vec<BPoint>>, Error> {
        let mut plain = Vec::new();
        let mut wrapping = Vec::new();
        for l in loops {
            match self.advance(&l) {
                (0.0, 0.0) => plain.push(l),
                (du, 0.0) => wrapping.push((l, du)),
                _ => return Err(self.fail("a loop goes round the surface along v")),
            }
        }
        match wrapping.len() {
            0 => {}
            1 => {
                let (l, du) = wrapping.pop().ok_or_else(|| self.fail("no loop"))?;
                let closed = self.close_at_pole(l, du, &plain)?;
                plain.push(closed);
            }
            2 => {
                let (b, db) = wrapping.pop().ok_or_else(|| self.fail("no loop"))?;
                let (a, da) = wrapping.pop().ok_or_else(|| self.fail("no loop"))?;
                if da.signum() == db.signum() {
                    return Err(self.fail("two loops go round the surface the same way"));
                }
                let band = self.join_band(a, da, b, db, &plain);
                plain.push(band);
            }
            _ => return Err(self.fail("more than two loops go round the surface")),
        }
        Ok(plain)
    }

    /// The loop, then its first point a period on, then the pole back to the start. The face
    /// lies left of the loop in parameter space (right when it runs against its surface), which
    /// picks the pole.
    fn close_at_pole(
        &mut self,
        l: Vec<BPoint>,
        du: f64,
        holes: &[Vec<BPoint>],
    ) -> Result<Vec<BPoint>, Error> {
        let [_, (v0, v1)] = self.s.param_bounds();
        let north = (du > 0.0) == self.same_sense;
        let v = if north { v1 } else { v0 };
        let mut l = restart_at(&l, du, seam_start(&l, v, holes, du.abs()));
        let first = l[0];
        if !v.is_finite() || !self.singular((first.uv.0, v)).0 {
            return Err(self.fail("a loop goes round the surface and no pole closes it"));
        }
        let end = first.uv.0 + du;
        let pole_p = self.s.point(first.uv.0, v);
        let pole_r = self.new_inner(pole_p);
        // The meridian from the loop to the pole is a seam inside the face: sampled once, its
        // points stand on both sides.
        let seam = self.seam_points(first.uv, (first.uv.0, v));
        l.push(BPoint {
            uv: (end, first.uv.1),
            ..first
        });
        l.extend(seam.iter().map(|b| BPoint {
            uv: (b.uv.0 + du, b.uv.1),
            ..*b
        }));
        let steps = self.pole_steps(du);
        for k in 0..=steps {
            let u = end - du * k as f64 / steps as f64;
            l.push(BPoint {
                uv: (u, v),
                p: pole_p,
                r: pole_r,
                pole: true,
            });
        }
        l.extend(seam.iter().rev().copied());
        Ok(l)
    }

    /// Two loops going round a periodic surface in opposite directions, joined at their nearest
    /// points into one polygon: all of `a`, its first point a period on, a seam to `b`, `b` from
    /// there round to the same point a period back, and the seam back. Both copies of the seam
    /// share their points. On a surface closed along v too (torus) the loops bound two bands:
    /// the face lies left of `a` in parameter space (right when it runs against its surface),
    /// which places `b` above or below it.
    fn join_band(
        &mut self,
        a: Vec<BPoint>,
        da: f64,
        mut b: Vec<BPoint>,
        db: f64,
        holes: &[Vec<BPoint>],
    ) -> Vec<BPoint> {
        if let Some(pv) = self.s.period().1 {
            let mean_v = |l: &[BPoint]| l.iter().map(|p| p.uv.1).sum::<f64>() / l.len() as f64;
            let gap = (mean_v(&a) - mean_v(&b)) / pv;
            let k = if (da > 0.0) == self.same_sense {
                gap.floor() + 1.0
            } else {
                gap.ceil() - 1.0
            };
            for p in &mut b {
                p.uv.1 += k * pv;
            }
        }
        let period = da.abs();
        // Straight where a straight seam passes, else winding between the holes.
        let (i, k, route) = match seam_pair(&a, &b, holes, period) {
            Some((i, k)) => (i, k, None),
            None => match seam_route(&a, &b, holes, period) {
                Some(r) => (r.from, r.to, Some(r.corners)),
                None => {
                    let mean_v = b.iter().map(|p| p.uv.1).sum::<f64>() / b.len() as f64;
                    let others: Vec<Vec<BPoint>> =
                        holes.iter().cloned().chain([b.clone()]).collect();
                    let i = seam_start(&a, mean_v, &others, period);
                    (i, nearest_along_u(&b, a[i].uv.0, period), None)
                }
            },
        };
        let a = restart_at(&a, da, i);
        // Corners of the seam from `a`'s first point to `b`'s point next to it.
        let route = match route {
            Some(r) => {
                let s = a[0].uv.0 - r[0].0;
                r.iter().map(|&(u, v)| (u + s, v)).collect()
            }
            None => {
                let shift = ((a[0].uv.0 + da - b[k].uv.0) / period).round() * period;
                vec![a[0].uv, (b[k].uv.0 + shift - da, b[k].uv.1)]
            }
        };
        let end = route[route.len() - 1];
        let shift = ((end.0 + da - b[k].uv.0) / period).round() * period;
        let seam = self.route_points(&route);
        let mut out = a.clone();
        out.push(BPoint {
            uv: (a[0].uv.0 + da, a[0].uv.1),
            ..a[0]
        });
        out.extend(seam.iter().map(|s| BPoint {
            uv: (s.uv.0 + da, s.uv.1),
            ..*s
        }));
        for i in 0..=b.len() {
            let j = (k + i) % b.len();
            // Points after the end of `b` continue past its wrap.
            let extra = if k + i >= b.len() { db } else { 0.0 };
            out.push(BPoint {
                uv: (b[j].uv.0 + shift + extra, b[j].uv.1),
                ..b[j]
            });
        }
        out.extend(seam.iter().rev().copied());
        out
    }

    /// Points of a seam along the straight parameter lines between its corners, the inner
    /// corners included (new face vertices).
    fn route_points(&mut self, route: &[(f64, f64)]) -> Vec<BPoint> {
        let mut out = Vec::new();
        for (k, w) in route.windows(2).enumerate() {
            if k > 0 {
                let p = self.s.point(w[0].0, w[0].1);
                let r = self.new_inner(p);
                out.push(BPoint {
                    uv: w[0],
                    p,
                    r,
                    pole: false,
                });
            }
            out.extend(self.seam_points(w[0], w[1]));
        }
        out
    }

    /// Points strictly between two parameter points, along the straight line joining them, at
    /// the edge tolerances (new face vertices).
    fn seam_points(&mut self, from: (f64, f64), to: (f64, f64)) -> Vec<BPoint> {
        let at = |t: f64| (from.0 + (to.0 - from.0) * t, from.1 + (to.1 - from.1) * t);
        let ts = self.segment_params(from, to);
        ts[1..ts.len() - 1]
            .iter()
            .map(|&t| {
                let uv = at(t);
                let p = self.s.point(uv.0, uv.1);
                let r = self.new_inner(p);
                BPoint {
                    uv,
                    p,
                    r,
                    pole: false,
                }
            })
            .collect()
    }

    /// Parameters (0 to 1) along the straight parameter line from `from` to `to` at the edge
    /// tolerances, starting from where it crosses the surface's knot lines.
    fn segment_params(&self, from: (f64, f64), to: (f64, f64)) -> Vec<f64> {
        let s = self.s;
        let (period, breaks) = (s.period(), s.breaks());
        let mut cuts = Vec::new();
        for (a, b, p, knots) in [
            (from.0, to.0, period.0, &breaks[0]),
            (from.1, to.1, period.1, &breaks[1]),
        ] {
            if a == b {
                continue;
            }
            let (lo, hi) = (a.min(b), a.max(b));
            for &k in knots {
                // A periodic parameter meets the knot once per period it runs through.
                let mut x = match p {
                    Some(p) => k + ((lo - k) / p).floor() * p,
                    None => k,
                };
                loop {
                    if x > lo && x < hi {
                        cuts.push((x - a) / (b - a));
                    }
                    match p {
                        Some(p) if x < hi => x += p,
                        _ => break,
                    }
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let curve = |t: f64| {
            let (u, v) = (from.0 + (to.0 - from.0) * t, from.1 + (to.1 - from.1) * t);
            let (p, su, sv) = s.derivs(u, v);
            (p, add(mul(su, to.0 - from.0), mul(sv, to.1 - from.1)))
        };
        sample_curve(&curve, 0.0, 1.0, &cuts, self.o)
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
        let (lo, hi) = bounds(&loops[reference]);
        // Along a direction the reference spans more than a period of (a helical strip), more
        // than one shift puts a loop beside it; the loop belongs where it lies inside it.
        let shifts = |period: Option<f64>, width: f64| -> Vec<f64> {
            match period {
                Some(p) if width > p => {
                    let n = (width / p).ceil() as i32 + 1;
                    let mut ks: Vec<i32> = (-n..=n).collect();
                    ks.sort_by_key(|k| k.abs());
                    ks.into_iter().map(|k| k as f64 * p).collect()
                }
                _ => vec![0.0],
            }
        };
        let (su, sv) = (shifts(pu, hi.0 - lo.0), shifts(pv, hi.1 - lo.1));
        let outer = if su.len() > 1 || sv.len() > 1 {
            loops[reference].clone()
        } else {
            Vec::new()
        };
        for (i, l) in loops.iter_mut().enumerate() {
            let (c, _) = centre(l);
            let du = pu.map_or(0.0, |p| ((rc.0 - c.0) / p).round() * p);
            let dv = pv.map_or(0.0, |p| ((rc.1 - c.1) / p).round() * p);
            let (du, dv) = if i == reference || (su.len() == 1 && sv.len() == 1) {
                (du, dv)
            } else {
                su.iter()
                    .flat_map(|&a| sv.iter().map(move |&b| (du + a, dv + b)))
                    .find(|&(a, b)| mostly_inside(l, (a, b), &outer))
                    .unwrap_or((du, dv))
            };
            for b in l.iter_mut() {
                b.uv = (b.uv.0 + du, b.uv.1 + dv);
            }
        }
    }

    /// Copies of a pole along a parameter span: one per angular step, at most `MAX_POLE_COPIES`
    /// (a B-spline's parameter span need not be an angle).
    fn pole_steps(&self, span: f64) -> usize {
        ((span.abs() / self.o.angular_deflection.min(0.5))
            .ceil()
            .max(1.0) as usize)
            .min(MAX_POLE_COPIES)
    }

    fn new_inner(&mut self, p: P3) -> VRef {
        self.inner.push(p);
        VRef::Inner(self.inner.len() - 1)
    }

    /// Boundary of a face that covers its whole closed surface (a sphere bounded only by a
    /// vertex loop): the parameter rectangle, with opposite periodic sides sharing their points
    /// and collapsed sides reduced to one point.
    /// `start` holds the parameters of its seams (v only on a surface closed along v).
    fn full_domain(&mut self, start: (f64, f64)) -> Result<Vec<Vec<BPoint>>, Error> {
        let (pu, pv) = self.s.period();
        let [(_, _), (vb0, vb1)] = self.s.param_bounds();
        let Some(pu) = pu else {
            return Err(self.fail("the surface is closed but not periodic"));
        };
        let (u0, u1) = (start.0, start.0 + pu);
        let (v0, v1) = match pv {
            Some(p) => (start.1, start.1 + p),
            None if vb0.is_finite() && vb1.is_finite() => (vb0, vb1),
            None => return Err(self.fail("the surface has no finite parameter range")),
        };
        let s = self.s;
        // Side samples (parameters) along u at v = v0 and along v at u = u0.
        let bottom_collapsed = self.singular((u0, v0)).0;
        let top_collapsed = self.singular((u0, v1)).0;
        let us: Vec<f64> = if bottom_collapsed {
            let n = self.pole_steps(pu);
            (0..=n)
                .map(|i| u0 + (u1 - u0) * i as f64 / n as f64)
                .collect()
        } else {
            self.segment_params((u0, v0), (u1, v0))
                .iter()
                .map(|t| u0 + (u1 - u0) * t)
                .collect()
        };
        let vs: Vec<f64> = self
            .segment_params((u0, v0), (u0, v1))
            .iter()
            .map(|t| v0 + (v1 - v0) * t)
            .collect();
        // Left side (u = u0) shared with the right side (u = u1).
        let mut side_v: Vec<VRef> = vs.iter().map(|&v| self.new_inner(s.point(u0, v))).collect();
        if pv.is_some() {
            // Closed along v too: the side ends where it starts.
            let last = side_v.len() - 1;
            side_v[last] = side_v[0];
        }
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
                .insert(cdt_point(b.uv, su, sv))
                .map_err(|_| Ctx::error(face, "a parameter point cannot be triangulated"))?;
            if h.index() == info.len() {
                info.push(b);
            } else if info[h.index()].r != b.r && norm(sub(info[h.index()].p, b.p)) > ON_SURFACE {
                return Err(Ctx::error(face, "two boundary points share parameters"));
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
                    let mut crossing: Vec<u64> = cdt
                        .get_conflicting_edges_between_vertices(a, b)
                        .flat_map(|e| e.vertices().map(|v| v.fix()))
                        .chain([a, b])
                        .filter_map(|v| match info[v.index()].r {
                            VRef::Edge { edge, .. } => Some(edge),
                            VRef::Inner(_) => None,
                        })
                        .collect();
                    crossing.sort_unstable();
                    crossing.dedup();
                    self.crossing = crossing;
                    return Err(self.fail("the boundary crosses itself in parameter space"));
                }
                if a != b && !cdt.exists_constraint(a, b) {
                    cdt.add_constraint(a, b);
                }
            }
        }
        let planar = matches!(self.s.surface, Surface::Plane(_));
        let mut checked = Checked::default();
        let mut converged = planar;
        // After the first round only triangles around the points just inserted are new.
        let mut fresh: Option<Vec<FixedVertexHandle>> = None;
        for _ in 0..if planar { 0 } else { MAX_REFINE_ROUNDS } {
            let inside = classify(&cdt);
            if let Some(limit) = self.o.max_triangles {
                let count = cdt.num_inner_faces();
                if count > limit {
                    return Err(Error::TriangleLimit { count, limit });
                }
            }
            let candidates =
                self.refinement_points(&cdt, &inside, &info, &mut checked, fresh.as_deref());
            if candidates.is_empty() {
                converged = true;
                break;
            }
            let mut inserted = Vec::new();
            for (uv, near) in candidates {
                // A point on a boundary edge would split it in this face only: a T-junction.
                let at = cdt_point(uv, su, sv);
                match cdt.locate_with_hint(at, near) {
                    PositionInTriangulation::OnVertex(_) => continue,
                    PositionInTriangulation::OnEdge(e)
                        if cdt.is_constraint_edge(e.as_undirected()) =>
                    {
                        continue;
                    }
                    // A hair off a boundary edge counts as on it: such points pile up into flat
                    // triangles the triangulation can no longer search.
                    PositionInTriangulation::OnFace(f)
                        if cdt.face(f).adjacent_edges().iter().any(|e| {
                            cdt.is_constraint_edge(e.fix().as_undirected())
                                && hair_off(at, e.from().position(), e.to().position())
                        }) =>
                    {
                        continue;
                    }
                    _ => {}
                }
                let p = self.s.point(uv.0, uv.1);
                let r = self.new_inner(p);
                inserted.push(insert(
                    &mut cdt,
                    &mut info,
                    BPoint {
                        uv,
                        p,
                        r,
                        pole: false,
                    },
                )?);
            }
            // What is left lies on boundary points or chords: the boundary polylines fix it.
            if inserted.is_empty() {
                converged = true;
                break;
            }
            fresh = Some(inserted);
        }
        if !converged {
            return Err(self.fail("the mesh does not reach the tolerances"));
        }
        let inside = classify(&cdt);
        let flip = !self.same_sense;
        let mut tris = Vec::new();
        for face in cdt.inner_faces() {
            if !inside[face.fix().index()] {
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

    /// Whether `p`, a point of a triangle or an edge whose longest side is `size` long, lies
    /// within the deflection of the surface. The surface point at the matching parameters `uv`
    /// is checked first; when that is too far, the nearest surface point is sought from there.
    /// The parameter line between boundary points need not follow the chord between them, so
    /// the first test alone flags slivers along the boundary forever; the nearest point must
    /// stay within `size` of the first one, not land on some other part of the surface.
    fn near_surface(&self, p: P3, uv: (f64, f64), size: f64) -> bool {
        let d = self.o.linear_deflection;
        let at = self.s.point(uv.0, uv.1);
        if norm(sub(at, p)) <= d {
            return true;
        }
        let Some(q) = self.s.foot(p, uv) else {
            return false;
        };
        let foot = self.s.point(q.0, q.1);
        norm(sub(foot, p)) <= d && norm(sub(foot, at)) <= size
    }

    /// Centroids and edge midpoints of inside triangles that are too far from the surface or
    /// that tilt too far from it, each with a vertex of its triangle or edge, where finding the
    /// point in the triangulation starts.
    fn refinement_points(
        &self,
        cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>,
        inside: &[bool],
        info: &[BPoint],
        checked: &mut Checked,
        fresh: Option<&[FixedVertexHandle]>,
    ) -> Vec<((f64, f64), FixedVertexHandle)> {
        let a = self.o.angular_deflection;
        let mut out = Vec::new();
        let cos_a = a.cos();
        checked.normals.resize(info.len(), None);
        let faces: Vec<FixedFaceHandle<InnerTag>> = match fresh {
            None => cdt.inner_faces().map(|f| f.fix()).collect(),
            Some(points) => {
                let mut seen = FxSet::default();
                points
                    .iter()
                    .flat_map(|&v| cdt.vertex(v).out_edges())
                    .filter_map(|e| e.face().as_inner())
                    .map(|f| f.fix())
                    .filter(|f| seen.insert(f.index()))
                    .collect()
            }
        };
        for face in faces.into_iter().map(|f| cdt.face(f)) {
            if !inside[face.fix().index()] {
                continue;
            }
            // A triangle with two copies of a pole is a sliver of the fan: nothing to refine. Two
            // copies of a seam point, though, mean the triangle spans the whole face.
            let hs = face.vertices().map(|h| h.fix().index());
            let v = hs.map(|h| info[h]);
            if (0..3).any(|i| v[i].pole && v[i].r == v[(i + 1) % 3].r) {
                continue;
            }
            let mut key = hs.map(|h| h as u32);
            key.sort_unstable();
            if checked.triangles.contains(&key) {
                continue;
            }
            let uv_c = (
                (v[0].uv.0 + v[1].uv.0 + v[2].uv.0) / 3.0,
                (v[0].uv.1 + v[1].uv.1 + v[2].uv.1) / 3.0,
            );
            let p_c = mul(add(add(v[0].p, v[1].p), v[2].p), 1.0 / 3.0);
            let size = (0..3)
                .map(|i| norm(sub(v[i].p, v[(i + 1) % 3].p)))
                .fold(0.0, f64::max);
            // Corners on the surface (or on its boundary polylines) within the deflection of
            // each other leave the triangle within the deflection too; refining such slivers
            // along the boundary only makes thinner ones.
            // A triangle flat in parameter space (a sliver of points a hair apart) has nothing
            // inside to insert, and points there break the triangulation's search.
            let [p0, p1, p2] = face.vertices().map(|h| h.position());
            let area = ((p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)).abs();
            let longest = [(p0, p1), (p1, p2), (p2, p0)]
                .iter()
                .map(|(x, y)| (x.x - y.x).powi(2) + (x.y - y.y).powi(2))
                .fold(0.0, f64::max);
            if size <= self.o.linear_deflection || area <= FLAT * longest {
                checked.triangles.insert(key);
                continue;
            }
            let near = self.near_surface(p_c, uv_c, size);
            // The surface turns too much between the corners (pole copies excluded: there the
            // normal depends on the side the pole is approached from).
            let normals: Vec<P3> = (0..3)
                .filter(|&i| !v[i].pole)
                .map(|i| {
                    *checked.normals[hs[i]]
                        .get_or_insert_with(|| self.s.normal(v[i].uv.0, v[i].uv.1))
                })
                .collect();
            // On a triangle this small the normals can only turn that far at a crease of the
            // surface (a C0 knot line), which no refinement smooths.
            let turns = size > CREASE_SIZE * self.o.linear_deflection
                && (0..normals.len())
                    .any(|i| (i + 1..normals.len()).any(|j| dot(normals[i], normals[j]) < cos_a));
            if !near || turns {
                out.push((uv_c, face.vertices()[0].fix()));
            } else {
                checked.triangles.insert(key);
            }
            for e in face.adjacent_edges() {
                if cdt.is_constraint_edge(e.fix().as_undirected()) {
                    continue;
                }
                let ends = e.vertices().map(|h| h.fix().index() as u32);
                let key = [ends[0].min(ends[1]), ends[0].max(ends[1])];
                if !checked.edges.insert(key) {
                    continue; // passed before, or already queued this round
                }
                let [x, y] = ends.map(|h| info[h as usize]);
                let uv_m = ((x.uv.0 + y.uv.0) / 2.0, (x.uv.1 + y.uv.1) / 2.0);
                if !self.near_surface(mul(add(x.p, y.p), 0.5), uv_m, norm(sub(x.p, y.p))) {
                    out.push((uv_m, e.from().fix()));
                }
            }
        }
        out
    }
}

/// What the refinement already knows: triangles and edges (by vertex handles) that met the
/// tolerances, and the surface normal at each vertex.
#[derive(Default)]
struct Checked {
    triangles: FxSet<[u32; 3]>,
    edges: FxSet<[u32; 2]>,
    normals: Vec<Option<P3>>,
}

/// Whether `p` lies within `ON_EDGE` of the segment's length of the segment `a`-`b`.
fn hair_off(p: Point2<f64>, a: Point2<f64>, b: Point2<f64>) -> bool {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    let cross = (p.x - a.x) * dy - (p.y - a.y) * dx;
    cross.abs() <= ON_EDGE * len2
}

/// Scaled parameters as a triangulation point. The triangulation rejects magnitudes below
/// 2^-142 other than zero; exporters write such values (1e-45) for zero.
fn cdt_point(uv: (f64, f64), su: f64, sv: f64) -> Point2<f64> {
    let snap = |x: f64| if x.abs() < 1e-30 { 0.0 } else { x };
    Point2::new(snap(uv.0 * su), snap(uv.1 * sv))
}

/// Inner faces inside the boundary: an odd number of constraint edges away from the outside.
fn classify(cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>) -> Vec<bool> {
    // Per face index: unvisited, even (outside) or odd (inside).
    let mut parity: Vec<Option<bool>> = vec![None; cdt.num_all_faces()];
    let mut queue = VecDeque::new();
    for face in cdt.inner_faces() {
        for e in face.adjacent_edges() {
            if e.rev().face().is_outer() && parity[face.fix().index()].is_none() {
                parity[face.fix().index()] = Some(cdt.is_constraint_edge(e.fix().as_undirected()));
                queue.push_back(face.fix());
            }
        }
    }
    while let Some(f) = queue.pop_front() {
        let odd = parity[f.index()] == Some(true);
        for e in cdt.face(f).adjacent_edges() {
            let Some(g) = e.rev().face().as_inner() else {
                continue;
            };
            if parity[g.fix().index()].is_some() {
                continue;
            }
            parity[g.fix().index()] = Some(odd != cdt.is_constraint_edge(e.fix().as_undirected()));
            queue.push_back(g.fix());
        }
    }
    parity.into_iter().map(|p| p == Some(true)).collect()
}

/// Lowest and highest parameters of a loop.
fn bounds(l: &[BPoint]) -> ((f64, f64), (f64, f64)) {
    l.iter().fold(
        ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN)),
        |(lo, hi), b| {
            (
                (lo.0.min(b.uv.0), lo.1.min(b.uv.1)),
                (hi.0.max(b.uv.0), hi.1.max(b.uv.1)),
            )
        },
    )
}

/// Points of a loop tested against another for containment.
const INSIDE_SAMPLES: usize = 8;

/// Whether most of a few points of `l`, moved by `shift`, lie inside `outer` (even-odd rule).
fn mostly_inside(l: &[BPoint], shift: (f64, f64), outer: &[BPoint]) -> bool {
    let step = (l.len() / INSIDE_SAMPLES).max(1);
    let samples: Vec<(f64, f64)> = l
        .iter()
        .step_by(step)
        .map(|b| (b.uv.0 + shift.0, b.uv.1 + shift.1))
        .collect();
    let inside = |(x, y): (f64, f64)| {
        let mut odd = false;
        for i in 0..outer.len() {
            let (a, b) = (outer[i].uv, outer[(i + 1) % outer.len()].uv);
            if (a.1 > y) != (b.1 > y) && x < a.0 + (y - a.1) / (b.1 - a.1) * (b.0 - a.0) {
                odd = !odd;
            }
        }
        odd
    };
    2 * samples.iter().filter(|&&p| inside(p)).count() > samples.len()
}

/// Shoelace area of a loop in parameter space, positive when counter-clockwise.
fn signed_area(l: &[BPoint]) -> f64 {
    let n = l.len();
    (0..n)
        .map(|i| {
            let (a, b) = (l[i].uv, l[(i + 1) % n].uv);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        / 2.0
}
