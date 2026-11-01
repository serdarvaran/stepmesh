//! Triangulation of one face in its (u, v) parameter space.

use std::collections::{HashMap, VecDeque};

use spade::handles::FixedVertexHandle;
use spade::{ConstrainedDelaunayTriangulation, Point2, PositionInTriangulation, Triangulation};

use super::edges::sample_curve;
use crate::fxhash::FxSet;
use crate::geom::Surface;
use crate::model::{Face, Loop, P3, add, dot, mul, norm, sub};
use crate::{Error, Options};

/// Refinement passes over the face; each inserts the centroid or edge midpoints of every
/// triangle that is still too far from the surface.
const MAX_REFINE_ROUNDS: usize = 40;
/// Ratio of the parameter derivatives below which a direction counts as collapsed.
const COLLAPSED: f64 = 1e-5;
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
    let mut loops = Vec::new();
    let edge_loops = f
        .loops
        .iter()
        .filter(|(l, _)| matches!(l, Loop::Edges(_)))
        .count();
    for (index, (l, orientation)) in f.loops.iter().enumerate() {
        if let Loop::Edges(list) = l {
            let pts = loop_points(list, *orientation, edges)
                .ok_or_else(|| Ctx::error(f.id, "a loop edge has no polyline"))?;
            let outer = match f.outer {
                Some(o) => Some(o == index),
                None => (edge_loops == 1).then_some(true),
            };
            loops.push(c.map_loop(&pts, outer)?);
        }
    }
    if loops.is_empty() {
        loops = c.full_domain((0.0, 0.0))?;
    } else {
        loops = c.close_wraps(loops)?;
        // A face can lie outside all its loops on a closed surface (a sphere with holes): the
        // whole parameter domain is then its outer boundary.
        let side = if f.same_sense { 1.0 } else { -1.0 };
        let (pu, pv) = c.s.period();
        if let Some(pu) = pu {
            if loops.iter().all(|l| side * signed_area(l) < 0.0) {
                // Its seams go where no loop reaches, so that the loops lie inside it.
                let start = (
                    free_param(&loops, 0, pu),
                    pv.map_or(0.0, |p| free_param(&loops, 1, p)),
                );
                loops.extend(c.full_domain(start)?);
            }
        }
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
    fn error(face: u64, reason: &'static str) -> Error {
        Error::Projection { face, reason }
    }

    fn fail(&self, reason: &'static str) -> Error {
        Ctx::error(self.f.id, reason)
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
        // not close, or closes the wrong way round, try the other choices.
        let Some(period) = self.s.period().0 else {
            return Ok(plain);
        };
        let mut segment = vec![0usize; pts.len()];
        let mut blocks = 0;
        for i in 0..pts.len() {
            if sing[i].0 && (i == 0 || !sing[i - 1].0) {
                blocks += 1;
            }
            segment[i] = blocks;
        }
        if blocks == 0 || blocks > 4 {
            return Ok(plain);
        }
        let want = outer.map(|o| if o == self.f.same_sense { 1.0 } else { -1.0 });
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
            let (a, b) = (neighbour(i, -1), neighbour(i, 1));
            let (from, to) = if sing[i].0 { (a.0, b.0) } else { (a.1, b.1) };
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
        // The seam starts at the first point of the (first) wrapping loop: there no hole may lie.
        let start = |l: Vec<BPoint>, du: f64, plain: &[Vec<BPoint>]| {
            if plain.is_empty() {
                l
            } else {
                restart_near(l, du, free_param(plain, 0, du.abs()))
            }
        };
        match wrapping.len() {
            0 => {}
            1 => {
                let (l, du) = wrapping.pop().ok_or_else(|| self.fail("no loop"))?;
                let l = start(l, du, &plain);
                plain.push(self.close_at_pole(l, du)?);
            }
            2 => {
                let (b, db) = wrapping.pop().ok_or_else(|| self.fail("no loop"))?;
                let (a, da) = wrapping.pop().ok_or_else(|| self.fail("no loop"))?;
                let a = start(a, da, &plain);
                if da.signum() == db.signum() {
                    return Err(self.fail("two loops go round the surface the same way"));
                }
                plain.push(self.join_band(a, da, b, db));
            }
            _ => return Err(self.fail("more than two loops go round the surface")),
        }
        Ok(plain)
    }

    /// The loop, then its first point a period on, then the pole back to the start. The face
    /// lies left of the loop in parameter space (right when it runs against its surface), which
    /// picks the pole.
    fn close_at_pole(&mut self, mut l: Vec<BPoint>, du: f64) -> Result<Vec<BPoint>, Error> {
        let [_, (v0, v1)] = self.s.param_bounds();
        let north = (du > 0.0) == self.f.same_sense;
        let v = if north { v1 } else { v0 };
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
    fn join_band(&mut self, a: Vec<BPoint>, da: f64, mut b: Vec<BPoint>, db: f64) -> Vec<BPoint> {
        if let Some(pv) = self.s.period().1 {
            let mean_v = |l: &[BPoint]| l.iter().map(|p| p.uv.1).sum::<f64>() / l.len() as f64;
            let gap = (mean_v(&a) - mean_v(&b)) / pv;
            let k = if (da > 0.0) == self.f.same_sense {
                gap.floor() + 1.0
            } else {
                gap.ceil() - 1.0
            };
            for p in &mut b {
                p.uv.1 += k * pv;
            }
        }
        let period = da.abs();
        let target = a[0].uv.0 + da;
        let gap = |x: f64| {
            (x - target)
                .rem_euclid(period)
                .min((target - x).rem_euclid(period))
        };
        let k = (0..b.len())
            .min_by(|&i, &j| gap(b[i].uv.0).total_cmp(&gap(b[j].uv.0)))
            .unwrap_or(0);
        let shift = ((target - b[k].uv.0) / period).round() * period;
        let bk = (b[k].uv.0 + shift, b[k].uv.1);
        let seam = self.seam_points(a[0].uv, (bk.0 - da, bk.1));
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
        for l in loops.iter_mut() {
            let (c, _) = centre(l);
            let du = pu.map_or(0.0, |p| ((rc.0 - c.0) / p).round() * p);
            let dv = pv.map_or(0.0, |p| ((rc.1 - c.1) / p).round() * p);
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
        let s = self.s.clone();
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
                    return Err(self.fail("the boundary crosses itself in parameter space"));
                }
                if a != b && !cdt.exists_constraint(a, b) {
                    cdt.add_constraint(a, b);
                }
            }
        }
        let planar = matches!(self.s, Surface::Plane(_));
        let mut checked = Checked::default();
        for _ in 0..if planar { 0 } else { MAX_REFINE_ROUNDS } {
            let inside = classify(&cdt);
            if let Some(limit) = self.o.max_triangles {
                let count = cdt.num_inner_faces();
                if count > limit {
                    return Err(Error::TriangleLimit { count, limit });
                }
            }
            let candidates = self.refinement_points(&cdt, &inside, &info, &mut checked);
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

    /// Centroids and edge midpoints of inside triangles that are too far from the surface or
    /// that tilt too far from the surface.
    fn refinement_points(
        &self,
        cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>,
        inside: &[bool],
        info: &[BPoint],
        checked: &mut Checked,
    ) -> Vec<(f64, f64)> {
        let (d, a) = (self.o.linear_deflection, self.o.angular_deflection);
        let mut out = Vec::new();
        let cos_a = a.cos();
        checked.normals.resize(info.len(), None);
        for face in cdt.inner_faces() {
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
            let off = norm(sub(self.s.point(uv_c.0, uv_c.1), p_c));
            // The surface turns too much between the corners (pole copies excluded: there the
            // normal depends on the side the pole is approached from).
            let normals: Vec<P3> = (0..3)
                .filter(|&i| !v[i].pole)
                .map(|i| {
                    *checked.normals[hs[i]]
                        .get_or_insert_with(|| self.s.normal(v[i].uv.0, v[i].uv.1))
                })
                .collect();
            let turns = (0..normals.len())
                .any(|i| (i + 1..normals.len()).any(|j| dot(normals[i], normals[j]) < cos_a));
            if off > d || turns {
                out.push(uv_c);
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
                if norm(sub(self.s.point(uv_m.0, uv_m.1), mul(add(x.p, y.p), 0.5))) > d {
                    out.push(uv_m);
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

/// A loop going round by `du` along u, restarted at its regular point nearest `u` (modulo the
/// period): the points before it follow the end, a period on.
fn restart_near(l: Vec<BPoint>, du: f64, u: f64) -> Vec<BPoint> {
    let p = du.abs();
    let dist = |x: f64| {
        let d = (x - u).rem_euclid(p);
        d.min(p - d)
    };
    let Some(i) = (0..l.len())
        .filter(|&i| !l[i].pole)
        .min_by(|&a, &b| dist(l[a].uv.0).total_cmp(&dist(l[b].uv.0)))
    else {
        return l;
    };
    let mut out = l[i..].to_vec();
    out.extend(l[..i].iter().map(|b| BPoint {
        uv: (b.uv.0 + du, b.uv.1),
        ..*b
    }));
    out
}

/// A parameter along `axis` (0: u, 1: v), of period `p`, that no loop reaches: the middle of
/// the widest gap between the loops' ranges.
fn free_param(loops: &[Vec<BPoint>], axis: usize, p: f64) -> f64 {
    let at = |b: &BPoint| if axis == 0 { b.uv.0 } else { b.uv.1 };
    let mut spans: Vec<(f64, f64)> = loops
        .iter()
        .map(|l| {
            let (lo, hi) = l
                .iter()
                .map(at)
                .fold((f64::MAX, f64::MIN), |(a, b), x| (a.min(x), b.max(x)));
            let s = lo.rem_euclid(p);
            (s, s + hi - lo)
        })
        .collect();
    spans.sort_by(|a, b| a.0.total_cmp(&b.0));
    let Some(&(first, mut reach)) = spans.first() else {
        return 0.0;
    };
    let (mut widest, mut middle) = (f64::MIN, first);
    for &(s, e) in &spans[1..] {
        if s - reach > widest {
            (widest, middle) = (s - reach, 0.5 * (s + reach));
        }
        reach = reach.max(e);
    }
    // The gap from the last loop round to the first.
    if first + p - reach > widest {
        middle = 0.5 * (reach + first + p);
    }
    middle
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
