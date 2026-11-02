use std::f64::consts::{FRAC_PI_2, TAU};

use super::curve::{Curve, angle, dist2, interior_knots};
use super::nurbs::Nurbs2;
use crate::model::{Frame, P3, add, cross, dot, mul, norm, sub, unit};

const NEWTON_STEPS: usize = 60;
/// Grid points the projection starts from when the hint does not lead to the point.
const MAX_STARTS: usize = 4;
/// Parameter step, as a share of the range, for the derivatives of an offset surface's normal.
const NORMAL_STEP: f64 = 1e-6;
/// Curve points sampled for the half-plane of a surface of revolution's curve.
const REVOLUTION_SAMPLES: usize = 16;
const INF: f64 = f64::INFINITY;

/// Face geometry in millimetres. Parametrisations follow ISO 10303-42.
#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    /// `origin + u·x + v·y`
    Plane(Frame),
    /// `origin + r(cos u·x + sin u·y) + v·z`
    Cylinder {
        frame: Frame,
        r: f64,
    },
    /// `origin + (r + v·tan α)(cos u·x + sin u·y) + v·z`
    Cone {
        frame: Frame,
        r: f64,
        half_angle: f64,
    },
    /// `origin + r·cos v(cos u·x + sin u·y) + r·sin v·z`
    Sphere {
        frame: Frame,
        r: f64,
    },
    /// `origin + (R + r·cos v)(cos u·x + sin u·y) + r·sin v·z`
    Torus {
        frame: Frame,
        major: f64,
        minor: f64,
    },
    Nurbs(Nurbs2),
    /// `curve(u) + v·dir`, `dir` of unit length.
    Extrusion {
        curve: Curve,
        dir: P3,
    },
    /// `curve(v)` turned by `u` about the axis (unit `axis_dir`).
    Revolution {
        curve: Curve,
        axis_origin: P3,
        axis_dir: P3,
    },
    /// `basis(u, v) + distance·normal(u, v)`, with the basis surface's parameters.
    Offset {
        basis: Box<Surface>,
        distance: f64,
    },
}

fn ring(f: &Frame, u: f64) -> (P3, P3) {
    let (s, c) = u.sin_cos();
    (
        add(mul(f.x, c), mul(f.y, s)),
        add(mul(f.x, -s), mul(f.y, c)),
    )
}

/// Rotation of `p` by angle `a` about the unit axis `k` through the origin (Rodrigues).
fn rotate(p: P3, k: P3, a: f64) -> P3 {
    let (s, c) = a.sin_cos();
    add(
        add(mul(p, c), mul(cross(k, p), s)),
        mul(k, dot(k, p) * (1.0 - c)),
    )
}

/// Shifts a periodic value to the representative nearest `near`.
fn nearest(x: f64, period: f64, near: f64) -> f64 {
    x + ((near - x) / period).round() * period
}

impl Surface {
    pub fn point(&self, u: f64, v: f64) -> P3 {
        self.derivs(u, v).0
    }

    /// Point and partial derivatives along u and v.
    pub fn derivs(&self, u: f64, v: f64) -> (P3, P3, P3) {
        match self {
            Surface::Plane(f) => (add(f.origin, add(mul(f.x, u), mul(f.y, v))), f.x, f.y),
            Surface::Cylinder { frame: f, r } => {
                let (e, de) = ring(f, u);
                (
                    add(f.origin, add(mul(e, *r), mul(f.z, v))),
                    mul(de, *r),
                    f.z,
                )
            }
            Surface::Cone {
                frame: f,
                r,
                half_angle,
            } => {
                let (e, de) = ring(f, u);
                let t = half_angle.tan();
                let rho = r + v * t;
                (
                    add(f.origin, add(mul(e, rho), mul(f.z, v))),
                    mul(de, rho),
                    add(mul(e, t), f.z),
                )
            }
            Surface::Sphere { frame: f, r } => {
                let (e, de) = ring(f, u);
                let (s, c) = v.sin_cos();
                (
                    add(f.origin, add(mul(e, r * c), mul(f.z, r * s))),
                    mul(de, r * c),
                    add(mul(e, -r * s), mul(f.z, r * c)),
                )
            }
            Surface::Torus {
                frame: f,
                major,
                minor,
            } => {
                let (e, de) = ring(f, u);
                let (s, c) = v.sin_cos();
                let rho = major + minor * c;
                (
                    add(f.origin, add(mul(e, rho), mul(f.z, minor * s))),
                    mul(de, rho),
                    add(mul(e, -minor * s), mul(f.z, minor * c)),
                )
            }
            Surface::Nurbs(n) => n.derivs(u, v),
            Surface::Extrusion { curve, dir } => {
                let (c, dc) = curve.derivs(u);
                (add(c, mul(*dir, v)), dc, *dir)
            }
            Surface::Revolution {
                curve,
                axis_origin,
                axis_dir,
            } => {
                let (c, dc) = curve.derivs(v);
                let q = rotate(sub(c, *axis_origin), *axis_dir, u);
                (
                    add(*axis_origin, q),
                    cross(*axis_dir, q),
                    rotate(dc, *axis_dir, u),
                )
            }
            Surface::Offset { basis, distance } => {
                let (p, du, dv) = basis.derivs(u, v);
                let n = basis.normal(u, v);
                let [(u0, u1), (v0, v1)] = basis.param_bounds();
                let (nu, nv) = (
                    basis.normal_change(u, u0, u1, |x| (x, v)),
                    basis.normal_change(v, v0, v1, |y| (u, y)),
                );
                (
                    add(p, mul(n, *distance)),
                    add(du, mul(nu, *distance)),
                    add(dv, mul(nv, *distance)),
                )
            }
        }
    }

    /// Derivative of the unit normal along one parameter at `x` (range `lo..hi`), by central
    /// differences kept inside the range: the bases give no second derivatives.
    fn normal_change(&self, x: f64, lo: f64, hi: f64, at: impl Fn(f64) -> (f64, f64)) -> P3 {
        let span = if (hi - lo).is_finite() {
            hi - lo
        } else {
            1.0 + x.abs()
        };
        let h = NORMAL_STEP * span;
        let (a, b) = ((x - h).max(lo), (x + h).min(hi));
        if b <= a {
            return [0.0; 3];
        }
        let (pa, pb) = (at(a), at(b));
        mul(
            sub(self.normal(pb.0, pb.1), self.normal(pa.0, pa.1)),
            1.0 / (b - a),
        )
    }

    /// Unit normal, `∂u × ∂v` direction.
    pub fn normal(&self, u: f64, v: f64) -> P3 {
        let analytic = match self {
            Surface::Plane(f) => Some(f.z),
            Surface::Cylinder { frame, .. } => Some(ring(frame, u).0),
            Surface::Sphere { frame, r } => {
                unit(sub(self.point(u, v), frame.origin)).map(|n| mul(n, r.signum()))
            }
            Surface::Offset { basis, .. } => Some(basis.normal(u, v)),
            _ => None,
        };
        if let Some(n) = analytic {
            return n;
        }
        let (_, du, dv) = self.derivs(u, v);
        if let Some(n) = unit(cross(du, dv)) {
            return n;
        }
        // Singular point (cone apex, pole of a B-spline patch): normal of a point next to it.
        let [(u0, u1), (v0, v1)] = self.param_bounds();
        let eps = 1e-7;
        let (u2, v2) = (
            (u + eps).clamp(u0, u1),
            if v + eps <= v1 {
                v + eps
            } else {
                (v - eps).max(v0)
            },
        );
        let (_, du, dv) = self.derivs(u2, v2);
        unit(cross(du, dv)).unwrap_or([0.0, 0.0, 1.0])
    }

    pub fn period(&self) -> (Option<f64>, Option<f64>) {
        match self {
            Surface::Plane(_) => (None, None),
            Surface::Cylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } => {
                (Some(TAU), None)
            }
            Surface::Torus { .. } => (Some(TAU), Some(TAU)),
            Surface::Nurbs(n) => {
                let [(u0, u1), (v0, v1)] = n.range();
                (
                    n.closed_along_u().then_some(u1 - u0),
                    n.closed_along_v().then_some(v1 - v0),
                )
            }
            Surface::Extrusion { curve, .. } => (curve.period(), None),
            Surface::Revolution { curve, .. } => (Some(TAU), curve.period()),
            Surface::Offset { basis, .. } => basis.period(),
        }
    }

    /// Interior knots along u and along v (see [`Curve::breaks`]).
    pub fn breaks(&self) -> [Vec<f64>; 2] {
        match self {
            Surface::Nurbs(n) => {
                let [ru, rv] = n.range();
                [interior_knots(&n.ku, ru), interior_knots(&n.kv, rv)]
            }
            Surface::Extrusion { curve, .. } => [curve.breaks(), Vec::new()],
            Surface::Revolution { curve, .. } => [Vec::new(), curve.breaks()],
            Surface::Offset { basis, .. } => basis.breaks(),
            _ => [Vec::new(), Vec::new()],
        }
    }

    pub fn param_bounds(&self) -> [(f64, f64); 2] {
        match self {
            Surface::Plane(_) => [(-INF, INF), (-INF, INF)],
            Surface::Cylinder { .. } => [(0.0, TAU), (-INF, INF)],
            Surface::Cone { r, half_angle, .. } => {
                let t = half_angle.tan();
                let apex = if t > 0.0 { -r / t } else { -INF };
                [(0.0, TAU), (apex, INF)]
            }
            Surface::Sphere { .. } => [(0.0, TAU), (-FRAC_PI_2, FRAC_PI_2)],
            Surface::Torus { .. } => [(0.0, TAU), (0.0, TAU)],
            Surface::Nurbs(n) => n.range(),
            Surface::Extrusion { curve, .. } => [curve.param_range(), (-INF, INF)],
            Surface::Revolution { curve, .. } => [(0.0, TAU), curve.param_range()],
            Surface::Offset { basis, .. } => basis.param_bounds(),
        }
    }

    /// Parameters of the surface point closest to `p`. Periodic parameters are returned next to
    /// `hint` when one is given. `None` when the iteration does not settle.
    pub fn project(&self, p: P3, hint: Option<(f64, f64)>) -> Option<(f64, f64)> {
        let (u, v) = match self {
            Surface::Plane(f) => {
                let d = sub(p, f.origin);
                (dot(d, f.x), dot(d, f.y))
            }
            Surface::Cylinder { frame: f, .. } => {
                (self.azimuth(f, p, hint), dot(sub(p, f.origin), f.z))
            }
            Surface::Cone {
                frame: f,
                r,
                half_angle,
            } => {
                let d = sub(p, f.origin);
                let h = dot(d, f.z);
                let rho = norm(sub(d, mul(f.z, h)));
                let t = half_angle.tan();
                let v = ((rho - r) * t + h) / (1.0 + t * t);
                (self.azimuth(f, p, hint), v.max(self.param_bounds()[1].0))
            }
            Surface::Sphere { frame: f, .. } => {
                let d = sub(p, f.origin);
                let h = dot(d, f.z);
                let rho = norm(sub(d, mul(f.z, h)));
                (self.azimuth(f, p, hint), h.atan2(rho))
            }
            Surface::Torus {
                frame: f, major, ..
            } => {
                let u = self.azimuth(f, p, hint);
                let (e, _) = ring(f, u);
                let q = sub(sub(p, f.origin), mul(e, *major));
                (u, dot(q, f.z).atan2(dot(q, e)).rem_euclid(TAU))
            }
            Surface::Offset { basis, .. } => {
                // A point of the offset surface is nearest the basis where it was offset from.
                let q = basis.project(p, hint)?;
                self.newton(p, q).unwrap_or(q)
            }
            _ => {
                // From the hint first; from a grid when that fails or stops short of the point.
                let from_hint = hint.and_then(|h| self.newton(p, h));
                let near = |q: (f64, f64)| dist2(self.point(q.0, q.1), p);
                match from_hint {
                    Some(q) if near(q) <= 1e-12 * (1.0 + dot(p, p)) => q,
                    _ => {
                        let from_grid = self
                            .starts(p)
                            .into_iter()
                            .filter_map(|s| self.newton(p, s))
                            .min_by(|a, b| near(*a).total_cmp(&near(*b)));
                        match (from_hint, from_grid) {
                            (Some(a), Some(b)) => {
                                if near(a) <= near(b) {
                                    a
                                } else {
                                    b
                                }
                            }
                            (a, b) => a.or(b)?,
                        }
                    }
                }
            }
        };
        let (pu, pv) = self.period();
        Some(match hint {
            Some((hu, hv)) => (
                pu.map_or(u, |per| nearest(u, per, hu)),
                pv.map_or(v, |per| nearest(v, per, hv)),
            ),
            None => (u, v),
        })
    }

    /// Angle around the frame's axis; on the axis itself the hint's angle (or 0) is used.
    fn azimuth(&self, f: &Frame, p: P3, hint: Option<(f64, f64)>) -> f64 {
        let d = sub(p, f.origin);
        let radial = sub(d, mul(f.z, dot(d, f.z)));
        if norm(radial) < 1e-12 {
            return hint.map_or(0.0, |h| h.0);
        }
        angle(f, p, 1.0, 1.0)
    }

    /// Starting parameters for the iteration: the nearest points of a grid over the bounds, at
    /// most one per position (all the parameters of a pole map to one point, and starting there
    /// leaves no direction to move in).
    fn starts(&self, p: P3) -> Vec<(f64, f64)> {
        if let Surface::Extrusion { curve, dir } = self {
            let u = curve.project(p);
            return vec![(u, dot(sub(p, curve.point(u)), *dir))];
        }
        if let Surface::Revolution {
            curve,
            axis_origin,
            axis_dir,
        } = self
        {
            // Turn `p` into the half-plane of the curve, project it onto the curve there, then
            // measure the angle between the two. The half-plane is that of the curve's point
            // farthest from the axis: the curve may start or end on it (a pole).
            let (t0, t1) = curve.param_range();
            let (a, b) = (t0.max(-1.0), t1.min(1.0).max(t0.max(-1.0)));
            let off_axis = |q: P3| {
                let d = sub(q, *axis_origin);
                norm(sub(d, mul(*axis_dir, dot(d, *axis_dir))))
            };
            let reference = (0..=REVOLUTION_SAMPLES)
                .map(|k| curve.point(a + (b - a) * k as f64 / REVOLUTION_SAMPLES as f64))
                .max_by(|x, y| off_axis(*x).total_cmp(&off_axis(*y)))
                .unwrap_or(*axis_origin);
            let frame = Frame::new(*axis_origin, *axis_dir, sub(reference, *axis_origin));
            let Some(frame) = frame else {
                return vec![(0.0, 0.0)];
            };
            let turn = angle(&frame, p, 1.0, 1.0);
            let v = curve.project(add(
                *axis_origin,
                rotate(sub(p, *axis_origin), *axis_dir, -turn),
            ));
            let u = turn - angle(&frame, curve.point(v), 1.0, 1.0);
            return vec![(u.rem_euclid(TAU), v)];
        }
        let [(u0, u1), (v0, v1)] = self.param_bounds();
        let n = match self {
            Surface::Nurbs(s) => (2 * s.ctrl.len().max(s.ctrl[0].len())).max(16),
            _ => 32,
        };
        let mut grid = Vec::with_capacity((n + 1) * (n + 1));
        for i in 0..=n {
            for j in 0..=n {
                let (u, v) = (
                    u0 + (u1 - u0) * i as f64 / n as f64,
                    v0 + (v1 - v0) * j as f64 / n as f64,
                );
                let q = self.point(u, v);
                grid.push((dist2(q, p), q, (u, v)));
            }
        }
        grid.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut out: Vec<(P3, (f64, f64))> = Vec::new();
        for (_, q, uv) in grid {
            if out.len() == MAX_STARTS {
                break;
            }
            if out
                .iter()
                .all(|(o, _)| dist2(*o, q) > 1e-18 * (1.0 + dot(q, q)))
            {
                out.push((q, uv));
            }
        }
        out.into_iter().map(|(_, uv)| uv).collect()
    }

    /// Gauss–Newton on the squared distance with step halving; non-periodic parameters stay in
    /// bounds, and a parameter held at its bound by the step leaves the other to move alone.
    fn newton(&self, p: P3, (mut u, mut v): (f64, f64)) -> Option<(f64, f64)> {
        let [(u0, u1), (v0, v1)] = self.param_bounds();
        let (pu, pv) = self.period();
        let keep = |x: f64, lo: f64, hi: f64, per: Option<f64>| {
            if per.is_some() { x } else { x.clamp(lo, hi) }
        };
        let pinned = |x: f64, step: f64, lo: f64, hi: f64, per: Option<f64>| {
            per.is_none() && ((x <= lo && step < 0.0) || (x >= hi && step > 0.0))
        };
        let close = 1e-12 * (1.0 + norm(p));
        for _ in 0..NEWTON_STEPS {
            let (s, su, sv) = self.derivs(u, v);
            let r = sub(p, s);
            if norm(r) <= close {
                return Some((u, v));
            }
            let (a, b, c) = (dot(su, su), dot(su, sv), dot(sv, sv));
            let (g1, g2) = (dot(su, r), dot(sv, r));
            let det = a * c - b * b;
            let (mut du, mut dv) = if det > 1e-24 * (a * c).max(1e-300) {
                ((c * g1 - b * g2) / det, (a * g2 - b * g1) / det)
            } else if a >= c && a > 0.0 {
                (g1 / a, 0.0)
            } else if c > 0.0 {
                (0.0, g2 / c)
            } else {
                return None;
            };
            if pinned(u, du, u0, u1, pu) {
                (du, dv) = (0.0, if c > 0.0 { g2 / c } else { 0.0 });
            }
            if pinned(v, dv, v0, v1, pv) {
                (du, dv) = (
                    if a > 0.0 && !pinned(u, g1, u0, u1, pu) {
                        g1 / a
                    } else {
                        0.0
                    },
                    0.0,
                );
            }
            let before = dot(r, r);
            let (mut nu, mut nv) = (u, v);
            for _ in 0..30 {
                nu = keep(u + du, u0, u1, pu);
                nv = keep(v + dv, v0, v1, pv);
                if dist2(self.point(nu, nv), p) <= before {
                    break;
                }
                du *= 0.5;
                dv *= 0.5;
            }
            let moved = (nu - u).abs().max((nv - v).abs());
            (u, v) = (nu, nv);
            if moved <= 1e-13 * (1.0 + u.abs().max(v.abs())) {
                break;
            }
        }
        // Settled: on the surface, or where the free directions are at right angles to the
        // remaining distance (a point off a curved surface or beyond its boundary).
        let (s, su, sv) = self.derivs(u, v);
        let r = sub(p, s);
        let nr = norm(r);
        // Within a micrometre: an edge sample a hair off its face (callers check the distance).
        if nr <= 1e-6 * (1.0 + norm(p)) {
            return Some((u, v));
        }
        // Otherwise stationary when the Gauss–Newton step, with the parameters held at a bound
        // dropped, would move the point by less than the arithmetic can resolve.
        let held =
            |x: f64, lo: f64, hi: f64, per: Option<f64>| per.is_none() && (x <= lo || x >= hi);
        let (hu, hv) = (held(u, u0, u1, pu), held(v, v0, v1, pv));
        let (a, b, c) = (dot(su, su), dot(su, sv), dot(sv, sv));
        let (g1, g2) = (dot(su, r), dot(sv, r));
        let det = a * c - b * b;
        let step = match (hu, hv) {
            (true, true) => [0.0; 3],
            (true, false) if c > 0.0 => mul(sv, g2 / c),
            (false, true) if a > 0.0 => mul(su, g1 / a),
            (false, false) if det > 0.0 => add(
                mul(su, (c * g1 - b * g2) / det),
                mul(sv, (a * g2 - b * g1) / det),
            ),
            _ => return None,
        };
        (norm(step) <= 1e-9 * (1.0 + norm(p))).then_some((u, v))
    }
}
