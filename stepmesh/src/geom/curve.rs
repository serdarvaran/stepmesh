use std::f64::consts::TAU;

use super::nurbs::Nurbs1;
use crate::model::{Frame, P3, add, dot, mul, sub};

const NEWTON_STEPS: usize = 50;

/// Edge geometry in millimetres.
#[derive(Debug, Clone, PartialEq)]
pub enum Curve {
    /// `origin + t·dir`, `dir` of unit length.
    Line {
        origin: P3,
        dir: P3,
    },
    Circle {
        frame: Frame,
        r: f64,
    },
    /// Semi-axis `a` along the frame's x, `b` along its y.
    Ellipse {
        frame: Frame,
        a: f64,
        b: f64,
    },
    Nurbs(Nurbs1),
}

impl Curve {
    pub fn point(&self, t: f64) -> P3 {
        self.derivs(t).0
    }

    pub fn derivative(&self, t: f64) -> P3 {
        self.derivs(t).1
    }

    pub fn derivs(&self, t: f64) -> (P3, P3) {
        match self {
            Curve::Line { origin, dir } => (add(*origin, mul(*dir, t)), *dir),
            Curve::Circle { frame, r } => conic(frame, *r, *r, t),
            Curve::Ellipse { frame, a, b } => conic(frame, *a, *b, t),
            Curve::Nurbs(n) => n.derivs(t),
        }
    }

    pub fn param_range(&self) -> (f64, f64) {
        match self {
            Curve::Line { .. } => (f64::NEG_INFINITY, f64::INFINITY),
            Curve::Circle { .. } | Curve::Ellipse { .. } => (0.0, TAU),
            Curve::Nurbs(n) => n.range(),
        }
    }

    /// Period of a closed curve: conics, and B-splines that end where they start.
    pub fn period(&self) -> Option<f64> {
        match self {
            Curve::Line { .. } => None,
            Curve::Circle { .. } | Curve::Ellipse { .. } => Some(TAU),
            Curve::Nurbs(n) => n.closed().then(|| n.range().1 - n.range().0),
        }
    }

    /// Parameter of the curve point closest to `p` (within the parameter range).
    pub fn project(&self, p: P3) -> f64 {
        match self {
            Curve::Line { origin, dir } => dot(sub(p, *origin), *dir),
            Curve::Circle { frame, .. } => angle(frame, p, 1.0, 1.0),
            Curve::Ellipse { frame, a, b } => self.refine(p, angle(frame, p, *a, *b)),
            Curve::Nurbs(n) => {
                let (t0, t1) = n.range();
                let samples = (4 * n.ctrl.len()).max(16);
                let start = (0..=samples)
                    .map(|i| t0 + (t1 - t0) * i as f64 / samples as f64)
                    .min_by(|&a, &b| dist2(self.point(a), p).total_cmp(&dist2(self.point(b), p)))
                    .unwrap_or(t0);
                self.refine(p, start)
            }
        }
    }

    /// Gauss–Newton on the squared distance, with step halving, kept inside the range
    /// (periodic curves wrap).
    fn refine(&self, p: P3, mut t: f64) -> f64 {
        let (t0, t1) = self.param_range();
        for _ in 0..NEWTON_STEPS {
            let (c, d) = self.derivs(t);
            let dd = dot(d, d);
            if dd == 0.0 {
                break;
            }
            let mut step = dot(sub(p, c), d) / dd;
            let before = dist2(c, p);
            let mut next = t;
            for _ in 0..20 {
                next = t + step;
                next = match self.period() {
                    Some(per) => t0 + (next - t0).rem_euclid(per),
                    None => next.clamp(t0, t1),
                };
                if dist2(self.point(next), p) <= before {
                    break;
                }
                step *= 0.5;
            }
            let moved = (next - t).abs();
            t = next;
            if moved <= 1e-15 * (1.0 + t.abs()) {
                break;
            }
        }
        t
    }
}

fn conic(f: &Frame, a: f64, b: f64, t: f64) -> (P3, P3) {
    let (s, c) = t.sin_cos();
    (
        add(f.origin, add(mul(f.x, a * c), mul(f.y, b * s))),
        add(mul(f.x, -a * s), mul(f.y, b * c)),
    )
}

/// Polar angle of `p` in the frame's xy plane with the axes scaled by `a` and `b`, in [0, 2π).
pub(crate) fn angle(f: &Frame, p: P3, a: f64, b: f64) -> f64 {
    let d = sub(p, f.origin);
    (dot(d, f.y) / b).atan2(dot(d, f.x) / a).rem_euclid(TAU)
}

pub(crate) fn dist2(a: P3, b: P3) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}
