//! Edge polylines. Each EDGE_CURVE is sampled once; both faces using it share the samples.

use crate::geom::Curve;
use crate::model::{Edge, P3, dot, norm, sub};
use crate::{Error, Options};

/// Bisection depth limit for B-spline edges (2^24 segments at most).
const MAX_DEPTH: usize = 24;

/// Polyline from `e.start` to `e.end` (both exact) within the chord and angle tolerances.
pub fn discretize_edge(e: &Edge, o: &Options) -> Result<Vec<P3>, Error> {
    let (t_start, t_end) = param_span(e);
    let mut ts = match &e.curve {
        Curve::Line { .. } => vec![t_start, t_end],
        Curve::Circle { r, .. } => uniform(t_start, t_end, *r, o),
        Curve::Ellipse { a, b, .. } => uniform(t_start, t_end, a.max(*b), o),
        Curve::Nurbs(n) => {
            let mut cuts = vec![t_start];
            let (lo, hi) = (t_start.min(t_end), t_start.max(t_end));
            let mut knots: Vec<f64> = n
                .knots
                .iter()
                .copied()
                .filter(|&k| k > lo && k < hi)
                .collect();
            knots.dedup();
            if t_end < t_start {
                knots.reverse();
            }
            cuts.extend(knots);
            cuts.push(t_end);
            let mut ts = vec![t_start];
            for w in cuts.windows(2) {
                adaptive(&|t| e.curve.derivs(t), w[0], w[1], o, 0, &mut ts);
            }
            ts
        }
    };
    if ts.len() < 2 {
        ts = vec![t_start, t_end];
    }
    let last = ts.len() - 1;
    Ok(ts
        .iter()
        .enumerate()
        .map(|(i, &t)| match i {
            0 => e.start,
            i if i == last => e.end,
            _ => e.curve.point(t),
        })
        .collect())
}

/// Curve parameters of the edge's start and end, in the direction the edge runs.
fn param_span(e: &Edge) -> (f64, f64) {
    let ts = e.curve.project(e.start);
    let closed = e.start_vertex == e.end_vertex;
    match e.curve.period() {
        Some(per) => {
            let te = if closed { ts } else { e.curve.project(e.end) };
            if e.same_sense {
                let te = if te <= ts { te + per } else { te };
                (ts, if closed { ts + per } else { te })
            } else {
                let te = if te >= ts { te - per } else { te };
                (ts, if closed { ts - per } else { te })
            }
        }
        None if closed => {
            let (t0, t1) = e.curve.param_range();
            if e.same_sense { (t0, t1) } else { (t1, t0) }
        }
        None => (ts, e.curve.project(e.end)),
    }
}

/// Equal parameter steps for a conic whose largest radius is `r`.
fn uniform(t0: f64, t1: f64, r: f64, o: &Options) -> Vec<f64> {
    let chord_step = if o.linear_deflection < r {
        2.0 * (1.0 - o.linear_deflection / r).acos()
    } else {
        f64::INFINITY
    };
    let step = chord_step
        .min(o.angular_deflection)
        .min(std::f64::consts::FRAC_PI_2);
    let n = ((t1 - t0).abs() / step).ceil().max(1.0) as usize;
    (0..=n)
        .map(|i| t0 + (t1 - t0) * i as f64 / n as f64)
        .collect()
}

/// Parameters from `a` to `b` (both included) of a curve given by point and derivative.
pub(crate) fn sample_curve(c: &dyn Fn(f64) -> (P3, P3), a: f64, b: f64, o: &Options) -> Vec<f64> {
    let mut out = vec![a];
    adaptive(c, a, b, o, 0, &mut out);
    out
}

/// Appends the parameters after `a` up to and including `b`, bisecting where the curve leaves
/// the chord by more than the deflection or turns by more than the angle.
fn adaptive(
    c: &dyn Fn(f64) -> (P3, P3),
    a: f64,
    b: f64,
    o: &Options,
    depth: usize,
    out: &mut Vec<f64>,
) {
    let (pa, da) = c(a);
    let (pb, db) = c(b);
    let chord = sub(pb, pa);
    let len2 = dot(chord, chord);
    let off_chord = |t: f64| {
        let q = sub(c(t).0, pa);
        let along = if len2 > 0.0 {
            dot(q, chord) / len2
        } else {
            0.0
        };
        norm(sub(
            q,
            [chord[0] * along, chord[1] * along, chord[2] * along],
        ))
    };
    let deviation = [0.25, 0.5, 0.75]
        .iter()
        .map(|f| off_chord(a + (b - a) * f))
        .fold(0.0, f64::max);
    let (na, nb) = (norm(da), norm(db));
    let turn = if na > 0.0 && nb > 0.0 {
        (dot(da, db) / (na * nb)).clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    if depth < MAX_DEPTH && (deviation > o.linear_deflection || turn > o.angular_deflection) {
        let m = 0.5 * (a + b);
        adaptive(c, a, m, o, depth + 1, out);
        adaptive(c, m, b, o, depth + 1, out);
    } else {
        out.push(b);
    }
}
