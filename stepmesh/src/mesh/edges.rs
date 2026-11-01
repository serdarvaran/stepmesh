//! Edge polylines. Each EDGE_CURVE is sampled once; both faces using it share the samples.

use crate::geom::Curve;
use crate::model::{Edge, P3, dot, norm, sub};
use crate::{Error, Options};

/// Bisection depth limit for B-spline edges (2^24 segments at most).
const MAX_DEPTH: usize = 24;
/// Samples of one conic edge beyond which the input is taken as out of range.
const MAX_EDGE_SAMPLES: usize = 10_000_000;

/// Polyline from `e.start` to `e.end` (both exact) within the chord and angle tolerances.
pub fn discretize_edge(e: &Edge, o: &Options) -> Result<Vec<P3>, Error> {
    let (t_start, t_end) = param_span(e);
    let mut ts = match &e.curve {
        Curve::Line { .. } => vec![t_start, t_end],
        Curve::Circle { r, .. } => uniform(t_start, t_end, *r, o)?,
        Curve::Ellipse { a, b, .. } => uniform(t_start, t_end, a.max(*b), o)?,
        Curve::Nurbs(_) => {
            let (lo, hi) = (t_start.min(t_end), t_start.max(t_end));
            let mut cuts: Vec<f64> = e
                .curve
                .breaks()
                .into_iter()
                .filter(|&k| k > lo && k < hi)
                .collect();
            if t_end < t_start {
                cuts.reverse();
            }
            sample_curve(&|t| e.curve.derivs(t), t_start, t_end, &cuts, o)
        }
    };
    if ts.len() < 2 {
        ts = vec![t_start, t_end];
    }
    // Interior samples that coincide with a neighbour (a knot cut a hair from an end) are
    // dropped: two mesh vertices at one place would leave a gap between the faces.
    let same = |a: P3, b: P3| norm(sub(a, b)) <= 1e-9 * (1.0 + norm(a));
    let mut pts = vec![e.start];
    for &t in &ts[1..ts.len() - 1] {
        let p = e.curve.point(t);
        if !same(p, pts[pts.len() - 1]) && !same(p, e.end) {
            pts.push(p);
        }
    }
    pts.push(e.end);
    Ok(pts)
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
fn uniform(t0: f64, t1: f64, r: f64, o: &Options) -> Result<Vec<f64>, Error> {
    let chord_step = if o.linear_deflection < r {
        2.0 * (1.0 - o.linear_deflection / r).acos()
    } else {
        f64::INFINITY
    };
    let step = chord_step
        .min(o.angular_deflection)
        .min(std::f64::consts::FRAC_PI_2);
    let n = ((t1 - t0).abs() / step).ceil().max(1.0);
    if n > MAX_EDGE_SAMPLES as f64 {
        return Err(Error::TriangleLimit {
            count: n as usize,
            limit: MAX_EDGE_SAMPLES,
        });
    }
    let n = n as usize;
    Ok((0..=n)
        .map(|i| t0 + (t1 - t0) * i as f64 / n as f64)
        .collect())
}

/// Parameters from `a` to `b` (both included) of a curve given by point and derivative.
/// `cuts`, strictly between `a` and `b` in order from `a`, are where the curve may stop being
/// smooth (knots): sampling starts from them, since a curve made of several pieces can return
/// to its chord at every test point of one bisection step.
pub(crate) fn sample_curve(
    c: &dyn Fn(f64) -> (P3, P3),
    a: f64,
    b: f64,
    cuts: &[f64],
    o: &Options,
) -> Vec<f64> {
    let mut out = vec![a];
    let mut from = a;
    for &t in cuts.iter().chain([&b]) {
        adaptive(c, from, t, o, 0, &mut out);
        from = t;
    }
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
