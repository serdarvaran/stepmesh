//! Seams inside a face: where a loop going round a periodic surface is joined to a pole or to
//! the other loop of a band, and where the seams of a whole parameter domain go.

use super::face::BPoint;

/// A loop going round by `du` along u, restarted at point `i`: the points before it follow the
/// end, a period on.
pub(super) fn restart_at(l: &[BPoint], du: f64, i: usize) -> Vec<BPoint> {
    let mut out = l[i..].to_vec();
    out.extend(l[..i].iter().map(|b| BPoint {
        uv: (b.uv.0 + du, b.uv.1),
        ..*b
    }));
    out
}

/// Candidates tried for the start of a seam before falling back to the loop's first point.
const SEAM_CANDIDATES: usize = 256;

/// The point of loop `l` (going round a surface of period `p` along u) where a seam straight
/// along v to the row `v` (the pole, or the other loop of a band) starts: the nearest to that
/// row whose seam crosses neither `l` nor `others`.
pub(super) fn seam_start(l: &[BPoint], v: f64, others: &[Vec<BPoint>], p: f64) -> usize {
    let mut order: Vec<usize> = (0..l.len()).filter(|&i| !l[i].pole).collect();
    order.sort_by(|&a, &b| (l[a].uv.1 - v).abs().total_cmp(&(l[b].uv.1 - v).abs()));
    let mut all = loop_segments(l, p);
    for o in others {
        all.extend(loop_segments(o, p));
    }
    order
        .into_iter()
        .take(SEAM_CANDIDATES)
        .find(|&i| {
            let (u, from) = l[i].uv;
            !all.iter().any(|&(a, b)| crosses_seam(u, from, v, a, b, p))
        })
        .unwrap_or(0)
}

/// Segments of a loop in parameter space; a loop going round (period `p` along u) closes on its
/// first point a period on.
fn loop_segments(c: &[BPoint], p: f64) -> Vec<((f64, f64), (f64, f64))> {
    let n = c.len();
    let turn = ((c[n - 1].uv.0 - c[0].uv.0) / p).round() * p;
    (0..n)
        .map(|i| {
            let next = c[(i + 1) % n].uv;
            let shift = if i + 1 == n { turn } else { 0.0 };
            (c[i].uv, (next.0 + shift, next.1))
        })
        .collect()
}

/// Index of the regular point of `l` nearest `u` along u, modulo the period `p`.
pub(super) fn nearest_along_u(l: &[BPoint], u: f64, p: f64) -> usize {
    let gap = |x: f64| {
        let d = (x - u).rem_euclid(p);
        d.min(p - d)
    };
    (0..l.len())
        .filter(|&i| !l[i].pole)
        .min_by(|&x, &y| gap(l[x].uv.0).total_cmp(&gap(l[y].uv.0)))
        .unwrap_or(0)
}

/// Points of the two loops of a band (`a` index, `b` index) joined by a seam straight along v
/// that crosses no loop or hole. Candidates go along the loop spanning more of u (a loop
/// winding several turns, as round a thread), nearest the other loop first, each with the
/// other loop's point at the same angle.
pub(super) fn seam_pair(
    a: &[BPoint],
    b: &[BPoint],
    holes: &[Vec<BPoint>],
    p: f64,
) -> Option<(usize, usize)> {
    let span = |l: &[BPoint]| {
        let (lo, hi) = l.iter().fold((f64::MAX, f64::MIN), |(x, y), q| {
            (x.min(q.uv.0), y.max(q.uv.0))
        });
        hi - lo
    };
    let swapped = span(a) > span(b);
    let (along, other) = if swapped { (a, b) } else { (b, a) };
    let other_v = other.iter().map(|q| q.uv.1).sum::<f64>() / other.len() as f64;
    let mut order: Vec<usize> = (0..along.len()).filter(|&j| !along[j].pole).collect();
    order.sort_by(|&x, &y| {
        (along[x].uv.1 - other_v)
            .abs()
            .total_cmp(&(along[y].uv.1 - other_v).abs())
    });
    let mut segments = Vec::new();
    for l in [a, b].into_iter().chain(holes.iter().map(Vec::as_slice)) {
        segments.extend(loop_segments(l, p));
    }
    order.into_iter().take(SEAM_CANDIDATES).find_map(|j| {
        let o = nearest_along_u(other, along[j].uv.0, p);
        let (from, to) = (along[j].uv, other[o].uv);
        let clear = !segments
            .iter()
            .any(|&(x, y)| crosses_seam(from.0, from.1, to.1, x, y, p));
        clear.then_some(if swapped { (j, o) } else { (o, j) })
    })
}

/// Whether segment `a`-`b` (or a copy a whole period `p` along u) meets the seam at `u` from
/// `from` to `to` strictly between its ends.
fn crosses_seam(u: f64, from: f64, to: f64, a: (f64, f64), b: (f64, f64), p: f64) -> bool {
    let (lo, hi) = (from.min(to), from.max(to));
    let eps = 1e-9 * (1.0 + hi.abs().max(lo.abs()));
    let inside = |y: f64| y > lo + eps && y < hi - eps;
    let (umin, umax) = (a.0.min(b.0), a.0.max(b.0));
    let k0 = ((u - umax) / p).ceil() as i64;
    let k1 = ((u - umin) / p).floor() as i64;
    (k0..=k1).any(|k| {
        let (a0, b0) = (a.0 + k as f64 * p, b.0 + k as f64 * p);
        if (b0 - a0).abs() <= eps {
            // Along the seam: overlapping it beyond its ends counts.
            (a0 - u).abs() <= eps && a.1.max(b.1) > lo + eps && a.1.min(b.1) < hi - eps
        } else {
            let t = (u - a0) / (b0 - a0);
            (0.0..=1.0).contains(&t) && inside(a.1 + t * (b.1 - a.1))
        }
    })
}

/// A parameter along `axis` (0: u, 1: v), of period `p`, that no loop reaches: the middle of
/// the widest gap between the loops' ranges.
pub(super) fn free_param(loops: &[Vec<BPoint>], axis: usize, p: f64) -> f64 {
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
