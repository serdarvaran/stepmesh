//! Certified bounds on the distance of edge chords and mesh triangles from B-spline, swept and
//! offset geometry.

use std::cell::RefCell;
use std::collections::HashMap;

use super::nurbs::find_span;
use super::{Curve, Nurbs1, Surface};
use crate::model::{P3, add, dot, mul, norm, sub, unit};
use crate::{Error, Options};

type H = [f64; 4];
type Uv = (f64, f64);
/// Control points of an error numerator and its weights.
type ErrorNet = (Vec<Vec<P3>>, Vec<Vec<f64>>);

#[derive(Clone, Copy)]
struct Section {
    span: usize,
    lo: f64,
    hi: f64,
    a: f64,
    b: f64,
}

fn sections(k: &[f64], p: usize, n: usize, closed: bool, lo: f64, hi: f64) -> Option<Vec<Section>> {
    if !lo.is_finite() || !hi.is_finite() || hi < lo {
        return None;
    }
    let (start, end) = (k[p], k[n]);
    let mut out = Vec::new();
    let mut x = lo;
    loop {
        if out.len() >= 10_000 {
            return None;
        }
        let mut offset = if closed {
            ((x - start) / (end - start)).floor() * (end - start)
        } else {
            0.
        };
        let mut span = find_span(n - 1, p, x - offset, k);
        // Outside the first period `x - offset` can round below the knot that ended the last
        // section, back into its span.
        while closed && k[span + 1] + offset <= x {
            if span + 1 < n {
                span += 1;
            } else {
                offset += end - start;
                span = find_span(n - 1, p, start, k);
            }
        }
        let a = if closed {
            (x - offset).max(k[span])
        } else {
            x - offset
        };
        let stop = if closed || span < n - 1 {
            k[span + 1] + offset
        } else {
            hi
        };
        let y = if stop > x { hi.min(stop) } else { hi };
        out.push(Section {
            span,
            lo: x,
            hi: y,
            a,
            b: y - offset,
        });
        if y >= hi {
            break;
        }
        if y <= x {
            return None;
        }
        x = y;
    }
    Some(out)
}

// Evaluating the blossom at a^(p-j), b^j extracts a span without knot insertion.
fn extract(ctrl: &[H], p: usize, k: &[f64], s: Section) -> Vec<H> {
    extract_local(&ctrl[s.span - p..=s.span], p, k, s)
}

fn extract_local(ctrl: &[H], p: usize, k: &[f64], s: Section) -> Vec<H> {
    let mut d = ctrl.to_vec();
    (0..=p)
        .map(|j| {
            d.copy_from_slice(ctrl);
            for r in 1..=p {
                let t = if r <= p - j { s.a } else { s.b };
                for q in (r..=p).rev() {
                    let i = s.span - p + q;
                    let alpha = (t - k[i]) / (k[i + p - r + 1] - k[i]);
                    d[q] = std::array::from_fn(|c| (1. - alpha) * d[q - 1][c] + alpha * d[q][c]);
                }
            }
            d[p]
        })
        .collect()
}

fn cart(h: H) -> P3 {
    [h[0] / h[3], h[1] / h[3], h[2] / h[3]]
}
fn valid(h: H) -> bool {
    h[3] > 0. && h.iter().all(|v| v.is_finite())
}

fn segment_distance(p: P3, a: P3, b: P3) -> f64 {
    let d = sub(b, a);
    let l = dot(d, d);
    let t = if l > 0. {
        (dot(sub(p, a), d) / l).clamp(0., 1.)
    } else {
        0.
    };
    norm(sub(sub(p, a), mul(d, t)))
}

fn curve_ok(c: &[H], o: &Options) -> bool {
    if !c.iter().copied().all(valid) {
        return false;
    }
    let pts: Vec<_> = c.iter().copied().map(cart).collect();
    let (a, b) = (pts[0], pts[pts.len() - 1]);
    // A piece a hair long (a knot an ulp from the end) is one point: the directions between its
    // control points are rounding noise. Its hull and chord then lie within twice the hair of
    // each other, so the hair stays below half the deflection.
    let hair = (1e-9 * (1. + norm(a))).min(0.25 * o.linear_deflection);
    if pts.iter().all(|&p| norm(sub(p, a)) <= hair) {
        return true;
    }
    if pts.iter().any(|&p| {
        let d = segment_distance(p, a, b);
        !d.is_finite() || d > o.linear_deflection
    }) {
        return false;
    }
    let dirs: Vec<_> = (0..pts.len())
        .flat_map(|i| ((i + 1)..pts.len()).map(move |j| (i, j)))
        .filter_map(|(i, j)| unit(sub(pts[j], pts[i])))
        .collect();
    // A positive rational Bézier tangent is a positive sum of forward control differences.
    let Some(axis) = unit(dirs.iter().copied().fold([0.; 3], add)) else {
        return dirs.is_empty();
    };
    let cos = (o.angular_deflection.min(std::f64::consts::PI) / 2.).cos();
    dirs.iter().all(|&d| dot(axis, d) >= cos)
}

fn split(c: &[H]) -> (Vec<H>, Vec<H>) {
    let mut d = c.to_vec();
    let (mut left, mut right) = (vec![d[0]], vec![d[d.len() - 1]]);
    for n in (1..d.len()).rev() {
        for j in 0..n {
            d[j] = std::array::from_fn(|k| 0.5 * d[j][k] + 0.5 * d[j + 1][k]);
        }
        left.push(d[0]);
        right.push(d[n - 1]);
    }
    right.reverse();
    (left, right)
}

fn sample(
    c: &[H],
    a: f64,
    b: f64,
    depth: usize,
    o: &Options,
    out: &mut Vec<f64>,
) -> Result<(), Error> {
    if out.len() >= 10_000_000 {
        return Err(Error::TriangleLimit {
            count: out.len() + 1,
            limit: 10_000_000,
        });
    }
    if curve_ok(c, o) {
        out.push(b);
        return Ok(());
    }
    let m = 0.5 * a + 0.5 * b;
    if depth == 24 || m == a || m == b {
        return Err(Error::Invalid {
            id: 0,
            msg: "spline edge does not reach the tolerances".into(),
        });
    }
    let (left, right) = split(c);
    sample(&left, a, m, depth + 1, o, out)?;
    sample(&right, m, b, depth + 1, o, out)
}

pub(crate) fn sample_nurbs(n: &Nurbs1, a: f64, b: f64, o: &Options) -> Result<Vec<f64>, Error> {
    let spans = sections(
        &n.knots,
        n.degree,
        n.ctrl.len(),
        n.closed(),
        a.min(b),
        a.max(b),
    )
    .ok_or_else(|| Error::Invalid {
        id: 0,
        msg: "invalid spline edge parameter span".into(),
    })?;
    let mut out = vec![a.min(b)];
    for s in spans {
        sample(
            &extract(&n.ctrl, n.degree, &n.knots, s),
            s.lo,
            s.hi,
            0,
            o,
            &mut out,
        )?;
    }
    if b < a {
        out.reverse();
    }
    Ok(out)
}

/// Derivative of a Bernstein net on the unit square, along its first index or its second.
fn diff<const K: usize>(c: &[Vec<[f64; K]>], along_u: bool) -> Vec<Vec<[f64; K]>> {
    let (p, q) = (c.len() - 1, c[0].len() - 1);
    let (di, dj, degree) = if along_u { (1, 0, p) } else { (0, 1, q) };
    if degree == 0 {
        // Constant along this direction: the zero polynomial, of degree 0 there.
        return vec![
            vec![[0.; K]; if along_u { q + 1 } else { 1 }];
            if along_u { 1 } else { p + 1 }
        ];
    }
    (0..=p - di)
        .map(|i| {
            (0..=q - dj)
                .map(|j| {
                    std::array::from_fn(|k| degree as f64 * (c[i + di][j + dj][k] - c[i][j][k]))
                })
                .collect()
        })
        .collect()
}

fn binomial(n: usize, k: usize) -> f64 {
    (0..k).fold(1., |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

/// Product of two Bernstein nets on the unit square, their values combined by the bilinear `f`.
fn product<T: Copy, U: Copy>(x: &[Vec<T>], y: &[Vec<U>], f: impl Fn(T, U) -> P3) -> Vec<Vec<P3>> {
    let (m1, n1, m2, n2) = (x.len() - 1, x[0].len() - 1, y.len() - 1, y[0].len() - 1);
    // B_i^m B_k^n = C(m,i) C(n,k) / C(m+n,i+k) B_(i+k)^(m+n), along each direction.
    let scales = |a: usize, b: usize| -> Vec<Vec<f64>> {
        (0..=a)
            .map(|i| {
                (0..=b)
                    .map(|k| binomial(a, i) * binomial(b, k) / binomial(a + b, i + k))
                    .collect()
            })
            .collect()
    };
    let (su, sv) = (scales(m1, m2), scales(n1, n2));
    let mut out = vec![vec![[0.; 3]; n1 + n2 + 1]; m1 + m2 + 1];
    for (i, k) in (0..=m1).flat_map(|i| (0..=m2).map(move |k| (i, k))) {
        for (j, l) in (0..=n1).flat_map(|j| (0..=n2).map(move |l| (j, l))) {
            let v = f(x[i][j], y[k][l]);
            out[i + k][j + l] = add(out[i + k][j + l], mul(v, su[i][k] * sv[j][l]));
        }
    }
    out
}

fn xyz(h: H) -> P3 {
    [h[0], h[1], h[2]]
}

/// Lower bound on |P| for a Bernstein net: its weights are positive and sum to one, so P·axis
/// is at least the smallest control's.
fn least_norm(p: &[Vec<P3>]) -> Option<f64> {
    let axis = unit(p.iter().flatten().copied().fold([0.; 3], add))?;
    let r = p
        .iter()
        .flatten()
        .map(|&c| dot(c, axis))
        .fold(f64::INFINITY, f64::min);
    (r > 0.).then_some(r)
}

fn largest(n: &[Vec<P3>]) -> f64 {
    n.iter().flatten().map(|&c| norm(c)).fold(0., f64::max)
}

struct Patch {
    c: Vec<Vec<H>>,
    u: (f64, f64),
    v: (f64, f64),
}

fn affine(uv: [Uv; 3], pts: [P3; 3], at: Uv) -> Option<P3> {
    let (a, b) = (
        (uv[1].0 - uv[0].0, uv[1].1 - uv[0].1),
        (uv[2].0 - uv[0].0, uv[2].1 - uv[0].1),
    );
    let det = a.0 * b.1 - a.1 * b.0;
    if det == 0. {
        return None;
    }
    let d = (at.0 - uv[0].0, at.1 - uv[0].1);
    let (s, t) = ((d.0 * b.1 - d.1 * b.0) / det, (a.0 * d.1 - a.1 * d.0) / det);
    Some(add(
        pts[0],
        add(mul(sub(pts[1], pts[0]), s), mul(sub(pts[2], pts[0]), t)),
    ))
}

/// The four halves-of-sides triangles; the affine map restricted to each.
fn quarters(uv: [Uv; 3], pts: [P3; 3]) -> [([Uv; 3], [P3; 3]); 4] {
    let mid_uv = |i: usize, j: usize| ((uv[i].0 + uv[j].0) / 2., (uv[i].1 + uv[j].1) / 2.);
    let mid_p = |i: usize, j: usize| mul(add(pts[i], pts[j]), 0.5);
    let us = [
        uv[0],
        uv[1],
        uv[2],
        mid_uv(0, 1),
        mid_uv(1, 2),
        mid_uv(2, 0),
    ];
    let ps = [
        pts[0],
        pts[1],
        pts[2],
        mid_p(0, 1),
        mid_p(1, 2),
        mid_p(2, 0),
    ];
    [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]].map(|ids| (ids.map(|i| us[i]), ids.map(|i| ps[i])))
}

impl Patch {
    fn position_bound(&self, uv: [Uv; 3], pts: [P3; 3]) -> Option<f64> {
        let (num, w) = self.error_net(uv, pts)?;
        let mut bound: f64 = 0.;
        for (n, w) in num.iter().flatten().zip(w.iter().flatten()) {
            let d = norm(*n) / w;
            if !d.is_finite() {
                return None;
            }
            bound = bound.max(d);
        }
        Some(bound)
    }

    /// Second derivative bounds of A - W·L (the patch's numerator less its weight times the
    /// affine map through `uv`, `pts`) along u, across u and v, and along v, with the
    /// smallest weight; `None` for a patch of zero width.
    fn second_derivatives(&self, uv: [Uv; 3], pts: [P3; 3]) -> Option<([f64; 3], f64)> {
        let (du, dv) = (self.u.1 - self.u.0, self.v.1 - self.v.0);
        if du <= 0. || dv <= 0. {
            return None;
        }
        let (num, w) = self.error_net(uv, pts)?;
        let (n, m) = (num.len() - 1, num[0].len() - 1);
        let mut d = [0f64; 3];
        for i in 0..=n {
            for j in 0..=m {
                let e = |a: usize, b: usize| num[i + a][j + b];
                if i + 2 <= n {
                    d[0] = d[0].max(norm(add(sub(e(2, 0), mul(e(1, 0), 2.)), e(0, 0))));
                }
                if i < n && j < m {
                    d[1] = d[1].max(norm(add(sub(e(1, 1), e(1, 0)), sub(e(0, 0), e(0, 1)))));
                }
                if j + 2 <= m {
                    d[2] = d[2].max(norm(add(sub(e(0, 2), mul(e(0, 1), 2.)), e(0, 0))));
                }
            }
        }
        let (n, m) = (n as f64, m as f64);
        let bounds = [
            n * (n - 1.) * d[0] / (du * du),
            n * m * d[1] / (du * dv),
            m * (m - 1.) * d[2] / (dv * dv),
        ];
        let w_min = w.iter().flatten().copied().fold(f64::INFINITY, f64::min);
        (bounds.iter().all(|b| b.is_finite()) && w_min > 0.).then_some((bounds, w_min))
    }

    /// Second derivative bounds of the unit normal N along u, across u and v, and along v. With
    /// the patch A/W, N = P/r for the polynomial P = (A_u W - A W_u) × (A_v W - A W_v), r = |P|.
    /// N_uu = (P_uu - r_uu N - 2 r_u N_u)/r with |N_u| = |ΠP_u|/r (Π across N): its part across
    /// N is at most |P_uu|/r + |P_u|²/r² (2|r_u||ΠP_u| ≤ |P_u|²), its part along N is |N_u|².
    /// So |N_uu| ≤ |P_uu|/r + 2|P_u|²/r², and likewise |N_uv| ≤ |P_uv|/r + 2|P_u||P_v|/r², for
    /// any lower bound r on |P|.
    fn normal_second_derivatives(&self) -> Option<[f64; 3]> {
        let (du, dv) = (self.u.1 - self.u.0, self.v.1 - self.v.0);
        if du <= 0. || dv <= 0. || self.c.len() < 2 || self.c[0].len() < 2 {
            return None;
        }
        if !self.c.iter().flatten().copied().all(valid) {
            return None;
        }
        let tangent = |d: H, h: H| sub(mul(xyz(d), h[3]), mul(xyz(h), d[3]));
        // With equal weights A = W·S, so A_u × A_v has the normal's direction at lower degree.
        let p = if self.c.iter().flatten().all(|h| h[3] == self.c[0][0][3]) {
            product(&diff(&self.c, true), &diff(&self.c, false), |a: H, b: H| {
                crate::model::cross(xyz(a), xyz(b))
            })
        } else {
            product(
                &product(&diff(&self.c, true), &self.c, tangent),
                &product(&diff(&self.c, false), &self.c, tangent),
                crate::model::cross,
            )
        };
        let r = least_norm(&p)?;
        let (pu, pv) = (diff(&p, true), diff(&p, false));
        let (a, b) = (largest(&pu), largest(&pv));
        let bounds = [
            (largest(&diff(&pu, true)) / r + 2. * a * a / (r * r)) / (du * du),
            (largest(&diff(&pu, false)) / r + 2. * a * b / (r * r)) / (du * dv),
            (largest(&diff(&pv, false)) / r + 2. * b * b / (r * r)) / (dv * dv),
        ];
        bounds.iter().all(|b| b.is_finite()).then_some(bounds)
    }

    /// Control net of A - W·L and of W, both raised by one degree in u and v (L is bilinear
    /// on the patch's rectangle).
    fn error_net(&self, uv: [Uv; 3], pts: [P3; 3]) -> Option<ErrorNet> {
        let (p, q) = (self.c.len() - 1, self.c[0].len() - 1);
        if !self.c.iter().flatten().copied().all(valid) {
            return None;
        }
        let mut l = [[[0.; 3]; 2]; 2];
        for (a, row) in l.iter_mut().enumerate() {
            for (b, v) in row.iter_mut().enumerate() {
                *v = affine(
                    uv,
                    pts,
                    (
                        if a == 0 { self.u.0 } else { self.u.1 },
                        if b == 0 { self.v.0 } else { self.v.1 },
                    ),
                )?;
            }
        }
        // S-L = (A-WL)/W. Degree elevation gives a rational control net for this error.
        let mut nums = vec![vec![[0.; 3]; q + 2]; p + 2];
        let mut ws = vec![vec![0.; q + 2]; p + 2];
        for i in 0..=p + 1 {
            for j in 0..=q + 1 {
                let (mut num, mut w) = ([0.; 3], 0.);
                for (a, row) in l.iter().enumerate() {
                    for (b, &value) in row.iter().enumerate() {
                        if i < a || i - a > p || j < b || j - b > q {
                            continue;
                        }
                        let f = (if a == 0 { p + 1 - i } else { i }) as f64 / (p + 1) as f64
                            * (if b == 0 { q + 1 - j } else { j }) as f64
                            / (q + 1) as f64;
                        let h = self.c[i - a][j - b];
                        num = add(num, mul(sub(cart(h), value), f * h[3]));
                        w += f * h[3];
                    }
                }
                if !num.iter().all(|x| x.is_finite()) || !w.is_finite() {
                    return None;
                }
                nums[i][j] = num;
                ws[i][j] = w;
            }
        }
        Some((nums, ws))
    }
}

fn rectangle(uv: [Uv; 3]) -> [(f64, f64); 2] {
    [uv.map(|p| p.0), uv.map(|p| p.1)].map(|a| {
        (
            a.into_iter().fold(f64::INFINITY, f64::min),
            a.into_iter().fold(f64::NEG_INFINITY, f64::max),
        )
    })
}

fn patches(s: &Surface, rect: [(f64, f64); 2]) -> Option<Vec<Patch>> {
    patches_over(s, rect, None)
}

/// With a margin, each section grown by it on both sides (and by at least a sliver of its span)
/// within its knot span: a superset of the clipped one in which no piece is a hair wide, whose
/// rounding would swamp its derivatives.
fn widened(sections: Vec<Section>, k: &[f64], margin: Option<f64>) -> Vec<Section> {
    let Some(margin) = margin else {
        return sections;
    };
    sections
        .into_iter()
        .map(|s| {
            let (lo, hi) = (k[s.span], k[s.span + 1]);
            let m = margin.max((hi - lo) / 1024.);
            let (a, b) = (s.a.min((s.a - m).max(lo)), s.b.max((s.b + m).min(hi)));
            Section {
                span: s.span,
                lo: s.lo + (a - s.a),
                hi: s.hi + (b - s.b),
                a,
                b,
            }
        })
        .collect()
}

fn widths([(u0, u1), (v0, v1)]: [(f64, f64); 2]) -> [f64; 2] {
    [u1 - u0, v1 - v0]
}

fn patches_over(
    s: &Surface,
    [(u0, u1), (v0, v1)]: [(f64, f64); 2],
    margin: Option<[f64; 2]>,
) -> Option<Vec<Patch>> {
    let (mu, mv) = (margin.map(|m| m[0]), margin.map(|m| m[1]));
    let mut out = Vec::new();
    match s {
        Surface::Nurbs(n) => {
            let us = sections(&n.ku, n.du, n.ctrl.len(), n.closed_along_u(), u0, u1)?;
            let vs = sections(&n.kv, n.dv, n.ctrl[0].len(), n.closed_along_v(), v0, v1)?;
            let (us, vs) = (widened(us, &n.ku, mu), widened(vs, &n.kv, mv));
            for u in us {
                for &v in &vs {
                    let cols: Vec<_> = (v.span - n.dv..=v.span)
                        .map(|j| {
                            extract_local(
                                &n.ctrl[u.span - n.du..=u.span]
                                    .iter()
                                    .map(|r| r[j])
                                    .collect::<Vec<_>>(),
                                n.du,
                                &n.ku,
                                u,
                            )
                        })
                        .collect();
                    let c = (0..=n.du)
                        .map(|i| {
                            extract_local(
                                &cols.iter().map(|r| r[i]).collect::<Vec<_>>(),
                                n.dv,
                                &n.kv,
                                v,
                            )
                        })
                        .collect();
                    out.push(Patch {
                        c,
                        u: (u.lo, u.hi),
                        v: (v.lo, v.hi),
                    });
                }
            }
        }
        Surface::Extrusion {
            curve: Curve::Nurbs(n),
            dir,
        } => {
            let us = sections(&n.knots, n.degree, n.ctrl.len(), n.closed(), u0, u1)?;
            for u in widened(us, &n.knots, mu) {
                let c = extract(&n.ctrl, n.degree, &n.knots, u)
                    .into_iter()
                    .map(|h| {
                        [v0, v1]
                            .map(|v| {
                                let p = add(cart(h), mul(*dir, v));
                                [p[0] * h[3], p[1] * h[3], p[2] * h[3], h[3]]
                            })
                            .to_vec()
                    })
                    .collect();
                out.push(Patch {
                    c,
                    u: (u.lo, u.hi),
                    v: (v0, v1),
                });
            }
        }
        _ => return None,
    }
    Some(out)
}

pub(crate) fn has_spline(s: &Surface) -> bool {
    match s {
        Surface::Nurbs(_)
        | Surface::Extrusion {
            curve: Curve::Nurbs(_),
            ..
        }
        | Surface::Revolution {
            curve: Curve::Nurbs(_),
            ..
        } => true,
        Surface::Offset { basis, .. } => has_spline(basis),
        _ => false,
    }
}

pub(crate) fn triangle_within(
    s: &Surface,
    uv: [Uv; 3],
    pts: [P3; 3],
    o: &Options,
    cache: &NormalBounds,
) -> Option<bool> {
    if !has_spline(s) {
        return None;
    }
    if visibly_off(s, uv, pts, o.linear_deflection) {
        return Some(false);
    }
    Some(match s {
        Surface::Revolution {
            curve: Curve::Nurbs(n),
            axis_origin,
            axis_dir,
        } => {
            revolution_hull_within(n, *axis_origin, *axis_dir, uv, pts, o.linear_deflection)
                || revolution_bound(n, *axis_origin, *axis_dir, uv, pts)
                    .is_some_and(|d| d <= o.linear_deflection)
        }
        Surface::Offset { basis, distance } => offset_within(basis, *distance, uv, pts, o, cache),
        _ => spline_within(s, uv, pts, o.linear_deflection, SPLIT_DEPTH),
    })
}

/// Times a triangle is split in four to certify it.
const SPLIT_DEPTH: usize = 3;

/// The triangle bound on the triangle or, failing that, on its quarters.
fn spline_within(s: &Surface, uv: [Uv; 3], pts: [P3; 3], tolerance: f64, depth: usize) -> bool {
    let rect = rectangle(uv);
    // The triangle bound is tight where the surface is smooth; across a crease, or where the
    // control net hugs the surface closer than its derivatives suggest, the rectangle's own.
    let certified = triangle_bound(s, rect, uv, pts).is_some_and(|d| d <= tolerance)
        || patches(s, rect).is_some_and(|ps| {
            ps.iter()
                .all(|p| p.position_bound(uv, pts).is_some_and(|d| d <= tolerance))
        });
    if certified || depth == 0 {
        return certified;
    }
    let quarters = quarters(uv, pts);
    !quarters
        .iter()
        .any(|&(uv, pts)| visibly_off(s, uv, pts, tolerance))
        && quarters
            .into_iter()
            .all(|(uv, pts)| spline_within(s, uv, pts, tolerance, depth - 1))
}

/// Whether the triangle's centroid or a side's midpoint already lies farther than `tolerance`
/// from the surface point at its parameters: no bound can then pass it.
fn visibly_off(s: &Surface, uv: [Uv; 3], pts: [P3; 3], tolerance: f64) -> bool {
    let third = 1. / 3.;
    [[0.5, 0.5, 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [third; 3]]
        .iter()
        .any(|w| {
            let at = (0..3).fold((0., 0.), |a, i| {
                (a.0 + w[i] * uv[i].0, a.1 + w[i] * uv[i].1)
            });
            let l = (0..3).fold([0.; 3], |a, i| add(a, mul(pts[i], w[i])));
            norm(sub(s.point(at.0, at.1), l)) > tolerance
        })
}

/// Bound on the distance between the surface and the triangle over the triangle alone (the
/// rectangle bound also covers its corners outside the triangle, which a thin triangle cannot
/// meet). With Ls the affine map through the surface's own corner points, |L - Ls| is largest
/// at a corner.
fn triangle_bound(s: &Surface, rect: [(f64, f64); 2], uv: [Uv; 3], pts: [P3; 3]) -> Option<f64> {
    let corners = uv.map(|p| s.point(p.0, p.1));
    let residual = (0..3)
        .map(|i| norm(sub(corners[i], pts[i])))
        .fold(0., f64::max);
    Some(residual + taylor(s, rect, uv, corners)?)
}

/// Bound on |S - Ls| over a triangle or segment whose parameter extents are `rect`, Ls the
/// affine map taking `uv` to `values` (surface points at least at the corners that end the
/// triangle or segment). E = A - W·Ls vanishes there, so a barycentric Taylor bound limits it
/// by its second derivatives. Needs a C1 surface over `rect`.
fn taylor(s: &Surface, rect: [(f64, f64); 2], uv: [Uv; 3], values: [P3; 3]) -> Option<f64> {
    if !smooth_over(s, rect, 1) {
        return None;
    }
    let (mut m, mut w_min) = ([0f64; 3], f64::INFINITY);
    for p in patches_over(s, rect, Some(widths(rect)))? {
        let (d, w) = p.second_derivatives(uv, values)?;
        for k in 0..3 {
            m[k] = m[k].max(d[k]);
        }
        w_min = w_min.min(w);
    }
    let bound = interpolation_error(m, rect) / w_min;
    bound.is_finite().then_some(bound)
}

/// Largest error of linear interpolation at a triangle's corners, over a triangle whose
/// parameter extents are `rect`, of a C1 function with second derivatives bounded by `m`
/// (along u, across u and v, along v).
fn interpolation_error(m: [f64; 3], [(u0, u1), (v0, v1)]: [(f64, f64); 2]) -> f64 {
    let (du, dv) = (u1 - u0, v1 - v0);
    m[0] * du * du / 8. + m[1] * du * dv / 4. + m[2] * dv * dv / 8.
}

/// Lines of constant u and of constant v inside `rect` where the spline surface is less than
/// C^`order`.
fn surface_breaks(s: &Surface, rect: [(f64, f64); 2], order: usize) -> Option<[Vec<f64>; 2]> {
    match s {
        Surface::Nurbs(n) => Some([
            breaks(
                &n.ku,
                n.du,
                n.ctrl.len(),
                n.closed_along_u(),
                rect[0],
                order,
            )?,
            breaks(
                &n.kv,
                n.dv,
                n.ctrl[0].len(),
                n.closed_along_v(),
                rect[1],
                order,
            )?,
        ]),
        Surface::Extrusion {
            curve: Curve::Nurbs(n),
            ..
        } => Some([
            breaks(&n.knots, n.degree, n.ctrl.len(), n.closed(), rect[0], order)?,
            Vec::new(),
        ]),
        Surface::Revolution {
            curve: Curve::Nurbs(n),
            ..
        } => Some([
            Vec::new(),
            breaks(&n.knots, n.degree, n.ctrl.len(), n.closed(), rect[1], order)?,
        ]),
        _ => None,
    }
}

fn smooth_over(s: &Surface, rect: [(f64, f64); 2], order: usize) -> bool {
    surface_breaks(s, rect, order).is_some_and(|b| b.iter().all(Vec::is_empty))
}

/// The triangle cut along the lines where the surface is less than C^`order`, as triangles,
/// each with the values of the affine map through `uv` and `pts` at its corners.
fn pieces(s: &Surface, uv: [Uv; 3], pts: [P3; 3], order: usize) -> Option<Vec<([Uv; 3], [P3; 3])>> {
    let lines = surface_breaks(s, rectangle(uv), order)?;
    if lines.iter().all(Vec::is_empty) {
        return Some(vec![(uv, pts)]);
    }
    let mut polygons = vec![uv.to_vec()];
    for (axis, cuts) in lines.iter().enumerate() {
        for &c in cuts {
            polygons = polygons
                .into_iter()
                .flat_map(|poly| cut(&poly, axis, c))
                .collect();
        }
    }
    let mut out = Vec::new();
    for poly in polygons {
        for i in 1..poly.len() - 1 {
            let piece = [poly[0], poly[i], poly[i + 1]];
            let [a, b, c] = piece.map(|x| affine(uv, pts, x));
            out.push((piece, [a?, b?, c?]));
        }
    }
    Some(out)
}

/// The two parts of a convex polygon on either side of the line where coordinate `axis` is `c`
/// (those with at least three corners); the corners on the line get `c` exactly.
fn cut(poly: &[Uv], axis: usize, c: f64) -> Vec<Vec<Uv>> {
    let at = |p: Uv| if axis == 0 { p.0 } else { p.1 };
    let mut sides = [Vec::new(), Vec::new()];
    for i in 0..poly.len() {
        let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
        let (a, b) = (at(p) - c, at(q) - c);
        if a <= 0. {
            sides[0].push(p);
        }
        if a >= 0. {
            sides[1].push(p);
        }
        if a * b < 0. {
            let t = a / (a - b);
            let x = (p.0 + (q.0 - p.0) * t, p.1 + (q.1 - p.1) * t);
            let x = if axis == 0 { (c, x.1) } else { (x.0, c) };
            sides[0].push(x);
            sides[1].push(x);
        }
    }
    sides.into_iter().filter(|s| s.len() >= 3).collect()
}

/// Bound on the distance between an offset surface and the triangle, the largest over its
/// pieces where the basis is C2.
fn offset_bound(
    basis: &Surface,
    distance: f64,
    uv: [Uv; 3],
    pts: [P3; 3],
    cache: &NormalBounds,
) -> Option<f64> {
    pieces(basis, uv, pts, 2)?
        .into_iter()
        .try_fold(0f64, |worst, (uv, pts)| {
            Some(worst.max(smooth_offset_bound(basis, distance, uv, pts, cache)?))
        })
}

/// `offset_bound` over a triangle where the basis is C2. S = B + d·N, and the affine maps
/// through the corners split the same way: the basis's Taylor bound plus |d| times that of its
/// unit normal N, whose derivatives are bounded over the cells around the triangle.
fn smooth_offset_bound(
    basis: &Surface,
    distance: f64,
    uv: [Uv; 3],
    pts: [P3; 3],
    cache: &NormalBounds,
) -> Option<f64> {
    let rect = rectangle(uv);
    if !smooth_over(basis, rect, 2) {
        return None;
    }
    if let Surface::Revolution {
        curve: Curve::Nurbs(n),
        axis_origin,
        axis_dir,
    } = basis
    {
        let normal = |v: f64| basis.normal(0., v);
        return revolution_offset_bound(n, *axis_origin, *axis_dir, distance, normal, uv, pts);
    }
    let corners = uv.map(|p| basis.point(p.0, p.1));
    let residual = (0..3)
        .map(|i| {
            let n = basis.normal(uv[i].0, uv[i].1);
            norm(sub(add(corners[i], mul(n, distance)), pts[i]))
        })
        .fold(0., f64::max);
    let base = residual + taylor(basis, rect, uv, corners)?;
    let m = cell_normal_bounds(basis, rect, cache)?;
    let bound = base + distance.abs() * interpolation_error(m, rect);
    bound.is_finite().then_some(bound)
}

/// Normal bounds of the cells met while meshing one face, by surface, knot spans, levels and
/// cells: neighbouring triangles share them.
#[derive(Default)]
pub(crate) struct NormalBounds(RefCell<HashMap<CellId, Option<[f64; 3]>>>);

/// A surface's address, then its cell's knot spans, levels and indices along u and v.
type CellId = (usize, [usize; 6]);

/// Times a knot span is halved at most for a cell.
const MAX_CELL_LEVEL: usize = 16;
/// Widths of the range a cell may span: wider cells serve more triangles.
const CELL_SHARE: f64 = 4.;

#[derive(Clone, Copy)]
struct Cell {
    span: usize,
    level: usize,
    index: usize,
    range: (f64, f64),
}

/// Normal bounds of the basis over cells covering `rect`, each a knot span halved until no wider
/// than `rect`: about as tight as over the triangle's own neighbourhood. An extrusion's normal
/// does not change along v, so its cells span v from 0 to 1.
fn cell_normal_bounds(
    basis: &Surface,
    rect: [(f64, f64); 2],
    cache: &NormalBounds,
) -> Option<[f64; 3]> {
    // One share of both ranges for both directions: a triangle thin along one gets the cells
    // of its neighbours, not cells of its own.
    let share = |k: &[f64], p: usize, n: usize, (lo, hi): (f64, f64)| (hi - lo) / (k[n] - k[p]);
    let (us, vs) = match basis {
        Surface::Nurbs(n) => {
            let (cu, cv) = (n.ctrl.len(), n.ctrl[0].len());
            let f =
                CELL_SHARE * share(&n.ku, n.du, cu, rect[0]).max(share(&n.kv, n.dv, cv, rect[1]));
            (
                cells(&n.ku, n.du, cu, n.closed_along_u(), rect[0], f)?,
                cells(&n.kv, n.dv, cv, n.closed_along_v(), rect[1], f)?,
            )
        }
        Surface::Extrusion {
            curve: Curve::Nurbs(n),
            ..
        } => (
            cells(
                &n.knots,
                n.degree,
                n.ctrl.len(),
                n.closed(),
                rect[0],
                CELL_SHARE * share(&n.knots, n.degree, n.ctrl.len(), rect[0]),
            )?,
            vec![Cell {
                span: 0,
                level: 0,
                index: 0,
                range: (0., 1.),
            }],
        ),
        _ => return None,
    };
    let surface = std::ptr::from_ref(basis) as usize;
    let mut m = [0f64; 3];
    for u in &us {
        for v in &vs {
            let id = (
                surface,
                [u.span, v.span, u.level, v.level, u.index, v.index],
            );
            let known = cache.0.borrow().get(&id).copied();
            let d = match known {
                Some(d) => d,
                None => {
                    let d = match patches(basis, [u.range, v.range])?.as_slice() {
                        [p] => p.normal_second_derivatives(),
                        _ => None,
                    };
                    cache.0.borrow_mut().insert(id, d);
                    d
                }
            }?;
            for k in 0..3 {
                m[k] = m[k].max(d[k]);
            }
        }
    }
    Some(m)
}

/// Cells covering `(lo, hi)`: in each knot span, the span halved until no wider than `share` of
/// the spline's whole range.
/// `None` for a range past the spline's ends (a cell there would extend it a whole cell) or
/// past `MAX_BREAKS` cells.
fn cells(
    k: &[f64],
    p: usize,
    n: usize,
    closed: bool,
    (lo, hi): (f64, f64),
    share: f64,
) -> Option<Vec<Cell>> {
    let mut out = Vec::new();
    for s in sections(k, p, n, closed, lo, hi)? {
        let (a, b) = (k[s.span], k[s.span + 1]);
        // A range a rounding error past an end (a face boundary on it) counts as at the end.
        let hair = 1e-9 * (b - a);
        if s.a < a - hair || s.b > b + hair || b <= a {
            return None;
        }
        let mut level = 0;
        while level < MAX_CELL_LEVEL && (b - a) / (1usize << level) as f64 > share * (k[n] - k[p]) {
            level += 1;
        }
        let parts = 1usize << level;
        let at = |i: usize| {
            if i == parts {
                b
            } else {
                a + (b - a) * i as f64 / parts as f64
            }
        };
        let w = (b - a) / parts as f64;
        let first = (((s.a - a) / w).floor().max(0.) as usize).min(parts - 1);
        let last = (((s.b - a) / w).ceil() as usize).clamp(first + 1, parts);
        let offset = s.lo - s.a;
        out.extend((first..last).map(|index| Cell {
            span: s.span,
            level,
            index,
            range: (offset + at(index), offset + at(index + 1)),
        }));
        if out.len() > MAX_BREAKS {
            return None;
        }
    }
    Some(out)
}

/// Most surface points sampled along one side of a thin triangle.
const MAX_FAN_SAMPLES: usize = 4096;
/// Most surface points sampled over a triangle by rows, in one of at most `ROW_PASSES`.
const MAX_ROW_SAMPLES: usize = 1 << 18;
const ROW_PASSES: usize = 3;

/// Whether a thin triangle lies within `tolerance` of the surface, measured to the nearest
/// surface point (a triangle thin in parameter space, or one along a boundary chord, follows the
/// point at the same parameters badly). The distance to the surface changes no faster than the
/// point moves, and the distance to any surface point bounds it, however well a search for the
/// nearest one converged.
pub(crate) fn fan_within(s: &Surface, uv: [Uv; 3], pts: [P3; 3], tolerance: f64) -> bool {
    sides_within(s, uv, pts, tolerance) || rows_within(s, uv, pts, tolerance)
}

/// `fan_within` by the sides at the corner opposite the longest side: split at its foot, every
/// point lies within `fan` of those sides (the distance to a segment is convex), so within `fan`
/// plus the bound of a side.
fn sides_within(s: &Surface, uv: [Uv; 3], pts: [P3; 3], tolerance: f64) -> bool {
    let side_length = |i: usize| norm(sub(pts[(i + 1) % 3], pts[i]));
    let long = (0..3)
        .max_by(|&i, &j| side_length(i).total_cmp(&side_length(j)))
        .unwrap();
    let (a, c, b) = (long, (long + 1) % 3, (long + 2) % 3);
    let base = sub(pts[c], pts[a]);
    let t = dot(sub(pts[b], pts[a]), base) / dot(base, base);
    if !t.is_finite() {
        return false;
    }
    let foot = add(pts[a], mul(base, t.clamp(0., 1.)));
    let fan = segment_distance(foot, pts[a], pts[b]).max(segment_distance(foot, pts[b], pts[c]));
    let budget = tolerance - fan;
    if budget.is_nan() || budget <= 0. {
        return false;
    }
    [(a, b), (b, c)]
        .into_iter()
        .all(|(x, y)| side_within(s, [uv[x], uv[y]], [pts[x], pts[y]], budget))
}

/// `fan_within` by rows of points parallel to the longest side, the lowest on it, each measured
/// to a surface point sought from its parameters. Every point of the triangle lies at most h/m
/// above a row (h the height, m the rows: the triangle narrows upwards) and Δ/2 along it from a
/// point (Δ their spacing), so within r = √((h/m)² + (Δ/2)²) of one: within its distance plus r.
fn rows_within(s: &Surface, uv: [Uv; 3], pts: [P3; 3], tolerance: f64) -> bool {
    // Rows for half the tolerance first; where the surface is farther, closer ones for the rest.
    let mut r = tolerance / 2.;
    for _ in 0..ROW_PASSES {
        let Some(worst) = rows_distance(s, uv, pts, r, tolerance) else {
            return false;
        };
        if worst + r <= tolerance {
            return true;
        }
        r = (tolerance - worst) / 2.;
    }
    false
}

/// Largest distance of the row points within `r` of every point of the triangle; `None` for
/// one at `tolerance` or farther, or past `MAX_ROW_SAMPLES` of them.
fn rows_distance(s: &Surface, uv: [Uv; 3], pts: [P3; 3], r: f64, tolerance: f64) -> Option<f64> {
    let side_length = |i: usize| norm(sub(pts[(i + 1) % 3], pts[i]));
    let long = (0..3)
        .max_by(|&i, &j| side_length(i).total_cmp(&side_length(j)))
        .unwrap();
    let (b, c, a) = (long, (long + 1) % 3, (long + 2) % 3);
    // h/m and Δ/2 at most r/√2 each.
    let rows = (segment_distance(pts[a], pts[b], pts[c]) * std::f64::consts::SQRT_2 / r)
        .ceil()
        .max(1.);
    let per_row = (side_length(long) / (r * std::f64::consts::SQRT_2))
        .ceil()
        .max(1.);
    let count = rows * (per_row + 1.);
    if !count.is_finite() || count > MAX_ROW_SAMPLES as f64 {
        return None;
    }
    let (m, n) = (rows as usize, per_row as usize);
    let mut worst: f64 = 0.;
    for (j, k) in (0..m).flat_map(|j| (0..=n).map(move |k| (j, k))) {
        let (f, g) = (j as f64 / m as f64, k as f64 / n as f64);
        let w = [((1. - f) * (1. - g), b), ((1. - f) * g, c), (f, a)];
        let point = w.iter().fold([0.; 3], |q, &(x, i)| add(q, mul(pts[i], x)));
        let hint = w.iter().fold((0., 0.), |q, &(x, i)| {
            (q.0 + x * uv[i].0, q.1 + x * uv[i].1)
        });
        let at = |(u, v): Uv| norm(sub(s.point(u, v), point));
        let d = at(s.approach(point, hint)).min(at(hint));
        if d.is_nan() || d >= tolerance {
            return None;
        }
        worst = worst.max(d);
    }
    Some(worst)
}

/// Whether every point of a mesh side lies within `budget` of the surface: by the Taylor
/// bound along its parameter segment, or else by samples spaced `step` (each within
/// `budget - step / 2`).
fn side_within(s: &Surface, uv: [Uv; 2], pts: [P3; 2], budget: f64) -> bool {
    let [p, q] = uv;
    // Any third point gives an affine extension; this one keeps it well conditioned.
    let tri = [p, q, (p.0 - (q.1 - p.1), p.1 + (q.0 - p.0))];
    let values = tri.map(|t| s.point(t.0, t.1));
    let residual = norm(sub(values[0], pts[0])).max(norm(sub(values[1], pts[1])));
    if taylor(s, rectangle([p, q, q]), tri, values).is_some_and(|d| residual + d <= budget) {
        return true;
    }
    let step = budget / 2.;
    let n = (norm(sub(pts[1], pts[0])) / step).ceil().max(1.);
    if !n.is_finite() || n > MAX_FAN_SAMPLES as f64 {
        return false;
    }
    (0..=n as usize).all(|k| {
        let f = k as f64 / n;
        let point = add(pts[0], mul(sub(pts[1], pts[0]), f));
        let hint = (p.0 + (q.0 - p.0) * f, p.1 + (q.1 - p.1) * f);
        // A local search is enough: any surface point bounds the distance.
        let at = |(u, v): Uv| norm(sub(s.point(u, v), point));
        let d = at(s.approach(point, hint)).min(at(hint));
        d.is_finite() && d + step / 2. <= budget
    })
}

/// Most lines a triangle is cut along to keep each piece smooth.
const MAX_BREAKS: usize = 32;

/// Parameters strictly inside `(lo, hi)` where the spline is less than C^`order`: its knots of
/// too high a multiplicity and a closed spline's seams. `None` past `MAX_BREAKS`.
fn breaks(
    k: &[f64],
    p: usize,
    n: usize,
    closed: bool,
    (lo, hi): (f64, f64),
    order: usize,
) -> Option<Vec<f64>> {
    let (start, end) = (k[p], k[n]);
    let period = end - start;
    let mut base = Vec::new();
    let mut i = 0;
    while i < k.len() {
        let mut j = i;
        while j < k.len() && k[j] == k[i] {
            j += 1;
        }
        // A knot of multiplicity m leaves the spline C^(p - m).
        if j - i + order > p {
            base.push(k[i]);
        }
        i = j;
    }
    if closed {
        base.push(start);
    }
    let mut out = Vec::new();
    for b in base {
        if !closed {
            if b > lo && b < hi {
                out.push(b);
            }
            continue;
        }
        let mut x = b + ((lo - b) / period).ceil() * period;
        while x < hi {
            if x > lo {
                out.push(x);
            }
            if out.len() > MAX_BREAKS {
                return None;
            }
            x += period;
        }
    }
    out.sort_by(f64::total_cmp);
    out.dedup();
    (out.len() <= MAX_BREAKS).then_some(out)
}

fn revolution_hull_within(
    n: &Nurbs1,
    origin: P3,
    axis: P3,
    uv: [Uv; 3],
    pts: [P3; 3],
    tolerance: f64,
) -> bool {
    let edge = (0..3)
        .max_by(|&i, &j| {
            norm(sub(pts[i], pts[(i + 1) % 3])).total_cmp(&norm(sub(pts[j], pts[(j + 1) % 3])))
        })
        .unwrap();
    let (a, b) = (pts[edge], pts[(edge + 1) % 3]);
    let height = segment_distance(pts[(edge + 2) % 3], a, b);
    if height >= tolerance {
        return false;
    }
    let [(u0, u1), (v0, v1)] = rectangle(uv);
    if u1 - u0 > std::f64::consts::FRAC_PI_2 {
        return false;
    }
    let rotate = |p: P3, u: f64| {
        let (sn, cs) = u.sin_cos();
        add(
            add(mul(p, cs), mul(crate::model::cross(axis, p), sn)),
            mul(axis, dot(axis, p) * (1. - cs)),
        )
    };
    let residual = (0..3)
        .map(|i| {
            norm(sub(
                add(origin, rotate(sub(n.derivs(uv[i].1).0, origin), uv[i].0)),
                pts[i],
            ))
        })
        .fold(0., f64::max);
    let Some(sections) = sections(&n.knots, n.degree, n.ctrl.len(), n.closed(), v0, v1) else {
        return false;
    };
    let cosine = ((u1 - u0) / 2.).cos();
    let mut deviation: f64 = 0.;
    for s in sections {
        for h in extract(&n.ctrl, n.degree, &n.knots, s) {
            if !valid(h) {
                return false;
            }
            let q = sub(cart(h), origin);
            let parallel = mul(axis, dot(axis, q));
            let radial = sub(q, parallel);
            for (u, scale) in [(u0, 1.), ((u0 + u1) / 2., 1. / cosine), (u1, 1.)] {
                let p = add(origin, add(parallel, mul(rotate(radial, u), scale)));
                let d = segment_distance(p, a, b);
                if !d.is_finite() {
                    return false;
                }
                deviation = deviation.max(d);
            }
        }
    }
    // The rational control hull bounds the whole patch about this chord. A surface path
    // joins its endpoints, so the chord and the thin triangle are covered in return.
    deviation + height + residual <= tolerance
}

fn revolution_bound(n: &Nurbs1, origin: P3, axis: P3, uv: [Uv; 3], pts: [P3; 3]) -> Option<f64> {
    let [_, (v0, v1)] = rectangle(uv);
    let ends = [n.derivs(v0).0, n.derivs(v1).0];
    rotated_line_bound(origin, axis, uv, pts, ends, profile_error(n, v0, v1, ends)?)
}

/// Largest distance between the profile over `[v0, v1]` and the line through `ends`, its points
/// at `v0` and `v1`.
fn profile_error(n: &Nurbs1, v0: f64, v1: f64, [c0, c1]: [P3; 2]) -> Option<f64> {
    if v1 == v0 {
        return None;
    }
    let line = |v: f64| add(c0, mul(sub(c1, c0), (v - v0) / (v1 - v0)));
    let mut error: f64 = 0.;
    for section in sections(&n.knots, n.degree, n.ctrl.len(), n.closed(), v0, v1)? {
        let c = extract(&n.ctrl, n.degree, &n.knots, section);
        if !c.iter().copied().all(valid) {
            return None;
        }
        let endpoints = [line(section.lo), line(section.hi)];
        let p = c.len() - 1;
        for i in 0..=p + 1 {
            let (mut num, mut w) = ([0.; 3], 0.);
            for (a, &value) in endpoints.iter().enumerate() {
                if i < a || i - a > p {
                    continue;
                }
                let f = (if a == 0 { p + 1 - i } else { i }) as f64 / (p + 1) as f64;
                let h = c[i - a];
                num = add(num, mul(sub(cart(h), value), f * h[3]));
                w += f * h[3];
            }
            let d = norm(num) / w;
            if !d.is_finite() {
                return None;
            }
            error = error.max(d);
        }
    }
    Some(error)
}

/// Bound on the distance between the triangle and the surface swept by turning a profile about
/// the axis, given the line through `ends` (its points at the rectangle's v0 and v1) and the
/// profile's largest distance `error` from it.
fn rotated_line_bound(
    origin: P3,
    axis: P3,
    uv: [Uv; 3],
    pts: [P3; 3],
    [c0, c1]: [P3; 2],
    error: f64,
) -> Option<f64> {
    let [(u0, u1), (v0, v1)] = rectangle(uv);
    if v1 == v0 {
        return None;
    }
    let line = |v: f64| add(c0, mul(sub(c1, c0), (v - v0) / (v1 - v0)));
    let radial = |p: P3| {
        let d = sub(p, origin);
        norm(sub(d, mul(axis, dot(d, axis))))
    };
    let radius = radial(c0).max(radial(c1));
    let du = u1 - u0;
    let rotate = |p: P3, u: f64| {
        let q = sub(p, origin);
        let (sn, cs) = u.sin_cos();
        add(
            origin,
            add(
                add(mul(q, cs), mul(crate::model::cross(axis, q), sn)),
                mul(axis, dot(axis, q) * (1. - cs)),
            ),
        )
    };
    let residual = (0..3)
        .map(|i| norm(sub(rotate(line(uv[i].1), uv[i].0), pts[i])))
        .fold(0., f64::max);
    // Linear interpolation of a rotated line: |S_uu| <= radius, and S_uv is the line's direction
    // turned about the axis, |S_uv| = |(c1-c0) across the axis|/dv.
    let across = radial(add(origin, sub(c1, c0)));
    let bound = error + residual + radius * du * du / 8. + across * du / 4.;
    bound.is_finite().then_some(bound)
}

/// `offset_bound` over a surface of revolution whose profile lies in a plane through the axis,
/// on one side of it (with the normal `normal` at the profile's own position). Its offset is
/// the revolution of the profile moved by `distance` along the profile's normal n = m × T, m
/// the plane's normal: so the profile's error plus |d| times n's linear interpolation error,
/// |n''| = |T''| ≤ |Q''|/r + 2|Q'|²/r² for Q along T and r ≤ |Q|.
fn revolution_offset_bound(
    n: &Nurbs1,
    origin: P3,
    axis: P3,
    distance: f64,
    normal: impl Fn(f64) -> P3,
    uv: [Uv; 3],
    pts: [P3; 3],
) -> Option<f64> {
    let [_, (v0, v1)] = rectangle(uv);
    if v1 <= v0 {
        return None;
    }
    let pieces: Vec<(Vec<H>, f64)> = widened(
        sections(&n.knots, n.degree, n.ctrl.len(), n.closed(), v0, v1)?,
        &n.knots,
        Some(v1 - v0),
    )
    .into_iter()
    .map(|s| (extract(&n.ctrl, n.degree, &n.knots, s), s.hi - s.lo))
    .collect();
    let controls: Vec<P3> = pieces
        .iter()
        .flat_map(|(c, _)| c.iter().map(|&h| sub(cart(h), origin)))
        .collect();
    let across = |q: P3| sub(q, mul(axis, dot(q, axis)));
    let far = controls
        .iter()
        .copied()
        .max_by(|&a, &b| norm(across(a)).total_cmp(&norm(across(b))))?;
    let out = unit(across(far))?;
    let plane = crate::model::cross(axis, out);
    // By the convex hull, the profile stays in the plane on the side of `out`.
    if !controls
        .iter()
        .all(|&q| dot(q, plane).abs() <= 1e-12 * (1. + norm(q)) && dot(q, out) > 0.)
    {
        return None;
    }
    let tangent = |d: H, h: H| sub(mul(xyz(d), h[3]), mul(xyz(h), d[3]));
    let mut curvature: f64 = 0.;
    for (c, width) in &pieces {
        if c.len() < 2 || !c.iter().copied().all(valid) {
            return None;
        }
        let net: Vec<Vec<H>> = c.iter().map(|&h| vec![h]).collect();
        let q = product(&diff(&net, true), &net, tangent);
        let r = least_norm(&q)?;
        let q1 = diff(&q, true);
        let (a, b) = (largest(&q1), largest(&diff(&q1, true)));
        curvature = curvature.max((b / r + 2. * a * a / (r * r)) / (width * width));
    }
    let on = [n.derivs(v0).0, n.derivs(v1).0];
    let ends = [
        add(on[0], mul(normal(v0), distance)),
        add(on[1], mul(normal(v1), distance)),
    ];
    let error =
        profile_error(n, v0, v1, on)? + distance.abs() * curvature * (v1 - v0) * (v1 - v0) / 8.;
    rotated_line_bound(origin, axis, uv, pts, ends, error)
}

#[derive(Clone, Copy)]
struct Interval(f64, f64);

impl Interval {
    fn minus(self, b: Self) -> Self {
        Self(self.0 - b.1, self.1 - b.0)
    }
    fn times(self, b: Self) -> Self {
        let values = [self.0 * b.0, self.0 * b.1, self.1 * b.0, self.1 * b.1];
        Self(
            values.into_iter().fold(f64::INFINITY, f64::min),
            values.into_iter().fold(f64::NEG_INFINITY, f64::max),
        )
    }
}

fn ranges(c: impl Iterator<Item = H>) -> [Interval; 4] {
    c.fold(
        [Interval(f64::INFINITY, f64::NEG_INFINITY); 4],
        |mut a, h| {
            for i in 0..4 {
                a[i].0 = a[i].0.min(h[i]);
                a[i].1 = a[i].1.max(h[i]);
            }
            a
        },
    )
}

fn cross_range(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        a[j].times(b[k]).minus(a[k].times(b[j]))
    })
}

fn cone(a: [Interval; 3]) -> Option<(P3, f64)> {
    let center = std::array::from_fn(|i| 0.5 * a[i].0 + 0.5 * a[i].1);
    let radius = norm(std::array::from_fn(|i| 0.5 * (a[i].1 - a[i].0)));
    let size = norm(center);
    if !radius.is_finite() || !size.is_finite() || radius >= size {
        return None;
    }
    Some((mul(center, 1. / size), (radius / size).clamp(0., 1.).asin()))
}

fn tangent_ranges(c: &[Vec<H>], along_u: bool) -> [Interval; 3] {
    let (p, q) = (c.len() - 1, c[0].len() - 1);
    let origin = cart(c[0][0]);
    let shifted: Vec<Vec<H>> = c
        .iter()
        .map(|r| {
            r.iter()
                .map(|&h| {
                    let v = mul(sub(cart(h), origin), h[3]);
                    [v[0], v[1], v[2], h[3]]
                })
                .collect()
        })
        .collect();
    let a = ranges(shifted.iter().flatten().copied());
    let mut derivative = Vec::new();
    for i in 0..=p - usize::from(along_u) {
        for j in 0..=q - usize::from(!along_u) {
            let (x, y) = (
                shifted[i][j],
                shifted[i + usize::from(along_u)][j + usize::from(!along_u)],
            );
            derivative.push(std::array::from_fn(|k| {
                (if along_u { p } else { q }) as f64 * (y[k] - x[k])
            }));
        }
    }
    let d = ranges(derivative.into_iter());
    std::array::from_fn(|i| d[i].times(a[3]).minus(a[i].times(d[3])))
}

fn normal_cones(s: &Surface, mut rect: [(f64, f64); 2]) -> Option<Vec<(P3, f64)>> {
    if let Surface::Offset { basis, .. } = s {
        return normal_cones(basis, rect);
    }
    let bounds = s.param_bounds();
    for i in 0..2 {
        let width = if bounds[i].0.is_finite() && bounds[i].1.is_finite() {
            bounds[i].1 - bounds[i].0
        } else {
            1.
        };
        let pad = (1e-5 * width - (rect[i].1 - rect[i].0)).max(0.) / 2.;
        rect[i] = (rect[i].0 - pad, rect[i].1 + pad);
    }
    if let Surface::Revolution {
        curve: Curve::Nurbs(n),
        axis_origin,
        axis_dir,
    } = s
    {
        let mut out = Vec::new();
        for section in sections(
            &n.knots,
            n.degree,
            n.ctrl.len(),
            n.closed(),
            rect[1].0,
            rect[1].1,
        )? {
            let controls = extract(&n.ctrl, n.degree, &n.knots, section);
            if !controls.iter().copied().all(valid) {
                return None;
            }
            let radial = ranges(controls.iter().map(|&h| {
                let p = sub(cart(h), *axis_origin);
                [p[0], p[1], p[2], 1.]
            }));
            let axis = axis_dir.map(|x| Interval(x, x));
            let tangent = tangent_ranges(
                &controls.into_iter().map(|h| vec![h]).collect::<Vec<_>>(),
                true,
            );
            let (dir, angle) = cone(cross_range(
                cross_range(axis, [radial[0], radial[1], radial[2]]),
                tangent,
            ))?;
            let u = 0.5 * rect[0].0 + 0.5 * rect[0].1;
            let (sn, cs) = u.sin_cos();
            let rotated = add(
                add(mul(dir, cs), mul(crate::model::cross(*axis_dir, dir), sn)),
                mul(*axis_dir, dot(*axis_dir, dir) * (1. - cs)),
            );
            out.push((
                rotated,
                (angle + 0.5 * (rect[0].1 - rect[0].0)).min(std::f64::consts::PI),
            ));
        }
        return Some(out);
    }
    patches(s, rect)?
        .into_iter()
        .map(|p| {
            if !p.c.iter().flatten().copied().all(valid) {
                return None;
            }
            cone(cross_range(
                tangent_ranges(&p.c, true),
                tangent_ranges(&p.c, false),
            ))
        })
        .collect()
}

fn offset_within(
    basis: &Surface,
    distance: f64,
    uv: [Uv; 3],
    pts: [P3; 3],
    o: &Options,
    cache: &NormalBounds,
) -> bool {
    if distance == 0. {
        return triangle_within(basis, uv, pts, o, cache) == Some(true);
    }
    if let Some(d) = offset_bound(basis, distance, uv, pts, cache) {
        return d <= o.linear_deflection;
    }
    // Without second derivatives (a crease, a degenerate normal), the normal cones bound
    // N - L_N to first order.
    let normals = uv.map(|p| basis.normal(p.0, p.1));
    let Some(cones) = normal_cones(basis, rectangle(uv)) else {
        return false;
    };
    let extra = cones
        .into_iter()
        .map(|(axis, angle)| {
            distance.abs()
                * (2. * (angle / 2.).sin()
                    + normals
                        .iter()
                        .map(|&n| norm(sub(n, axis)))
                        .fold(0., f64::max))
        })
        .fold(0., f64::max);
    if !extra.is_finite() || extra >= o.linear_deflection {
        return false;
    }
    let mut tighter = o.clone();
    tighter.linear_deflection -= extra;
    triangle_within(
        basis,
        uv,
        std::array::from_fn(|i| sub(pts[i], mul(normals[i], distance))),
        &tighter,
        cache,
    ) == Some(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Nurbs2;

    fn bezier_at(c: &[H], t: f64) -> P3 {
        let mut d = c.to_vec();
        for n in (1..d.len()).rev() {
            for j in 0..n {
                d[j] = std::array::from_fn(|k| (1. - t) * d[j][k] + t * d[j + 1][k]);
            }
        }
        cart(d[0])
    }

    #[test]
    fn sections_of_a_closed_spline_stay_in_their_spans_in_every_period() {
        let h = std::f64::consts::FRAC_PI_4;
        let k: Vec<f64> = [0., 0., 0., 0., 1., 2., 3., 4., 5., 5., 5., 5.]
            .map(|i| i * h)
            .to_vec();
        let period = 5. * h;
        for m in -2..=3 {
            for j in 0..10 {
                let lo = j as f64 * 0.37 * h + f64::from(m) * period;
                for hi in [lo + 0.9 * h, lo + 2.3 * h, lo + period] {
                    for s in sections(&k, 3, 8, true, lo, hi).unwrap() {
                        let tol = 1e-12 * period;
                        assert!(
                            s.a >= k[s.span] - tol && s.b <= k[s.span + 1] + tol,
                            "period {m}, ({lo}, {hi}): span {} gets [{}, {}]",
                            s.span,
                            s.a,
                            s.b
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn extraction_matches_unclamped_rational_splines_and_extensions() {
        for degree in [1, 2, 3, 7, 16] {
            let n = degree + 4;
            let spline = Nurbs1 {
                degree,
                knots: (0..n + degree + 1)
                    .map(|i| i as f64 - degree as f64)
                    .collect(),
                ctrl: (0..n)
                    .map(|i| {
                        let x = i as f64;
                        let w = 1. + x / 10.;
                        [x * w, x.sin() * w, x.cos() * w, w]
                    })
                    .collect(),
            };
            for piece in sections(&spline.knots, degree, n, false, -0.05, 4.05).unwrap() {
                let c = extract(&spline.ctrl, degree, &spline.knots, piece);
                for j in 0..=10 {
                    let t = j as f64 / 10.;
                    let expected = spline.derivs(piece.a + (piece.b - piece.a) * t).0;
                    assert!(
                        norm(sub(bezier_at(&c, t), expected)) < 1e-10,
                        "degree {degree}"
                    );
                }
            }
        }
    }

    fn surface() -> Surface {
        Surface::Nurbs(Nurbs2 {
            du: 3,
            dv: 2,
            ku: vec![0., 0., 0., 0., 0.5, 1., 1., 1., 1.],
            kv: vec![0., 0., 0., 1., 1., 1.],
            closed_u: false,
            closed_v: false,
            ctrl: (0..5)
                .map(|i| {
                    (0..3)
                        .map(|j| {
                            let w = 0.5 + (i + j) as f64 / 4.;
                            [i as f64 * w, j as f64 * w, (i * j) as f64 * w / 10., w]
                        })
                        .collect()
                })
                .collect(),
        })
    }

    #[test]
    fn rational_surface_bounds_include_dense_triangle_samples() {
        let s = surface();
        let uv = [(0.12, 0.1), (0.91, 0.2), (0.35, 0.95)];
        let pts = uv.map(|p| s.point(p.0, p.1));
        let pieces = patches(&s, rectangle(uv)).unwrap();
        assert_eq!(pieces.len(), 2);
        for patch in &pieces {
            let bound = patch.position_bound(uv, pts).unwrap();
            for i in 0..=40 {
                for j in 0..=40 - i {
                    let (a, b) = (i as f64 / 40., j as f64 / 40.);
                    let at = (
                        uv[0].0 + a * (uv[1].0 - uv[0].0) + b * (uv[2].0 - uv[0].0),
                        uv[0].1 + a * (uv[1].1 - uv[0].1) + b * (uv[2].1 - uv[0].1),
                    );
                    if at.0 < patch.u.0 || at.0 > patch.u.1 || at.1 < patch.v.0 || at.1 > patch.v.1
                    {
                        continue;
                    }
                    assert!(
                        norm(sub(s.point(at.0, at.1), affine(uv, pts, at).unwrap()))
                            <= bound + 1e-12
                    );
                }
            }
        }
        assert_eq!(
            triangle_within(&s, uv, pts, &Options::default(), &NormalBounds::default()),
            Some(false)
        );
    }

    #[test]
    fn a_hair_long_spline_piece_is_one_point() {
        // The last section of a closed edge, an ulp long: its control points differ by rounding.
        let c: Vec<H> = [[0., 0.], [-1., -1.2], [-2., -2.4], [-1., -3.6]]
            .iter()
            .map(|&[x, z]| [-7.25 + x * 1e-15, -39.34, z * 1e-15, 1.])
            .collect();
        assert!(curve_ok(&c, &Options::default()));
    }

    #[test]
    fn periodic_sampling_splits_at_the_seam_in_both_directions() {
        let n = Nurbs1 {
            degree: 1,
            knots: vec![0., 0., 1., 2., 3., 4., 4.],
            ctrl: vec![
                [0., 0., 0., 1.],
                [1., 0., 0., 1.],
                [1., 1., 0., 1.],
                [0., 1., 0., 1.],
                [0., 0., 0., 1.],
            ],
        };
        let o = Options::default();
        assert_eq!(sample_nurbs(&n, 3.5, 4.5, &o).unwrap(), vec![3.5, 4., 4.5]);
        assert_eq!(sample_nurbs(&n, 4.5, 3.5, &o).unwrap(), vec![4.5, 4., 3.5]);
    }

    #[test]
    fn periodic_surface_bounds_wrap_u_and_v() {
        let rows = vec![
            vec![[0., 0., 0., 1.], [0., 0., 1., 1.]],
            vec![[1., 0., 0., 1.], [1., 0., 1., 1.]],
            vec![[1., 1., 0., 1.], [1., 1., 1., 1.]],
            vec![[0., 1., 0., 1.], [0., 1., 1., 1.]],
            vec![[0., 0., 0., 1.], [0., 0., 1., 1.]],
        ];
        let mut n = Nurbs2 {
            du: 1,
            dv: 1,
            ku: vec![0., 0., 1., 2., 3., 4., 4.],
            kv: vec![0., 0., 1., 1.],
            ctrl: rows,
            closed_u: true,
            closed_v: false,
        };
        for swapped in [false, true] {
            if swapped {
                n.ctrl = (0..2)
                    .map(|j| n.ctrl.iter().map(|r| r[j]).collect())
                    .collect();
                std::mem::swap(&mut n.ku, &mut n.kv);
            }
            let uv = [(3.75, 0.1), (4.25, 0.1), (3.75, 0.9)]
                .map(|(u, v)| if swapped { (v, u) } else { (u, v) });
            let s = Surface::Nurbs(n.clone());
            let pts = uv.map(|p| s.point(p.0, p.1));
            let o = Options {
                linear_deflection: 1.,
                ..Options::default()
            };
            assert_eq!(
                triangle_within(&s, uv, pts, &o, &NormalBounds::default()),
                Some(true)
            );
            assert_eq!(
                triangle_within(&s, uv, pts, &Options::default(), &NormalBounds::default()),
                Some(false)
            );
        }
    }

    #[test]
    fn swept_and_offset_splines_have_bounds() {
        let curve = Curve::Nurbs(Nurbs1 {
            degree: 3,
            knots: vec![0., 0., 0., 0., 1., 1., 1., 1.],
            ctrl: vec![
                [2., 0., 0., 1.],
                [3., 0., 0.3, 1.],
                [1., 0., 0.7, 1.],
                [2., 0., 1., 1.],
            ],
        });
        for s in [
            Surface::Extrusion {
                curve: curve.clone(),
                dir: [0., 1., 0.],
            },
            Surface::Revolution {
                curve,
                axis_origin: [0.; 3],
                axis_dir: [0., 0., 1.],
            },
            Surface::Offset {
                basis: Box::new(surface()),
                distance: 0.,
            },
            Surface::Offset {
                basis: Box::new(surface()),
                distance: 0.1,
            },
        ] {
            let uv = [(0., 0.), (0.8, 0.), (0.2, 1.)];
            let pts = uv.map(|p| s.point(p.0, p.1));
            assert_eq!(
                triangle_within(&s, uv, pts, &Options::default(), &NormalBounds::default()),
                Some(false)
            );
        }
    }

    #[test]
    fn flat_offset_spline_can_be_certified() {
        let basis = Surface::Nurbs(Nurbs2 {
            du: 1,
            dv: 1,
            ku: vec![0., 0., 1., 1.],
            kv: vec![0., 0., 1., 1.],
            closed_u: false,
            closed_v: false,
            ctrl: vec![
                vec![[0., 0., 0., 1.], [0., 1., 0., 1.]],
                vec![[1., 0., 0., 1.], [1., 1., 0., 1.]],
            ],
        });
        let s = Surface::Offset {
            basis: Box::new(basis),
            distance: 1.5,
        };
        let uv = [(0., 0.), (1., 0.), (0., 1.)];
        let pts = uv.map(|p| s.point(p.0, p.1));
        assert_eq!(
            triangle_within(&s, uv, pts, &Options::default(), &NormalBounds::default()),
            Some(true)
        );
    }

    fn plane() -> Surface {
        Surface::Nurbs(Nurbs2 {
            du: 1,
            dv: 1,
            ku: vec![0., 0., 10., 10.],
            kv: vec![0., 0., 10., 10.],
            closed_u: false,
            closed_v: false,
            ctrl: vec![
                vec![[0., 0., 0., 1.], [0., 10., 0., 1.]],
                vec![[10., 0., 0., 1.], [10., 10., 0., 1.]],
            ],
        })
    }

    /// Largest distance between the surface and the triangle at the same parameters.
    fn sampled_error(s: &Surface, uv: [Uv; 3], pts: [P3; 3]) -> f64 {
        let mut worst: f64 = 0.;
        for i in 0..=60 {
            for j in 0..=60 - i {
                let (a, b) = (i as f64 / 60., j as f64 / 60.);
                let at = (
                    uv[0].0 + a * (uv[1].0 - uv[0].0) + b * (uv[2].0 - uv[0].0),
                    uv[0].1 + a * (uv[1].1 - uv[0].1) + b * (uv[2].1 - uv[0].1),
                );
                let l = add(
                    pts[0],
                    add(mul(sub(pts[1], pts[0]), a), mul(sub(pts[2], pts[0]), b)),
                );
                worst = worst.max(norm(sub(s.point(at.0, at.1), l)));
            }
        }
        worst
    }

    #[test]
    fn triangle_bound_holds_for_wide_and_thin_triangles() {
        let s = surface();
        for uv in [
            [(0.12, 0.1), (0.91, 0.2), (0.35, 0.95)],
            [(0.1, 0.5), (0.9, 0.52), (0.5, 0.511)],
            [(0.6, 0.1), (0.61, 0.9), (0.605, 0.4)],
        ] {
            // corners a little off the surface, as boundary points are
            let pts = uv.map(|p| add(s.point(p.0, p.1), [1e-4, -1e-4, 2e-4]));
            let bound = triangle_bound(&s, rectangle(uv), uv, pts).unwrap();
            assert!(sampled_error(&s, uv, pts) <= bound + 1e-12, "{uv:?}");
        }
    }

    #[test]
    fn offset_bound_holds_and_is_second_order() {
        // Only C1 at u = 0.5, where its curvature jumps: the bound is taken on either side.
        let crease = Surface::Nurbs(Nurbs2 {
            du: 3,
            dv: 2,
            ku: vec![0., 0., 0., 0., 0.5, 0.5, 1., 1., 1., 1.],
            kv: vec![0., 0., 0., 1., 1., 1.],
            closed_u: false,
            closed_v: false,
            ctrl: (0..6)
                .map(|i| {
                    (0..3)
                        .map(|j| {
                            let w = 0.5 + (i + j) as f64 / 4.;
                            let z = (i * j) as f64 / 10. + if i == 3 { 0.4 } else { 0. };
                            [i as f64 * w, j as f64 * w, z * w, w]
                        })
                        .collect()
                })
                .collect(),
        });
        let basis = surface();
        for (b, distance) in [(&basis, 0.3), (&basis, -0.05), (&basis, 5.), (&crease, 2.)] {
            let s = Surface::Offset {
                basis: Box::new(b.clone()),
                distance,
            };
            for uv in [
                [(0.12, 0.1), (0.91, 0.2), (0.35, 0.95)],
                [(0.1, 0.5), (0.9, 0.52), (0.5, 0.511)],
                [(0.3, 0.3), (0.32, 0.3), (0.3, 0.32)],
            ] {
                let pts = uv.map(|p| add(s.point(p.0, p.1), [1e-4, -1e-4, 2e-4]));
                let bound = offset_bound(b, distance, uv, pts, &NormalBounds::default()).unwrap();
                assert!(sampled_error(&s, uv, pts) <= bound + 1e-12, "{uv:?}");
            }
        }
        // Halving a triangle quarters the normal's share of the bound.
        let s = Surface::Offset {
            basis: Box::new(basis.clone()),
            distance: 1.,
        };
        let [big, small] = [0.1, 0.05].map(|h| {
            let uv = [(0.3, 0.3), (0.3 + h, 0.3), (0.3, 0.3 + h)];
            let on = |s: &Surface| uv.map(|p| s.point(p.0, p.1));
            offset_bound(&basis, 1., uv, on(&s), &NormalBounds::default()).unwrap()
                - offset_bound(&basis, 0., uv, on(&basis), &NormalBounds::default()).unwrap()
        });
        assert!(small < big / 3., "{big} {small}");
    }

    #[test]
    fn a_revolved_line_along_the_axis_bounds_long_thin_triangles() {
        // A cylinder as the revolution of a line: the line does not turn as it is revolved.
        let line = Nurbs1 {
            degree: 1,
            knots: vec![0., 0., 1., 1.],
            ctrl: vec![[0.5, 0., 0., 1.], [0.5, 0., 14., 1.]],
        };
        let s = Surface::Revolution {
            curve: Curve::Nurbs(line.clone()),
            axis_origin: [0.; 3],
            axis_dir: [0., 0., 1.],
        };
        let uv = [(0., 0.), (0.01, 0.5), (0., 1.)];
        let pts = uv.map(|p| s.point(p.0, p.1));
        let bound = revolution_bound(&line, [0.; 3], [0., 0., 1.], uv, pts).unwrap();
        assert!(bound <= Options::default().linear_deflection, "{bound}");
    }

    #[test]
    fn offset_of_a_revolution_has_a_second_order_bound() {
        let basis = Surface::Revolution {
            curve: Curve::Nurbs(Nurbs1 {
                degree: 3,
                knots: vec![0., 0., 0., 0., 1., 1., 1., 1.],
                ctrl: vec![
                    [2., 0., 0., 1.],
                    [3., 0., 0.3, 1.],
                    [1., 0., 0.7, 1.],
                    [2., 0., 1., 1.],
                ],
            }),
            axis_origin: [0.; 3],
            axis_dir: [0., 0., 1.],
        };
        for distance in [0.1, -0.1] {
            let s = Surface::Offset {
                basis: Box::new(basis.clone()),
                distance,
            };
            for uv in [
                [(0., 0.1), (0.5, 0.2), (0.2, 0.9)],
                [(1., 0.4), (1.05, 0.4), (1., 0.45)],
            ] {
                let pts = uv.map(|p| add(s.point(p.0, p.1), [1e-4, -1e-4, 2e-4]));
                let bound =
                    offset_bound(&basis, distance, uv, pts, &NormalBounds::default()).unwrap();
                assert!(sampled_error(&s, uv, pts) <= bound + 1e-12, "{uv:?}");
            }
        }
    }

    #[test]
    fn creases_and_seams_leave_the_triangle_bound() {
        let k = [0., 0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1., 1.];
        assert_eq!(breaks(&k, 3, 7, false, (0.4, 0.6), 1), Some(vec![0.5]));
        assert_eq!(breaks(&k, 3, 7, false, (0.1, 0.45), 1), Some(vec![]));
        let c1 = [0., 0., 0., 0., 0.5, 0.5, 1., 1., 1., 1.];
        assert_eq!(breaks(&c1, 3, 6, false, (0.4, 0.6), 1), Some(vec![]));
        assert_eq!(breaks(&c1, 3, 6, false, (0.4, 0.6), 2), Some(vec![0.5]));
        let periodic = [0., 0., 1., 2., 3., 4., 4.];
        assert_eq!(breaks(&periodic, 1, 5, true, (3.5, 4.5), 1), Some(vec![4.]));
    }

    #[test]
    fn thin_triangles_are_measured_to_the_nearest_surface_point() {
        let s = plane();
        let o = Options::default();
        // On the plane, but given parameters 0.01 off: only the nearest point is close.
        let pts = [[1., 1., 0.], [9., 1.001, 0.], [5., 1.0004, 0.]];
        let uv = pts.map(|p| (p[0] + 0.01, p[1]));
        assert_eq!(
            triangle_within(&s, uv, pts, &o, &NormalBounds::default()),
            Some(false)
        );
        assert!(fan_within(&s, uv, pts, o.linear_deflection));
        // Off the plane by more than the tolerance.
        let lifted = pts.map(|p| [p[0], p[1], 0.004]);
        assert!(!fan_within(&s, uv, lifted, o.linear_deflection));
        // Higher than the tolerance, so measured by rows of points instead of its fan.
        let strip = [[1., 1., 0.], [2., 1., 0.], [1.5, 1.005, 0.]];
        let strip_uv = strip.map(|p| (p[0] + 0.01, p[1]));
        assert!(!sides_within(&s, strip_uv, strip, o.linear_deflection));
        assert!(fan_within(&s, strip_uv, strip, o.linear_deflection));
        // Off the plane by more than half the tolerance: closer rows.
        let lifted = strip.map(|p| [p[0], p[1], 0.002]);
        assert!(fan_within(&s, strip_uv, lifted, o.linear_deflection));
        let lifted = strip.map(|p| [p[0], p[1], 0.0035]);
        assert!(!fan_within(&s, strip_uv, lifted, o.linear_deflection));
        // Too wide for either.
        let wide = [[1., 1., 0.], [9., 1., 0.], [5., 3., 0.]];
        assert!(!fan_within(
            &s,
            wide.map(|p| (p[0], p[1])),
            wide,
            o.linear_deflection
        ));
    }
}
