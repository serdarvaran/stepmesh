//! Rational B-splines (Piegl & Tiller, *The NURBS Book*, A2.1, A2.3).

/// Knot vector from distinct knots and their multiplicities.
pub fn expand_knots(multiplicities: &[usize], knots: &[f64]) -> Vec<f64> {
    multiplicities
        .iter()
        .zip(knots)
        .flat_map(|(&m, &k)| std::iter::repeat_n(k, m))
        .collect()
}

/// Curve with homogeneous control points `[w·x, w·y, w·z, w]`.
#[derive(Debug, Clone, PartialEq)]
pub struct Nurbs1 {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub ctrl: Vec<[f64; 4]>,
}

/// Surface; `ctrl[i][j]` is the point at u index `i`, v index `j`.
#[derive(Debug, Clone, PartialEq)]
pub struct Nurbs2 {
    pub du: usize,
    pub dv: usize,
    pub ku: Vec<f64>,
    pub kv: Vec<f64>,
    pub ctrl: Vec<Vec<[f64; 4]>>,
    pub closed_u: bool,
    pub closed_v: bool,
}

/// Knot span index containing `u` (A2.1); `n` is the last control point index.
fn find_span(n: usize, p: usize, u: f64, k: &[f64]) -> usize {
    if u >= k[n + 1] {
        // the end of the range belongs to the last non-empty span
        let mut s = n;
        while s > p && k[s] >= k[n + 1] {
            s -= 1;
        }
        return s;
    }
    if u <= k[p] {
        let mut s = p;
        while s < n && k[s + 1] <= k[p] {
            s += 1;
        }
        return s;
    }
    let (mut lo, mut hi) = (p, n + 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if u < k[mid] { hi = mid } else { lo = mid }
    }
    lo
}

/// Highest degree read; the basis is computed on the stack up to it.
pub const MAX_DEGREE: usize = 16;
const N: usize = MAX_DEGREE + 1;

/// Non-zero basis functions and their first derivatives at `u` (A2.3 with n = 1); entries
/// `0..=p` are meaningful.
fn basis(span: usize, u: f64, p: usize, k: &[f64]) -> ([f64; N], [f64; N]) {
    // Small work arrays for the usual degrees: zeroing the full table dominated the cost.
    match p {
        0..=3 => basis_in::<4>(span, u, p, k),
        4..=7 => basis_in::<8>(span, u, p, k),
        _ => basis_in::<N>(span, u, p, k),
    }
}

fn basis_in<const M: usize>(span: usize, u: f64, p: usize, k: &[f64]) -> ([f64; N], [f64; N]) {
    let mut ndu = [[0.0; M]; M];
    let (mut left, mut right) = ([0.0; M], [0.0; M]);
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - k[span + 1 - j];
        right[j] = k[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r] == 0.0 {
                0.0
            } else {
                ndu[r][j - 1] / ndu[j][r]
            };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let (mut n, mut d) = ([0.0; N], [0.0; N]);
    for r in 0..=p {
        n[r] = ndu[r][p];
        if p == 0 {
            continue;
        }
        let mut v = 0.0;
        if r >= 1 {
            v += ndu[r - 1][p - 1] / ndu[p][r - 1];
        }
        if r < p {
            v -= ndu[r][p - 1] / ndu[p][r];
        }
        let v = v * p as f64;
        d[r] = if v.is_finite() { v } else { 0.0 };
    }
    (n, d)
}

fn euclidean(c: [f64; 4]) -> [f64; 3] {
    [c[0] / c[3], c[1] / c[3], c[2] / c[3]]
}

/// Points this close (millimetres, plus a relative part) close a curve or surface on itself:
/// files leave seams apart by their modelling tolerance.
const CLOSING: f64 = 1e-5;

fn same_point(a: [f64; 4], b: [f64; 4]) -> bool {
    near(euclidean(a), euclidean(b))
}

fn near(p: [f64; 3], q: [f64; 3]) -> bool {
    let d: f64 = (0..3).map(|i| (p[i] - q[i]).powi(2)).sum();
    let s: f64 = p.iter().map(|x| x * x).sum();
    d.sqrt() <= CLOSING + 1e-9 * s.sqrt()
}

/// End knots repeated `degree + 1` times: the curve starts and ends at its end control points.
fn clamped(k: &[f64], p: usize) -> bool {
    k.len() > 2 * p + 1
        && k[..=p].iter().all(|&x| x == k[0])
        && k[k.len() - 1 - p..].iter().all(|&x| x == k[k.len() - 1])
}

/// `t` moved into `[t0, t1)` by whole periods.
fn wrap(t: f64, (t0, t1): (f64, f64)) -> f64 {
    if t < t0 || t > t1 {
        t0 + (t - t0).rem_euclid(t1 - t0)
    } else {
        t
    }
}

fn rational(a: [f64; 4], da: [f64; 4]) -> ([f64; 3], [f64; 3]) {
    let w = a[3];
    let p = [a[0] / w, a[1] / w, a[2] / w];
    let d = [
        (da[0] - da[3] * p[0]) / w,
        (da[1] - da[3] * p[1]) / w,
        (da[2] - da[3] * p[2]) / w,
    ];
    (p, d)
}

impl Nurbs1 {
    pub fn range(&self) -> (f64, f64) {
        (self.knots[self.degree], self.knots[self.ctrl.len()])
    }

    /// Whether the curve ends where it starts (then parameters outside the range wrap round).
    pub fn closed(&self) -> bool {
        let (first, last) = (self.ctrl[0], self.ctrl[self.ctrl.len() - 1]);
        if clamped(&self.knots, self.degree) {
            return same_point(first, last);
        }
        let (t0, t1) = self.range();
        near(self.at(t0).0, self.at(t1).0)
    }

    /// Point and first derivative.
    pub fn derivs(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        let (t0, t1) = self.range();
        if (t < t0 || t > t1) && self.closed() {
            return self.at(wrap(t, (t0, t1)));
        }
        self.at(t)
    }

    fn at(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        let n = self.ctrl.len() - 1;
        let s = find_span(n, self.degree, t, &self.knots);
        let (b, db) = basis(s, t, self.degree, &self.knots);
        let (mut a, mut da) = ([0.0; 4], [0.0; 4]);
        for j in 0..=self.degree {
            let c = self.ctrl[s - self.degree + j];
            for i in 0..4 {
                a[i] += b[j] * c[i];
                da[i] += db[j] * c[i];
            }
        }
        rational(a, da)
    }
}

impl Nurbs2 {
    pub fn range(&self) -> [(f64, f64); 2] {
        let (nu, nv) = (self.ctrl.len(), self.ctrl[0].len());
        [
            (self.ku[self.du], self.ku[nu]),
            (self.kv[self.dv], self.kv[nv]),
        ]
    }

    /// Whether the surface closes on itself along u (first and last rows of control points
    /// coincide on a clamped knot vector).
    pub fn closed_along_u(&self) -> bool {
        if clamped(&self.ku, self.du) {
            let last = &self.ctrl[self.ctrl.len() - 1];
            return self.ctrl[0]
                .iter()
                .zip(last)
                .all(|(&a, &b)| same_point(a, b));
        }
        let [(u0, u1), (v0, v1)] = self.range();
        (0..5).all(|i| {
            let v = v0 + (v1 - v0) * i as f64 / 4.0;
            near(self.at(u0, v).0, self.at(u1, v).0)
        })
    }

    pub fn closed_along_v(&self) -> bool {
        if clamped(&self.kv, self.dv) {
            return self
                .ctrl
                .iter()
                .all(|row| same_point(row[0], row[row.len() - 1]));
        }
        let [(u0, u1), (v0, v1)] = self.range();
        (0..5).all(|i| {
            let u = u0 + (u1 - u0) * i as f64 / 4.0;
            near(self.at(u, v0).0, self.at(u, v1).0)
        })
    }

    /// Point and the partial derivatives along u and v; parameters outside the range of a closed
    /// direction wrap round.
    pub fn derivs(&self, u: f64, v: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let [ru, rv] = self.range();
        let u = if (u < ru.0 || u > ru.1) && self.closed_along_u() {
            wrap(u, ru)
        } else {
            u
        };
        let v = if (v < rv.0 || v > rv.1) && self.closed_along_v() {
            wrap(v, rv)
        } else {
            v
        };
        self.at(u, v)
    }

    fn at(&self, u: f64, v: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let (nu, nv) = (self.ctrl.len() - 1, self.ctrl[0].len() - 1);
        let su = find_span(nu, self.du, u, &self.ku);
        let sv = find_span(nv, self.dv, v, &self.kv);
        let (bu, dbu) = basis(su, u, self.du, &self.ku);
        let (bv, dbv) = basis(sv, v, self.dv, &self.kv);
        let (mut a, mut au, mut av) = ([0.0; 4], [0.0; 4], [0.0; 4]);
        for (i, (bi, dbi)) in bu[..=self.du].iter().zip(&dbu).enumerate() {
            for (j, (bj, dbj)) in bv[..=self.dv].iter().zip(&dbv).enumerate() {
                let c = self.ctrl[su - self.du + i][sv - self.dv + j];
                for k in 0..4 {
                    a[k] += bi * bj * c[k];
                    au[k] += dbi * bj * c[k];
                    av[k] += bi * dbj * c[k];
                }
            }
        }
        let (p, du) = rational(a, au);
        let (_, dv) = rational(a, av);
        (p, du, dv)
    }
}
