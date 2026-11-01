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

/// Non-zero basis functions and their first derivatives at `u` (A2.3 with n = 1).
fn basis(span: usize, u: f64, p: usize, k: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let (mut left, mut right) = (vec![0.0; p + 1], vec![0.0; p + 1]);
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
    let n: Vec<f64> = (0..=p).map(|j| ndu[j][p]).collect();
    let d = (0..=p)
        .map(|r| {
            if p == 0 {
                return 0.0;
            }
            let mut v = 0.0;
            if r >= 1 {
                v += ndu[r - 1][p - 1] / ndu[p][r - 1];
            }
            if r < p {
                v -= ndu[r][p - 1] / ndu[p][r];
            }
            let v = v * p as f64;
            if v.is_finite() { v } else { 0.0 }
        })
        .collect();
    (n, d)
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

    /// Point and first derivative.
    pub fn derivs(&self, t: f64) -> ([f64; 3], [f64; 3]) {
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

    /// Point and the partial derivatives along u and v.
    pub fn derivs(&self, u: f64, v: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let (nu, nv) = (self.ctrl.len() - 1, self.ctrl[0].len() - 1);
        let su = find_span(nu, self.du, u, &self.ku);
        let sv = find_span(nv, self.dv, v, &self.kv);
        let (bu, dbu) = basis(su, u, self.du, &self.ku);
        let (bv, dbv) = basis(sv, v, self.dv, &self.kv);
        let (mut a, mut au, mut av) = ([0.0; 4], [0.0; 4], [0.0; 4]);
        for (i, (bi, dbi)) in bu.iter().zip(&dbu).enumerate() {
            for (j, (bj, dbj)) in bv.iter().zip(&dbv).enumerate() {
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
