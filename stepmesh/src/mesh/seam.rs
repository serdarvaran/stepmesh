//! Seams inside a face: where a loop going round a periodic surface is joined to a pole or to
//! the other loop of a band, and where the seams of a whole parameter domain go.

use std::collections::VecDeque;

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

/// Cells along each side of the grids a winding seam is sought on, coarsest first.
const ROUTE_GRIDS: [usize; 5] = [64, 128, 256, 512, 1024];
/// Points of a loop tried as the end of a winding seam, nearest first.
const ROUTE_ENDS: usize = 16;

type Segment = ((f64, f64), (f64, f64));

/// A seam winding between holes: the points of the two loops it joins and its corners from the
/// first one on, continuous in u; the last corner is the second point moved by whole periods.
pub(super) struct Route {
    pub(super) from: usize,
    pub(super) to: usize,
    pub(super) corners: Vec<(f64, f64)>,
}

/// A seam from loop `a` to loop `b` of a band that winds between the holes where no straight one
/// passes (holes in staggered rows): the shortest way over the free cells of a grid laid over
/// one period of the band.
pub(super) fn seam_route(
    a: &[BPoint],
    b: &[BPoint],
    holes: &[Vec<BPoint>],
    p: f64,
) -> Option<Route> {
    let u0 = a.iter().map(|q| q.uv.0).fold(f64::MAX, f64::min);
    let (v0, v1) = [a, b]
        .into_iter()
        .chain(holes.iter().map(Vec::as_slice))
        .flatten()
        .fold((f64::MAX, f64::MIN), |(lo, hi), q| {
            (lo.min(q.uv.1), hi.max(q.uv.1))
        });
    if v1 <= v0 {
        return None;
    }
    // Segments by owner: 0 for `a`, 1 for `b`, 2 for a hole; copied a period either way so that
    // every part of the window [u0, u0 + p] is covered.
    let mut segments: Vec<(Segment, u8)> = Vec::new();
    for (owner, l) in [(0u8, a), (1, b)]
        .into_iter()
        .chain(holes.iter().map(|h| (2u8, h.as_slice())))
    {
        for (x, y) in loop_segments(l, p) {
            let (lo, hi) = (x.0.min(y.0), x.0.max(y.0));
            let k0 = ((u0 - hi) / p).floor() as i64 - 1;
            let k1 = ((u0 + p - lo) / p).ceil() as i64 + 1;
            for k in k0..=k1 {
                let s = k as f64 * p;
                segments.push((((x.0 + s, x.1), (y.0 + s, y.1)), owner));
            }
        }
    }
    ROUTE_GRIDS.iter().find_map(|&n| {
        let grid = Grid::new(&segments, (u0, v0), (p / n as f64, (v1 - v0) / n as f64), n);
        let cells = grid.shortest(0, 1)?;
        let (from, start) = grid.join(a, cells[0], &segments, p)?;
        let (to, end) = grid.join(b, cells[cells.len() - 1], &segments, p)?;
        let mut corners = vec![start];
        corners.extend(straighten(cells.iter().map(|&c| grid.centre(c))));
        corners.push(end);
        Some(Route { from, to, corners })
    })
}

/// Cells of a grid over one period of a band, each with the segments that touch it; columns
/// wrap round the period.
struct Grid {
    origin: (f64, f64),
    size: (f64, f64),
    n: usize,
    touching: Vec<Vec<u32>>,
    owners: Vec<u8>,
}

impl Grid {
    fn new(segments: &[(Segment, u8)], origin: (f64, f64), size: (f64, f64), n: usize) -> Grid {
        let mut touching = vec![Vec::new(); n * n];
        let mut owners = vec![0u8; n * n];
        for (id, &((x, y), owner)) in segments.iter().enumerate() {
            let col = |u: f64| ((u - origin.0) / size.0).floor() as i64;
            let row = |v: f64| ((v - origin.1) / size.1).floor() as i64;
            let last = n as i64 - 1;
            let (c0, c1) = (
                col(x.0.min(y.0)).clamp(0, last),
                col(x.0.max(y.0)).min(last),
            );
            let (r0, r1) = (
                row(x.1.min(y.1)).clamp(0, last),
                row(x.1.max(y.1)).min(last),
            );
            for c in c0..=c1 {
                for r in r0..=r1 {
                    let lo = (origin.0 + c as f64 * size.0, origin.1 + r as f64 * size.1);
                    if meets_box(x, y, lo, (lo.0 + size.0, lo.1 + size.1)) {
                        let cell = r as usize * n + c as usize;
                        touching[cell].push(id as u32);
                        owners[cell] |= 1 << owner;
                    }
                }
            }
        }
        Grid {
            origin,
            size,
            n,
            touching,
            owners,
        }
    }

    /// Neighbours of a cell given by unwrapped column and row.
    fn neighbours(&self, (c, r): (i64, usize)) -> impl Iterator<Item = (i64, usize)> + '_ {
        let n = self.n;
        [(c - 1, r), (c + 1, r), (c, r.wrapping_sub(1)), (c, r + 1)]
            .into_iter()
            .filter(move |&(_, r)| r < n)
    }

    fn index(&self, (c, r): (i64, usize)) -> usize {
        r * self.n + c.rem_euclid(self.n as i64) as usize
    }

    fn free(&self, cell: (i64, usize)) -> bool {
        self.touching[self.index(cell)].is_empty()
    }

    fn next_to(&self, cell: (i64, usize), owner: u8) -> bool {
        self.neighbours(cell)
            .any(|x| self.owners[self.index(x)] & (1 << owner) != 0)
    }

    /// Free cells from one next to segments of owner `from` to one next to segments of owner
    /// `to`, fewest steps.
    fn shortest(&self, from: u8, to: u8) -> Option<Vec<(i64, usize)>> {
        let n = self.n;
        let mut parent: Vec<Option<(i64, usize)>> = vec![None; n * n];
        let mut seen = vec![false; n * n];
        let mut queue = VecDeque::new();
        for r in 0..n {
            for c in 0..n as i64 {
                if self.free((c, r)) && self.next_to((c, r), from) {
                    seen[self.index((c, r))] = true;
                    queue.push_back((c, r));
                }
            }
        }
        while let Some(cell) = queue.pop_front() {
            if self.next_to(cell, to) {
                let mut path = vec![cell];
                while let Some(prev) = parent[self.index(path[path.len() - 1])] {
                    path.push(prev);
                }
                path.reverse();
                // Columns continue across the wrap from the first cell on.
                let mut col = path[0].0;
                for k in 1..path.len() {
                    let step = (path[k].0 - path[k - 1].0).rem_euclid(n as i64);
                    col += if step == n as i64 - 1 { -1 } else { step };
                    path[k].0 = col;
                }
                return Some(path);
            }
            for next in self.neighbours(cell) {
                let i = self.index(next);
                if !seen[i] && self.free(next) {
                    seen[i] = true;
                    parent[i] = Some(cell);
                    queue.push_back((next.0.rem_euclid(n as i64), next.1));
                }
            }
        }
        None
    }

    fn centre(&self, (c, r): (i64, usize)) -> (f64, f64) {
        (
            self.origin.0 + (c as f64 + 0.5) * self.size.0,
            self.origin.1 + (r as f64 + 0.5) * self.size.1,
        )
    }

    /// The point of `l` the cell's centre sees across no segment, nearest first, moved by whole
    /// periods next to the cell.
    fn join(
        &self,
        l: &[BPoint],
        cell: (i64, usize),
        segments: &[(Segment, u8)],
        p: f64,
    ) -> Option<(usize, (f64, f64))> {
        let c = self.centre(cell);
        let near = |q: &BPoint| {
            let u = q.uv.0 + ((c.0 - q.uv.0) / p).round() * p;
            (u, q.uv.1)
        };
        let dist = |q: (f64, f64)| (q.0 - c.0).powi(2) + (q.1 - c.1).powi(2);
        let mut order: Vec<usize> = (0..l.len()).filter(|&i| !l[i].pole).collect();
        order.sort_by(|&x, &y| dist(near(&l[x])).total_cmp(&dist(near(&l[y]))));
        // The segments cover the window a period either way; the test is made there.
        let back = ((c.0 - self.origin.0) / p).floor() * p;
        order.into_iter().take(ROUTE_ENDS).find_map(|i| {
            let q = near(&l[i]);
            let (qw, cw) = ((q.0 - back, q.1), (c.0 - back, c.1));
            let clear = !segments
                .iter()
                .any(|&((x, y), _)| cross_strictly(qw, cw, x, y));
            clear.then_some((i, q))
        })
    }
}

/// Centres of a path of cells without those in the middle of a straight run.
fn straighten(points: impl Iterator<Item = (f64, f64)>) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::new();
    for q in points {
        if let [.., a, b] = out[..] {
            let turn = (b.0 - a.0) * (q.1 - a.1) - (b.1 - a.1) * (q.0 - a.0);
            if turn.abs() <= 1e-12 * ((b.0 - a.0).abs() + (b.1 - a.1).abs()).powi(2) {
                out.pop();
            }
        }
        out.push(q);
    }
    out
}

/// Whether segment `x`-`y` meets the box from `lo` to `hi`.
fn meets_box(x: (f64, f64), y: (f64, f64), lo: (f64, f64), hi: (f64, f64)) -> bool {
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for (start, delta, min, max) in [(x.0, y.0 - x.0, lo.0, hi.0), (x.1, y.1 - x.1, lo.1, hi.1)] {
        if delta == 0.0 {
            if start < min || start > max {
                return false;
            }
            continue;
        }
        let (a, b) = ((min - start) / delta, (max - start) / delta);
        t0 = t0.max(a.min(b));
        t1 = t1.min(a.max(b));
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// Whether segments `a`-`b` and `c`-`d` cross at a point inside both (touching ends and lying
/// along each other do not count; ends moved by whole periods match only to rounding).
fn cross_strictly(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    let side = |p: (f64, f64), q: (f64, f64), r: (f64, f64)| {
        let (x, y) = ((q.0 - p.0, q.1 - p.1), (r.0 - p.0, r.1 - p.1));
        let turn = x.0 * y.1 - x.1 * y.0;
        let scale = x.0.hypot(x.1) * y.0.hypot(y.1);
        if turn.abs() <= 1e-9 * scale {
            0.0
        } else {
            turn.signum()
        }
    };
    let (s1, s2) = (side(a, b, c), side(a, b, d));
    let (s3, s4) = (side(c, d, a), side(c, d, b));
    s1 * s2 < 0.0 && s3 * s4 < 0.0
}
