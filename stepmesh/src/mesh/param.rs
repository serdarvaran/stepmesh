//! A surface seen with its parameters in either order: the face code closes loops that go round
//! along u, and meshes faces whose loops go round along v with u and v swapped.

use crate::Options;
use crate::geom::Surface;
use crate::geom::bounds::NormalBounds;
use crate::model::{P3, mul};

#[derive(Clone, Copy)]
pub(super) struct Param<'a> {
    pub surface: &'a Surface,
    /// u and v exchanged: the normal of the view (∂u × ∂v) is the surface's reversed.
    pub swapped: bool,
}

impl Param<'_> {
    pub fn triangle_within(
        &self,
        uv: [(f64, f64); 3],
        pts: [P3; 3],
        o: &Options,
        cache: &NormalBounds,
    ) -> Option<bool> {
        crate::geom::bounds::triangle_within(self.surface, uv.map(|p| self.order(p)), pts, o, cache)
    }
    pub fn fan_within(&self, uv: [(f64, f64); 3], pts: [P3; 3], o: &Options) -> bool {
        crate::geom::bounds::fan_within(
            self.surface,
            uv.map(|p| self.order(p)),
            pts,
            o.linear_deflection,
        )
    }
    fn order<T>(&self, (a, b): (T, T)) -> (T, T) {
        if self.swapped { (b, a) } else { (a, b) }
    }

    pub fn point(&self, u: f64, v: f64) -> P3 {
        let (u, v) = self.order((u, v));
        self.surface.point(u, v)
    }

    pub fn derivs(&self, u: f64, v: f64) -> (P3, P3, P3) {
        let (u, v) = self.order((u, v));
        let (p, du, dv) = self.surface.derivs(u, v);
        let (du, dv) = self.order((du, dv));
        (p, du, dv)
    }

    pub fn normal(&self, u: f64, v: f64) -> P3 {
        let (u, v) = self.order((u, v));
        let n = self.surface.normal(u, v);
        if self.swapped { mul(n, -1.0) } else { n }
    }

    pub fn project(&self, p: P3, hint: Option<(f64, f64)>) -> Option<(f64, f64)> {
        self.surface
            .project(p, hint.map(|h| self.order(h)))
            .map(|q| self.order(q))
    }

    pub fn foot(&self, p: P3, start: (f64, f64)) -> Option<(f64, f64)> {
        self.surface
            .foot(p, self.order(start))
            .map(|q| self.order(q))
    }

    pub fn period(&self) -> (Option<f64>, Option<f64>) {
        self.order(self.surface.period())
    }

    pub fn param_bounds(&self) -> [(f64, f64); 2] {
        let [a, b] = self.surface.param_bounds();
        let (a, b) = self.order((a, b));
        [a, b]
    }

    pub fn breaks(&self) -> [Vec<f64>; 2] {
        let [a, b] = self.surface.breaks();
        let (a, b) = self.order((a, b));
        [a, b]
    }
}
