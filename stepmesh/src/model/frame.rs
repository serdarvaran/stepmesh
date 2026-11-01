use super::Decoder;
use crate::Error;

pub type P3 = [f64; 3];

pub(crate) fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(crate) fn add(a: P3, b: P3) -> P3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(crate) fn mul(a: P3, s: f64) -> P3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(crate) fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(crate) fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}
pub(crate) fn unit(a: P3) -> Option<P3> {
    let n = norm(a);
    (n > 0.0 && n.is_finite()).then(|| mul(a, 1.0 / n))
}

/// Right-handed orthonormal frame (AXIS2_PLACEMENT_3D), location in millimetres.
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    pub origin: P3,
    pub x: P3,
    pub y: P3,
    pub z: P3,
}

impl Frame {
    /// ISO 10303-42 `build_axes`: `x` is the part of `ref_dir` orthogonal to `z`.
    pub fn new(origin: P3, z: P3, ref_dir: P3) -> Option<Frame> {
        let z = unit(z)?;
        let mut x = sub(ref_dir, mul(z, dot(ref_dir, z)));
        if norm(x) < 1e-12 {
            // ref_dir missing or parallel to z: any perpendicular direction
            let a = if z[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            x = sub(a, mul(z, dot(a, z)));
        }
        let x = unit(x)?;
        Some(Frame {
            origin,
            x,
            y: cross(z, x),
            z,
        })
    }

    /// Local coordinates to the frame's parent space.
    pub fn to_world(&self, p: P3) -> P3 {
        add(
            self.origin,
            add(add(mul(self.x, p[0]), mul(self.y, p[1])), mul(self.z, p[2])),
        )
    }
}

/// Rigid transform; row `i` gives output component `i` as `m[i][0..3]·p + m[i][3]`.
#[derive(Debug, Clone, PartialEq)]
pub struct Transform {
    pub m: [[f64; 4]; 3],
}

impl Transform {
    pub fn identity() -> Transform {
        Transform::from(&Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        })
    }

    /// Placement with location `origin`, axis `z` and reference direction `x` (unnormalized).
    pub fn from_frame(origin: P3, z: P3, x: P3) -> Transform {
        Frame::new(origin, z, x).map_or_else(Transform::identity, |f| Transform::from(&f))
    }

    pub fn apply(&self, p: P3) -> P3 {
        let r = |i: usize| {
            self.m[i][0] * p[0] + self.m[i][1] * p[1] + self.m[i][2] * p[2] + self.m[i][3]
        };
        [r(0), r(1), r(2)]
    }

    /// Applies `inner` first, then `self`.
    pub fn compose(&self, inner: &Transform) -> Transform {
        let mut m = [[0.0; 4]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = (0..3).map(|k| self.m[i][k] * inner.m[k][j]).sum::<f64>()
                    + if j == 3 { self.m[i][3] } else { 0.0 };
            }
        }
        Transform { m }
    }

    /// Inverse of a rigid transform (the rotation part is orthonormal).
    pub fn inverse(&self) -> Transform {
        let mut m = [[0.0; 4]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate().take(3) {
                *v = self.m[j][i];
            }
            row[3] = -(0..3).map(|k| self.m[k][i] * self.m[k][3]).sum::<f64>();
        }
        Transform { m }
    }
}

impl From<&Frame> for Transform {
    fn from(f: &Frame) -> Transform {
        let c = |i: usize| [f.x[i], f.y[i], f.z[i], f.origin[i]];
        Transform {
            m: [c(0), c(1), c(2)],
        }
    }
}

impl Decoder<'_> {
    pub(crate) fn point(&self, id: u64, scale: f64) -> Result<P3, Error> {
        let r = self.record(id, "CARTESIAN_POINT")?;
        let c = self.reals(id, r.params.get(1))?;
        let p = match c[..] {
            [x, y, z] => [x * scale, y * scale, z * scale],
            [x, y] => [x * scale, y * scale, 0.0],
            _ => return self.invalid(id, "point needs 2 or 3 coordinates"),
        };
        if p.iter().all(|v| v.is_finite()) {
            Ok(p)
        } else {
            self.invalid(id, "coordinate out of range")
        }
    }

    pub(crate) fn direction(&self, id: u64) -> Result<P3, Error> {
        let r = self.record(id, "DIRECTION")?;
        let c = self.reals(id, r.params.get(1))?;
        let d = match c[..] {
            [x, y, z] => [x, y, z],
            [x, y] => [x, y, 0.0],
            _ => return self.invalid(id, "direction needs 2 or 3 components"),
        };
        unit(d).map_or_else(|| self.invalid(id, "zero direction"), Ok)
    }

    /// AXIS2_PLACEMENT_3D (location, axis, ref_direction), location scaled to millimetres.
    pub(crate) fn frame(&self, id: u64, scale: f64) -> Result<Frame, Error> {
        let r = self.record(id, "AXIS2_PLACEMENT_3D")?;
        let origin = self.point(self.reference(id, r.params.get(1))?, scale)?;
        let z = match self.optional_reference(id, r.params.get(2))? {
            Some(a) => self.direction(a)?,
            None => [0.0, 0.0, 1.0],
        };
        let x = match self.optional_reference(id, r.params.get(3))? {
            Some(a) => self.direction(a)?,
            None => [1.0, 0.0, 0.0],
        };
        Frame::new(origin, z, x).map_or_else(|| self.invalid(id, "degenerate placement"), Ok)
    }
}
