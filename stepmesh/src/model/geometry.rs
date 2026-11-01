//! Curves and surfaces of the exchange structure, in millimetres and radians.

use super::{Decoder, P3, Units, unit};
use crate::Error;
use crate::geom::{Curve, Nurbs1, Nurbs2, Surface, expand_knots};
use crate::p21::{Record, Value};

/// Records that give a B-spline curve its knots, explicitly or by form.
const CURVE_FORMS: [&str; 4] = [
    "B_SPLINE_CURVE_WITH_KNOTS",
    "UNIFORM_CURVE",
    "QUASI_UNIFORM_CURVE",
    "BEZIER_CURVE",
];
const SURFACE_FORMS: [&str; 4] = [
    "B_SPLINE_SURFACE_WITH_KNOTS",
    "UNIFORM_SURFACE",
    "QUASI_UNIFORM_SURFACE",
    "BEZIER_SURFACE",
];

impl Decoder<'_> {
    fn unsupported<T>(&self, id: u64) -> Result<T, Error> {
        Err(Error::Unsupported {
            entity: self.first(id)?.name.clone(),
            id,
        })
    }

    fn integer(&self, id: u64, v: Option<&Value>) -> Result<usize, Error> {
        match v {
            Some(Value::Int(i)) if *i >= 0 => Ok(*i as usize),
            _ => self.invalid(id, "expected a non-negative integer"),
        }
    }

    fn integers(&self, id: u64, v: Option<&Value>) -> Result<Vec<usize>, Error> {
        match v {
            Some(Value::List(items)) => items.iter().map(|i| self.integer(id, Some(i))).collect(),
            _ => self.invalid(id, "expected a list of integers"),
        }
    }

    fn logical(&self, v: Option<&Value>) -> bool {
        matches!(v, Some(Value::Enum(e)) if e == "T")
    }

    /// VECTOR (orientation, magnitude) as a unit direction; the magnitude only scales the parameter.
    fn vector_direction(&self, id: u64) -> Result<P3, Error> {
        let r = self.record(id, "VECTOR")?;
        self.direction(self.reference(id, r.params.get(1))?)
    }

    pub(crate) fn curve(&self, id: u64, units: Units) -> Result<Curve, Error> {
        let e = self.entity(id)?;
        if e.records
            .iter()
            .any(|r| r.name == "B_SPLINE_CURVE" || CURVE_FORMS.contains(&r.name.as_str()))
        {
            return self.bspline_curve(id, units);
        }
        let r = &e.records[0];
        let s = units.length_mm;
        match r.name.as_str() {
            "LINE" => {
                let origin = self.point(self.reference(id, r.params.get(1))?, s)?;
                let dir = self.vector_direction(self.reference(id, r.params.get(2))?)?;
                Ok(Curve::Line { origin, dir })
            }
            "CIRCLE" => Ok(Curve::Circle {
                frame: self.frame(self.reference(id, r.params.get(1))?, s)?,
                r: self.real(id, r.params.get(2))? * s,
            }),
            "ELLIPSE" => Ok(Curve::Ellipse {
                frame: self.frame(self.reference(id, r.params.get(1))?, s)?,
                a: self.real(id, r.params.get(2))? * s,
                b: self.real(id, r.params.get(3))? * s,
            }),
            // The 3D curve; pcurves are not read (boundary points are projected onto the surface).
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" => {
                self.curve(self.reference(id, r.params.get(1))?, units)
            }
            // Edges are trimmed by their vertices, so the basis curve is enough.
            "TRIMMED_CURVE" => self.curve(self.reference(id, r.params.get(1))?, units),
            _ => self.unsupported(id),
        }
    }

    /// The B-spline's base record, the record that gives its knots (one of `forms`), and the
    /// index of the degree in the base record. A simple instance carries everything on one
    /// record after the name; a complex one splits it by supertype.
    fn bspline_record(
        &self,
        id: u64,
        base: &str,
        forms: &[&str],
    ) -> Result<(&Record, &Record, usize), Error> {
        let e = self.entity(id)?;
        if e.records.len() == 1 {
            return Ok((&e.records[0], &e.records[0], 1));
        }
        let b = e
            .records
            .iter()
            .find(|r| r.name == base)
            .map_or_else(|| self.invalid(id, base), Ok)?;
        let knots = e
            .records
            .iter()
            .find(|r| forms.contains(&r.name.as_str()))
            .map_or_else(|| self.invalid(id, "B-spline without a knot form"), Ok)?;
        Ok((b, knots, 0))
    }

    /// Knot vector of `n` control points of degree `p`: explicit (`*_WITH_KNOTS`, from the
    /// multiplicities and values at `at` and `at + offset`) or implied by the form (ISO 10303-42).
    fn knot_vector(
        &self,
        id: u64,
        form: &Record,
        at: usize,
        offset: usize,
        n: usize,
        p: usize,
    ) -> Result<Vec<f64>, Error> {
        let kind = form.name.as_str();
        if kind.ends_with("_WITH_KNOTS") {
            let mults = self.integers(id, form.params.get(at))?;
            let ks = self.reals(id, form.params.get(at + offset))?;
            return Ok(expand_knots(&mults, &ks));
        }
        if n <= p {
            return self.invalid(id, "fewer control points than the degree needs");
        }
        Ok(if kind.starts_with("UNIFORM") {
            (0..n + p + 1).map(|i| i as f64 - p as f64).collect()
        } else if kind.starts_with("QUASI_UNIFORM") {
            let last = (n - p) as f64;
            let mut k = vec![0.0; p + 1];
            k.extend((1..n - p).map(|i| i as f64));
            k.extend(std::iter::repeat_n(last, p + 1));
            k
        } else {
            // BEZIER: piecewise Bézier segments of degree p
            if (n - 1) % p != 0 {
                return self.invalid(id, "control points do not form whole Bézier segments");
            }
            let segments = (n - 1) / p;
            let mut mults = vec![p; segments + 1];
            mults[0] = p + 1;
            mults[segments] = p + 1;
            expand_knots(
                &mults,
                &(0..=segments).map(|i| i as f64).collect::<Vec<_>>(),
            )
        })
    }

    fn weights(&self, id: u64, record: &str) -> Result<Option<&Value>, Error> {
        Ok(self
            .entity(id)?
            .records
            .iter()
            .find(|r| r.name == record)
            .and_then(|r| r.params.first()))
    }

    fn bspline_curve(&self, id: u64, units: Units) -> Result<Curve, Error> {
        let (base, form, off) = self.bspline_record(id, "B_SPLINE_CURVE", &CURVE_FORMS)?;
        let degree = self.integer(id, base.params.get(off))?;
        let pts = self.references(id, base.params.get(off + 1))?;
        let k0 = if off == 1 { 6 } else { 0 };
        let w = match self.weights(id, "RATIONAL_B_SPLINE_CURVE")? {
            Some(v) => self.reals(id, Some(v))?,
            None => vec![1.0; pts.len()],
        };
        let knots = self.knot_vector(id, form, k0, 1, pts.len(), degree)?;
        if degree == 0 || w.len() != pts.len() || knots.len() != pts.len() + degree + 1 {
            return self.invalid(id, "inconsistent B-spline curve");
        }
        let ctrl = pts
            .iter()
            .zip(&w)
            .map(|(&p, &w)| {
                self.point(p, units.length_mm)
                    .map(|q| [q[0] * w, q[1] * w, q[2] * w, w])
            })
            .collect::<Result<_, _>>()?;
        Ok(Curve::Nurbs(Nurbs1 {
            degree,
            knots,
            ctrl,
        }))
    }

    pub(crate) fn surface(&self, id: u64, units: Units) -> Result<Surface, Error> {
        let e = self.entity(id)?;
        if e.records
            .iter()
            .any(|r| r.name == "B_SPLINE_SURFACE" || SURFACE_FORMS.contains(&r.name.as_str()))
        {
            return self.bspline_surface(id, units);
        }
        let r = &e.records[0];
        let s = units.length_mm;
        let frame = |i: usize| self.frame(self.reference(id, r.params.get(i))?, s);
        let len = |i: usize| self.real(id, r.params.get(i)).map(|v| v * s);
        match r.name.as_str() {
            "PLANE" => Ok(Surface::Plane(frame(1)?)),
            "CYLINDRICAL_SURFACE" => Ok(Surface::Cylinder {
                frame: frame(1)?,
                r: len(2)?,
            }),
            "CONICAL_SURFACE" => Ok(Surface::Cone {
                frame: frame(1)?,
                r: len(2)?,
                half_angle: self.real(id, r.params.get(3))? * units.angle_rad,
            }),
            "SPHERICAL_SURFACE" => Ok(Surface::Sphere {
                frame: frame(1)?,
                r: len(2)?,
            }),
            "TOROIDAL_SURFACE" => Ok(Surface::Torus {
                frame: frame(1)?,
                major: len(2)?,
                minor: len(3)?,
            }),
            "SURFACE_OF_LINEAR_EXTRUSION" => Ok(Surface::Extrusion {
                curve: self.curve(self.reference(id, r.params.get(1))?, units)?,
                dir: self.vector_direction(self.reference(id, r.params.get(2))?)?,
            }),
            "SURFACE_OF_REVOLUTION" => {
                let axis = self.reference(id, r.params.get(2))?;
                let a = self.record(axis, "AXIS1_PLACEMENT")?;
                let axis_origin = self.point(self.reference(axis, a.params.get(1))?, s)?;
                let axis_dir = match self.optional_reference(axis, a.params.get(2))? {
                    Some(d) => self.direction(d)?,
                    None => [0.0, 0.0, 1.0],
                };
                Ok(Surface::Revolution {
                    curve: self.curve(self.reference(id, r.params.get(1))?, units)?,
                    axis_origin,
                    axis_dir: unit(axis_dir).unwrap_or([0.0, 0.0, 1.0]),
                })
            }
            _ => self.unsupported(id),
        }
    }

    fn bspline_surface(&self, id: u64, units: Units) -> Result<Surface, Error> {
        let (base, form, off) = self.bspline_record(id, "B_SPLINE_SURFACE", &SURFACE_FORMS)?;
        let du = self.integer(id, base.params.get(off))?;
        let dv = self.integer(id, base.params.get(off + 1))?;
        let Some(Value::List(rows)) = base.params.get(off + 2) else {
            return self.invalid(id, "B-spline surface without control points");
        };
        let closed_u = self.logical(base.params.get(off + 4));
        let closed_v = self.logical(base.params.get(off + 5));
        let k0 = if off == 1 { 8 } else { 0 };
        let nu = rows.len();
        let nv = match rows.first() {
            Some(Value::List(r)) => r.len(),
            _ => return self.invalid(id, "B-spline surface without control points"),
        };
        let ku = self.knot_vector(id, form, k0, 2, nu, du)?;
        let kv = self.knot_vector(id, form, k0 + 1, 2, nv, dv)?;
        let weights: Option<Vec<Vec<f64>>> = match self.weights(id, "RATIONAL_B_SPLINE_SURFACE")? {
            Some(Value::List(ws)) => Some(
                ws.iter()
                    .map(|w| self.reals(id, Some(w)))
                    .collect::<Result<_, _>>()?,
            ),
            Some(_) => return self.invalid(id, "bad weights"),
            None => None,
        };
        let mut ctrl = Vec::with_capacity(rows.len());
        for (i, row) in rows.iter().enumerate() {
            let pts = self.references(id, Some(row))?;
            let mut out = Vec::with_capacity(pts.len());
            for (j, &p) in pts.iter().enumerate() {
                let w = match &weights {
                    Some(ws) => *ws
                        .get(i)
                        .and_then(|r| r.get(j))
                        .map_or_else(|| self.invalid(id, "weights"), Ok)?,
                    None => 1.0,
                };
                let q = self.point(p, units.length_mm)?;
                out.push([q[0] * w, q[1] * w, q[2] * w, w]);
            }
            ctrl.push(out);
        }
        let (nu, nv) = (ctrl.len(), ctrl.first().map_or(0, Vec::len));
        if du == 0
            || dv == 0
            || nu == 0
            || ctrl.iter().any(|r| r.len() != nv)
            || ku.len() != nu + du + 1
            || kv.len() != nv + dv + 1
        {
            return self.invalid(id, "inconsistent B-spline surface");
        }
        Ok(Surface::Nurbs(Nurbs2 {
            du,
            dv,
            ku,
            kv,
            ctrl,
            closed_u,
            closed_v,
        }))
    }
}
