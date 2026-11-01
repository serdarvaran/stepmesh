//! Curves and surfaces of the exchange structure, in millimetres and radians.

use super::{Decoder, P3, Units, unit};
use crate::Error;
use crate::geom::{Curve, Nurbs1, Nurbs2, Surface, expand_knots};
use crate::p21::{Record, Value};

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
            .any(|r| r.name == "B_SPLINE_CURVE_WITH_KNOTS")
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

    fn bspline_record(
        &self,
        id: u64,
        simple: &str,
        complex: &str,
    ) -> Result<(&Record, &Record, usize), Error> {
        // Simple instance: all attributes on one record after the name. Complex: split by supertype.
        let e = self.entity(id)?;
        if e.records.len() == 1 {
            return Ok((&e.records[0], &e.records[0], 1));
        }
        let base = e
            .records
            .iter()
            .find(|r| r.name == simple)
            .map_or_else(|| self.invalid(id, simple), Ok)?;
        let knots = e
            .records
            .iter()
            .find(|r| r.name == complex)
            .map_or_else(|| self.invalid(id, complex), Ok)?;
        Ok((base, knots, 0))
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
        let (base, knots, off) =
            self.bspline_record(id, "B_SPLINE_CURVE", "B_SPLINE_CURVE_WITH_KNOTS")?;
        let degree = self.integer(id, base.params.get(off))?;
        let pts = self.references(id, base.params.get(off + 1))?;
        let k0 = if off == 1 { 6 } else { 0 };
        let mults = self.integers(id, knots.params.get(k0))?;
        let ks = self.reals(id, knots.params.get(k0 + 1))?;
        let w = match self.weights(id, "RATIONAL_B_SPLINE_CURVE")? {
            Some(v) => self.reals(id, Some(v))?,
            None => vec![1.0; pts.len()],
        };
        let knots = expand_knots(&mults, &ks);
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
            .any(|r| r.name == "B_SPLINE_SURFACE_WITH_KNOTS")
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
        let (base, knots, off) =
            self.bspline_record(id, "B_SPLINE_SURFACE", "B_SPLINE_SURFACE_WITH_KNOTS")?;
        let du = self.integer(id, base.params.get(off))?;
        let dv = self.integer(id, base.params.get(off + 1))?;
        let Some(Value::List(rows)) = base.params.get(off + 2) else {
            return self.invalid(id, "B-spline surface without control points");
        };
        let closed_u = self.logical(base.params.get(off + 4));
        let closed_v = self.logical(base.params.get(off + 5));
        let k0 = if off == 1 { 8 } else { 0 };
        let ku = expand_knots(
            &self.integers(id, knots.params.get(k0))?,
            &self.reals(id, knots.params.get(k0 + 2))?,
        );
        let kv = expand_knots(
            &self.integers(id, knots.params.get(k0 + 1))?,
            &self.reals(id, knots.params.get(k0 + 3))?,
        );
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
