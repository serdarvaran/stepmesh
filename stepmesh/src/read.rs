//! `read_step`: from file bytes to placed solid meshes.

use std::collections::HashMap;
use std::rc::Rc;

use crate::mesh::{LocalMesh, mesh_brep, signed_volume, triangle_area};
use crate::model::{Decoder, Model, Transform, Units, decode_model};
use crate::{Error, Options, p21};

/// One placed solid, in millimetres. Triangles are counter-clockwise seen from outside.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidMesh {
    /// Name of the product the solid belongs to, when the file gives one.
    pub name: Option<String>,
    pub positions: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
}

impl SolidMesh {
    /// Enclosed volume in mm³ (voids subtracted).
    pub fn volume(&self) -> f64 {
        signed_volume(&self.positions, &self.triangles)
    }

    /// Surface area in mm².
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.positions[i as usize]);
                triangle_area(a, b, c)
            })
            .sum()
    }

    /// Smallest and largest coordinates.
    pub fn bounds(&self) -> [[f64; 3]; 2] {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &self.positions {
            for i in 0..3 {
                lo[i] = lo[i].min(p[i]);
                hi[i] = hi[i].max(p[i]);
            }
        }
        [lo, hi]
    }
}

/// Every solid of a STEP file, one per placed instance.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StepMesh {
    pub solids: Vec<SolidMesh>,
}

/// Reads a STEP file and meshes each of its B-rep solids (assembly placements and units applied).
///
/// A solid whose mesh would not be closed and consistently oriented, or that uses an entity this
/// crate does not read, fails the whole call: no partial or approximate result is returned.
pub fn read_step(bytes: &[u8], options: &Options) -> Result<StepMesh, Error> {
    let positive = |v: f64| v.is_finite() && v > 0.0;
    if !positive(options.linear_deflection) {
        return Err(Error::InvalidOptions(
            "linear_deflection must be finite and positive",
        ));
    }
    if !positive(options.angular_deflection) || options.angular_deflection > std::f64::consts::PI {
        return Err(Error::InvalidOptions(
            "angular_deflection must be in (0, pi]",
        ));
    }
    let x = p21::parse(bytes)?;
    let model = decode_model(&x)?;
    let d = Decoder::new(&x);
    let mut placed = Vec::new();
    for root in &model.roots {
        collect(&model, root.part, &root.transform, &mut placed);
    }
    let mut cache: HashMap<u64, Rc<LocalMesh>> = HashMap::new();
    let mut out = StepMesh::default();
    let mut total = 0usize;
    for (solid, (brep, units, transform, name)) in placed.into_iter().enumerate() {
        let local = match cache.get(&brep) {
            Some(m) => m.clone(),
            None => {
                let b = d.brep(brep, units)?;
                let budget = options.max_triangles.map(|l| l.saturating_sub(total));
                let m = Rc::new(mesh_brep(&b, options, solid, budget)?);
                cache.insert(brep, m.clone());
                m
            }
        };
        total += local.triangles.len();
        if let Some(limit) = options.max_triangles {
            if total > limit {
                return Err(Error::TriangleLimit {
                    count: total,
                    limit,
                });
            }
        }
        out.solids.push(SolidMesh {
            name,
            positions: local
                .positions
                .iter()
                .map(|&p| transform.apply(p))
                .collect(),
            triangles: local.triangles.clone(),
        });
    }
    Ok(out)
}

type Placed = (u64, Units, Transform, Option<String>);

fn collect(m: &Model, part: usize, t: &Transform, out: &mut Vec<Placed>) {
    let p = &m.parts[part];
    out.extend(
        p.breps
            .iter()
            .map(|&(id, units)| (id, units, t.clone(), p.name.clone())),
    );
    for c in &p.children {
        collect(m, c.part, &t.compose(&c.transform), out);
    }
}
