use super::Decoder;
use crate::Error;
use crate::p21::Value;

/// Conversion chains (inch → millimetre → …) longer than this are treated as cycles.
const MAX_UNIT_CHAIN: usize = 8;

fn si_prefix(p: &str) -> Option<f64> {
    Some(match p {
        "EXA" => 1e18,
        "PETA" => 1e15,
        "TERA" => 1e12,
        "GIGA" => 1e9,
        "MEGA" => 1e6,
        "KILO" => 1e3,
        "HECTO" => 1e2,
        "DECA" => 1e1,
        "DECI" => 1e-1,
        "CENTI" => 1e-2,
        "MILLI" => 1e-3,
        "MICRO" => 1e-6,
        "NANO" => 1e-9,
        "PICO" => 1e-12,
        "FEMTO" => 1e-15,
        "ATTO" => 1e-18,
        _ => return None,
    })
}

/// Length and plane angle units of a representation context.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Units {
    pub length_mm: f64,
    pub angle_rad: f64,
}

impl Decoder<'_> {
    /// Units of a representation's context. A context without a length or angle unit is read
    /// as millimetres and radians, as OpenCascade does.
    pub(crate) fn units(&self, representation: u64) -> Result<Units, Error> {
        let rep = self.first(representation)?;
        let ctx = self.reference(representation, rep.params.get(2))?;
        let mut out = Units {
            length_mm: 1.0,
            angle_rad: 1.0,
        };
        let Some(assigned) = self
            .entity(ctx)?
            .records
            .iter()
            .find(|r| r.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT")
        else {
            return Ok(out);
        };
        for u in self.references(ctx, assigned.params.first())? {
            let e = self.entity(u)?;
            if e.records.iter().any(|r| r.name == "LENGTH_UNIT") {
                out.length_mm = self.unit_factor(u, "METRE", 1000.0, 0)?;
            } else if e.records.iter().any(|r| r.name == "PLANE_ANGLE_UNIT") {
                out.angle_rad = self.unit_factor(u, "RADIAN", 1.0, 0)?;
            }
        }
        Ok(out)
    }

    /// Size of unit `id` in the target unit, where one SI `base` unit is `si` target units.
    fn unit_factor(&self, id: u64, base: &str, si: f64, depth: usize) -> Result<f64, Error> {
        if depth > MAX_UNIT_CHAIN {
            return self.invalid(id, "unit conversion chain too long or cyclic");
        }
        let e = self.entity(id)?;
        if let Some(r) = e.records.iter().find(|r| r.name == "SI_UNIT") {
            if !matches!(r.params.get(1), Some(Value::Enum(n)) if n == base) {
                return self.invalid(
                    id,
                    &format!("unit is not based on the {}", base.to_lowercase()),
                );
            }
            let prefix = match r.params.first() {
                Some(Value::Enum(p)) => {
                    si_prefix(p).map_or_else(|| self.invalid(id, "unknown SI prefix"), Ok)?
                }
                _ => 1.0,
            };
            return Ok(prefix * si);
        }
        if let Some(c) = e.records.iter().find(|r| r.name == "CONVERSION_BASED_UNIT") {
            let m = self.reference(id, c.params.get(1))?;
            let mr = self.first(m)?;
            let value = self.real(m, mr.params.first())?;
            let unit = self.reference(m, mr.params.get(1))?;
            return Ok(value * self.unit_factor(unit, base, si, depth + 1)?);
        }
        Err(Error::Unsupported {
            entity: e.records[0].name.clone(),
            id,
        })
    }
}
