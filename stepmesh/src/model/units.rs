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

impl Decoder<'_> {
    /// Millimetres per length unit of a representation's context. A context without a length unit is
    /// read as millimetres, as OpenCascade does.
    pub(crate) fn length_scale_mm(&self, representation: u64) -> Result<f64, Error> {
        let rep = self.first(representation)?;
        let ctx = self.reference(representation, rep.params.get(2))?;
        let Some(units) = self
            .entity(ctx)?
            .records
            .iter()
            .find(|r| r.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT")
        else {
            return Ok(1.0);
        };
        let Some(Value::List(units)) = units.params.first() else {
            return self.invalid(ctx, "GLOBAL_UNIT_ASSIGNED_CONTEXT without units");
        };
        for u in units {
            let Value::Ref(u) = *u else { continue };
            if self
                .entity(u)?
                .records
                .iter()
                .any(|r| r.name == "LENGTH_UNIT")
            {
                return self.unit_mm(u, 0);
            }
        }
        Ok(1.0)
    }

    fn unit_mm(&self, id: u64, depth: usize) -> Result<f64, Error> {
        if depth > MAX_UNIT_CHAIN {
            return self.invalid(id, "unit conversion chain too long or cyclic");
        }
        let e = self.entity(id)?;
        if let Some(si) = e.records.iter().find(|r| r.name == "SI_UNIT") {
            if !matches!(si.params.get(1), Some(Value::Enum(n)) if n == "METRE") {
                return self.invalid(id, "length unit is not based on the metre");
            }
            let factor = match si.params.first() {
                Some(Value::Enum(p)) => {
                    si_prefix(p).map_or_else(|| self.invalid(id, "unknown SI prefix"), Ok)?
                }
                _ => 1.0,
            };
            return Ok(factor * 1000.0);
        }
        if let Some(c) = e.records.iter().find(|r| r.name == "CONVERSION_BASED_UNIT") {
            let m = self.reference(id, c.params.get(1))?;
            let mr = self.first(m)?;
            let value = match mr.params.first() {
                Some(Value::Typed(_, v)) => self.real(m, Some(v))?,
                v => self.real(m, v)?,
            };
            let base = self.reference(m, mr.params.get(1))?;
            return Ok(value * self.unit_mm(base, depth + 1)?);
        }
        Err(Error::Unsupported {
            entity: e.records[0].name.clone(),
            id,
        })
    }
}
