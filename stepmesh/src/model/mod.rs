//! Typed reading of the exchange structure: units, placements, product structure.

mod assembly;
mod frame;
mod geometry;
mod units;

pub use assembly::decode_model;
pub use frame::{Frame, P3, Transform};
pub(crate) use frame::{add, cross, dot, mul, norm, sub, unit};
pub use units::Units;

use crate::Error;
use crate::p21::{Entity, Exchange, Record, Value};

/// One placed use of a part.
#[derive(Debug, Clone, PartialEq)]
pub struct Instance {
    pub part: usize,
    pub transform: Transform,
}

/// A product definition (or a mapped representation): its solids and its placed sub-parts.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Part {
    pub name: Option<String>,
    /// B-rep solid instances with the units of the representation holding them.
    pub breps: Vec<(u64, Units)>,
    pub children: Vec<Instance>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Model {
    pub parts: Vec<Part>,
    pub roots: Vec<Instance>,
}

/// Typed access to the instances of one exchange structure.
pub(crate) struct Decoder<'a> {
    pub x: &'a Exchange,
}

impl<'a> Decoder<'a> {
    pub fn new(x: &'a Exchange) -> Self {
        Decoder { x }
    }

    pub fn invalid<T>(&self, id: u64, msg: &str) -> Result<T, Error> {
        Err(Error::Invalid {
            id,
            msg: msg.into(),
        })
    }

    pub fn entity(&self, id: u64) -> Result<&'a Entity, Error> {
        self.x.entities.get(&id).ok_or(Error::Invalid {
            id,
            msg: "undefined instance".into(),
        })
    }

    /// The record of a simple instance, or the first one of a complex instance.
    pub fn first(&self, id: u64) -> Result<&'a Record, Error> {
        Ok(&self.entity(id)?.records[0])
    }

    /// The record named `name` (simple or complex instance), or `Invalid`.
    pub fn record(&self, id: u64, name: &str) -> Result<&'a Record, Error> {
        match self.entity(id)?.records.iter().find(|r| r.name == name) {
            Some(r) => Ok(r),
            None => self.invalid(id, &format!("expected {name}")),
        }
    }

    pub fn reference(&self, id: u64, v: Option<&Value>) -> Result<u64, Error> {
        match v {
            Some(Value::Ref(r)) => Ok(*r),
            _ => self.invalid(id, "expected an instance reference"),
        }
    }

    pub fn optional_reference(&self, id: u64, v: Option<&Value>) -> Result<Option<u64>, Error> {
        match v {
            Some(Value::Ref(r)) => Ok(Some(*r)),
            Some(Value::Unset) | None => Ok(None),
            _ => self.invalid(id, "expected an instance reference or $"),
        }
    }

    pub fn references(&self, id: u64, v: Option<&Value>) -> Result<Vec<u64>, Error> {
        match v {
            Some(Value::List(items)) => items.iter().map(|i| self.reference(id, Some(i))).collect(),
            _ => self.invalid(id, "expected a list of instance references"),
        }
    }

    pub fn real(&self, id: u64, v: Option<&Value>) -> Result<f64, Error> {
        match v {
            Some(Value::Real(r)) if r.is_finite() => Ok(*r),
            Some(Value::Int(i)) => Ok(*i as f64),
            Some(Value::Typed(_, inner)) => self.real(id, Some(inner)),
            _ => self.invalid(id, "expected a number"),
        }
    }

    pub fn reals(&self, id: u64, v: Option<&Value>) -> Result<Vec<f64>, Error> {
        match v {
            Some(Value::List(items)) => items.iter().map(|i| self.real(id, Some(i))).collect(),
            _ => self.invalid(id, "expected a list of numbers"),
        }
    }

    pub fn string(&self, v: Option<&Value>) -> Option<String> {
        match v {
            Some(Value::Str(s)) if !s.is_empty() => Some(s.clone()),
            _ => None,
        }
    }

    /// Ids of all instances with a record named `name`, sorted for a deterministic order.
    pub fn all(&self, name: &str) -> Vec<u64> {
        let mut ids: Vec<u64> = self
            .x
            .entities
            .iter()
            .filter(|(_, e)| e.records.iter().any(|r| r.name == name))
            .map(|(id, _)| *id)
            .collect();
        ids.sort_unstable();
        ids
    }
}

/// Curve instance `id` (tests and tools).
pub fn decode_curve(x: &Exchange, id: u64, units: Units) -> Result<crate::geom::Curve, Error> {
    Decoder::new(x).curve(id, units)
}

/// Surface instance `id` (tests and tools).
pub fn decode_surface(x: &Exchange, id: u64, units: Units) -> Result<crate::geom::Surface, Error> {
    Decoder::new(x).surface(id, units)
}
