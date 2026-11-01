//! B-rep topology: shells, faces, loops and edges of one solid.

use std::collections::HashMap;

use super::{Decoder, P3, Units};
use crate::Error;
use crate::geom::{Curve, Surface};
use crate::p21::Value;

/// EDGE_CURVE in millimetres. `same_sense` tells whether the curve runs from `start` to `end`.
#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    pub curve: Curve,
    pub start: P3,
    pub end: P3,
    pub start_vertex: u64,
    pub end_vertex: u64,
    pub same_sense: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Loop {
    /// EDGE_CURVE ids with their ORIENTED_EDGE orientation.
    Edges(Vec<(u64, bool)>),
    /// VERTEX_LOOP: a degenerate boundary such as the pole of a sphere.
    Vertex(P3),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Face {
    pub id: u64,
    pub surface: Surface,
    pub same_sense: bool,
    /// Index of the FACE_OUTER_BOUND in `loops`, if the file marks one.
    pub outer: Option<usize>,
    /// Bounds with their FACE_BOUND orientation.
    pub loops: Vec<(Loop, bool)>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Shell {
    pub id: u64,
    pub faces: Vec<Face>,
    /// An inner boundary of a BREP_WITH_VOIDS.
    pub void: bool,
    /// ORIENTED_CLOSED_SHELL with orientation `.F.`: the faces are used reversed.
    pub reversed: bool,
}

/// One solid in millimetres, in the coordinates of its representation.
#[derive(Debug, Clone, PartialEq)]
pub struct Brep {
    pub id: u64,
    pub shells: Vec<Shell>,
    pub edges: HashMap<u64, Edge>,
}

impl Decoder<'_> {
    fn flag(&self, id: u64, v: Option<&Value>) -> Result<bool, Error> {
        match v {
            Some(Value::Enum(e)) if e == "T" => Ok(true),
            Some(Value::Enum(e)) if e == "F" => Ok(false),
            _ => self.invalid(id, "expected .T. or .F."),
        }
    }

    pub(crate) fn brep(&self, id: u64, units: Units) -> Result<Brep, Error> {
        let e = self.entity(id)?;
        let r = e
            .records
            .iter()
            .find(|r| r.name == "MANIFOLD_SOLID_BREP" || r.name == "BREP_WITH_VOIDS")
            .map_or_else(|| self.invalid(id, "expected a B-rep solid"), Ok)?;
        let mut b = Brep {
            id,
            shells: Vec::new(),
            edges: HashMap::new(),
        };
        let outer = self.reference(id, r.params.get(1))?;
        let shell = self.shell(outer, false, false, units, &mut b.edges)?;
        b.shells.push(shell);
        // BREP_WITH_VOIDS adds the voids after the outer shell (simple or complex instance).
        let voids = match e.records.iter().find(|r| r.name == "BREP_WITH_VOIDS") {
            Some(v) if v.params.len() >= 3 => Some(&v.params[2]),
            Some(v) => v.params.first(),
            None => None,
        };
        if let Some(voids) = voids {
            for v in self.references(id, Some(voids))? {
                let shell = self.shell(v, true, false, units, &mut b.edges)?;
                b.shells.push(shell);
            }
        }
        Ok(b)
    }

    fn shell(
        &self,
        id: u64,
        void: bool,
        reversed: bool,
        units: Units,
        edges: &mut HashMap<u64, Edge>,
    ) -> Result<Shell, Error> {
        let r = self.first(id)?;
        match r.name.as_str() {
            "CLOSED_SHELL" => {
                let faces = self
                    .references(id, r.params.get(1))?
                    .into_iter()
                    .map(|f| self.face(f, units, edges))
                    .collect::<Result<_, _>>()?;
                Ok(Shell {
                    id,
                    faces,
                    void,
                    reversed,
                })
            }
            "ORIENTED_CLOSED_SHELL" => {
                let inner = self.reference(id, r.params.get(2))?;
                let orientation = self.flag(id, r.params.get(3))?;
                if self.first(inner)?.name != "CLOSED_SHELL" {
                    return self.invalid(id, "an oriented shell must orient a CLOSED_SHELL");
                }
                self.shell(inner, void, reversed != !orientation, units, edges)
                    .map(|s| Shell { id, ..s })
            }
            "OPEN_SHELL" | "ORIENTED_OPEN_SHELL" => Err(Error::OpenShell { id }),
            other => Err(Error::Unsupported {
                entity: other.to_string(),
                id,
            }),
        }
    }

    fn face(&self, id: u64, units: Units, edges: &mut HashMap<u64, Edge>) -> Result<Face, Error> {
        let r = self.first(id)?;
        if r.name != "ADVANCED_FACE" && r.name != "FACE_SURFACE" {
            return Err(Error::Unsupported {
                entity: r.name.clone(),
                id,
            });
        }
        let surface = self.surface(self.reference(id, r.params.get(2))?, units)?;
        let same_sense = self.flag(id, r.params.get(3))?;
        let mut loops = Vec::new();
        let mut outer = None;
        for bound in self.references(id, r.params.get(1))? {
            let b = self.first(bound)?;
            match b.name.as_str() {
                "FACE_OUTER_BOUND" => outer = Some(loops.len()),
                "FACE_BOUND" => {}
                other => {
                    return Err(Error::Unsupported {
                        entity: other.to_string(),
                        id: bound,
                    });
                }
            }
            let l = self.edge_loop(self.reference(bound, b.params.get(1))?, units, edges)?;
            loops.push((l, self.flag(bound, b.params.get(2))?));
        }
        Ok(Face {
            id,
            surface,
            same_sense,
            outer,
            loops,
        })
    }

    fn edge_loop(
        &self,
        id: u64,
        units: Units,
        edges: &mut HashMap<u64, Edge>,
    ) -> Result<Loop, Error> {
        let r = self.first(id)?;
        match r.name.as_str() {
            "EDGE_LOOP" => {
                let mut out = Vec::new();
                for oe in self.references(id, r.params.get(1))? {
                    let o = self.record(oe, "ORIENTED_EDGE")?;
                    let edge = self.reference(oe, o.params.get(3))?;
                    let orientation = self.flag(oe, o.params.get(4))?;
                    if let std::collections::hash_map::Entry::Vacant(slot) = edges.entry(edge) {
                        slot.insert(self.edge(edge, units)?);
                    }
                    out.push((edge, orientation));
                }
                if out.is_empty() {
                    return self.invalid(id, "empty edge loop");
                }
                Ok(Loop::Edges(out))
            }
            "VERTEX_LOOP" => {
                let v = self.reference(id, r.params.get(1))?;
                Ok(Loop::Vertex(self.vertex(v, units)?))
            }
            other => Err(Error::Unsupported {
                entity: other.to_string(),
                id,
            }),
        }
    }

    fn vertex(&self, id: u64, units: Units) -> Result<P3, Error> {
        let r = self.record(id, "VERTEX_POINT")?;
        self.point(self.reference(id, r.params.get(1))?, units.length_mm)
    }

    fn edge(&self, id: u64, units: Units) -> Result<Edge, Error> {
        let r = self.record(id, "EDGE_CURVE")?;
        let start_vertex = self.reference(id, r.params.get(1))?;
        let end_vertex = self.reference(id, r.params.get(2))?;
        Ok(Edge {
            curve: self.curve(self.reference(id, r.params.get(3))?, units)?,
            start: self.vertex(start_vertex, units)?,
            end: self.vertex(end_vertex, units)?,
            start_vertex,
            end_vertex,
            same_sense: self.flag(id, r.params.get(4))?,
        })
    }
}
