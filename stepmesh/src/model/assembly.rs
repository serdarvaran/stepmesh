//! Product structure: which solids each product holds and where its sub-assemblies sit.

use std::collections::{BTreeMap, HashMap, HashSet};

use super::{Decoder, Instance, Model, Part, Transform};
use crate::Error;
use crate::p21::{Exchange, Value};

/// Deeper assembly trees than this are treated as malformed.
const MAX_ASSEMBLY_DEPTH: usize = 64;

const SOLID_ITEMS: [&str; 2] = ["MANIFOLD_SOLID_BREP", "BREP_WITH_VOIDS"];
const SURFACE_ITEMS: [&str; 2] = ["SHELL_BASED_SURFACE_MODEL", "FACE_BASED_SURFACE_MODEL"];
const UNSUPPORTED_SOLIDS: [&str; 6] = [
    "FACETED_BREP",
    "EXTRUDED_AREA_SOLID",
    "REVOLVED_AREA_SOLID",
    "CSG_SOLID",
    "EXTRUDED_FACE_SOLID",
    "SWEPT_DISK_SOLID",
];

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Key {
    Product(u64),
    Representation(u64),
}

struct Builder<'a> {
    d: Decoder<'a>,
    /// Representations of each product definition, from SHAPE_DEFINITION_REPRESENTATION.
    product_reps: BTreeMap<u64, Vec<u64>>,
    /// Owning product definition of each such representation.
    rep_owner: HashMap<u64, u64>,
    /// Representation pairs joined without a transformation (same coordinate space).
    links: HashMap<u64, Vec<u64>>,
    parts: Vec<Part>,
    index: HashMap<Key, usize>,
    surfaces: Vec<u64>,
}

pub fn decode_model(x: &Exchange) -> Result<Model, Error> {
    let mut b = Builder {
        d: Decoder::new(x),
        product_reps: BTreeMap::new(),
        rep_owner: HashMap::new(),
        links: HashMap::new(),
        parts: Vec::new(),
        index: HashMap::new(),
        surfaces: Vec::new(),
    };
    b.collect_shape_definitions()?;
    b.collect_links()?;
    let roots = if b.product_reps.is_empty() {
        b.loose_roots()?
    } else {
        b.product_roots()?
    };
    let model = Model {
        parts: b.parts,
        roots,
    };
    let mut solids = 0;
    for r in &model.roots {
        solids += count_solids(&model, r.part, &mut Vec::new())?;
    }
    if solids == 0 {
        return match b.surfaces.first() {
            Some(&id) => Err(Error::OpenShell { id }),
            None => Err(Error::NoGeometry),
        };
    }
    Ok(model)
}

/// Solid count of a part with its sub-parts; an assembly that contains itself is an error.
fn count_solids(m: &Model, part: usize, stack: &mut Vec<usize>) -> Result<usize, Error> {
    if stack.contains(&part) || stack.len() > MAX_ASSEMBLY_DEPTH {
        return Err(Error::Invalid {
            id: 0,
            msg: "assembly contains itself or is nested too deeply".into(),
        });
    }
    stack.push(part);
    let mut n = m.parts[part].breps.len();
    for c in &m.parts[part].children {
        n += count_solids(m, c.part, stack)?;
    }
    stack.pop();
    Ok(n)
}

impl Builder<'_> {
    fn collect_shape_definitions(&mut self) -> Result<(), Error> {
        for sdr in self.d.all("SHAPE_DEFINITION_REPRESENTATION") {
            let r = self.d.record(sdr, "SHAPE_DEFINITION_REPRESENTATION")?;
            let pds = self.d.reference(sdr, r.params.first())?;
            let rep = self.d.reference(sdr, r.params.get(1))?;
            let def = self.d.reference(pds, self.d.first(pds)?.params.get(2))?;
            if self
                .d
                .entity(def)?
                .records
                .iter()
                .any(|r| r.name == "PRODUCT_DEFINITION")
            {
                self.product_reps.entry(def).or_default().push(rep);
                self.rep_owner.insert(rep, def);
            }
        }
        Ok(())
    }

    /// (rep_1, rep_2) of a representation relationship, simple or complex.
    fn related_reps(&self, id: u64) -> Result<(u64, u64), Error> {
        let e = self.d.entity(id)?;
        let Some(r) = e
            .records
            .iter()
            .find(|r| r.name.ends_with("REPRESENTATION_RELATIONSHIP") && r.params.len() >= 4)
        else {
            return self
                .d
                .invalid(id, "representation relationship without representations");
        };
        Ok((
            self.d.reference(id, r.params.get(2))?,
            self.d.reference(id, r.params.get(3))?,
        ))
    }

    fn collect_links(&mut self) -> Result<(), Error> {
        let mut ids = self.d.all("SHAPE_REPRESENTATION_RELATIONSHIP");
        ids.extend(self.d.all("REPRESENTATION_RELATIONSHIP"));
        ids.sort_unstable();
        ids.dedup();
        for id in ids {
            let e = self.d.entity(id)?;
            if e.records
                .iter()
                .any(|r| r.name == "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
            {
                continue;
            }
            let (a, b) = self.related_reps(id)?;
            self.links.entry(a).or_default().push(b);
            self.links.entry(b).or_default().push(a);
        }
        Ok(())
    }

    /// Representations reachable from `start` without a transformation and without entering
    /// a representation that belongs to another product.
    fn closure(&self, start: &[u64], owner: Option<u64>) -> Vec<u64> {
        let mut seen: HashSet<u64> = start.iter().copied().collect();
        let mut todo = start.to_vec();
        let mut out = Vec::new();
        while let Some(r) = todo.pop() {
            out.push(r);
            for &n in self.links.get(&r).into_iter().flatten() {
                let foreign = matches!((self.rep_owner.get(&n), owner), (Some(o), Some(p)) if *o != p)
                    || (owner.is_none() && self.rep_owner.contains_key(&n));
                if !foreign && seen.insert(n) {
                    todo.push(n);
                }
            }
        }
        out.sort_unstable();
        out
    }

    fn part(&mut self, key: Key) -> Result<usize, Error> {
        if let Some(&i) = self.index.get(&key) {
            return Ok(i);
        }
        let i = self.parts.len();
        self.index.insert(key, i);
        self.parts.push(Part::default());
        let (reps, owner) = match key {
            Key::Product(pd) => {
                self.parts[i].name = self.product_name(pd);
                (
                    self.product_reps.get(&pd).cloned().unwrap_or_default(),
                    Some(pd),
                )
            }
            Key::Representation(rep) => (vec![rep], None),
        };
        for rep in self.closure(&reps, owner) {
            self.add_items(i, rep)?;
        }
        Ok(i)
    }

    fn product_name(&self, pd: u64) -> Option<String> {
        let pdf = self
            .d
            .reference(pd, self.d.first(pd).ok()?.params.get(2))
            .ok()?;
        let product = self
            .d
            .reference(pdf, self.d.first(pdf).ok()?.params.get(2))
            .ok()?;
        let p = self.d.first(product).ok()?;
        self.d
            .string(p.params.get(1))
            .or_else(|| self.d.string(p.params.first()))
    }

    fn add_items(&mut self, part: usize, rep: u64) -> Result<(), Error> {
        let r = self.d.first(rep)?;
        let scale = self.d.length_scale_mm(rep)?;
        for item in self.d.references(rep, r.params.get(1))? {
            let e = self.d.entity(item)?;
            let has = |names: &[&str]| {
                e.records
                    .iter()
                    .find(|r| names.contains(&r.name.as_str()))
                    .map(|r| r.name.clone())
            };
            if has(&SOLID_ITEMS).is_some() {
                self.parts[part].breps.push((item, scale));
            } else if has(&SURFACE_ITEMS).is_some() {
                self.surfaces.push(item);
            } else if let Some(name) = has(&UNSUPPORTED_SOLIDS) {
                return Err(Error::Unsupported {
                    entity: name,
                    id: item,
                });
            } else if e.records[0].name == "MAPPED_ITEM" {
                let child = self.mapped_item(item, scale)?;
                self.parts[part].children.push(child);
            }
        }
        Ok(())
    }

    /// MAPPED_ITEM: the mapped representation placed by `target ∘ origin⁻¹`.
    fn mapped_item(&mut self, item: u64, scale: f64) -> Result<Instance, Error> {
        let r = self.d.first(item)?;
        let map = self.d.reference(item, r.params.get(1))?;
        let target = self.d.reference(item, r.params.get(2))?;
        let m = self.d.record(map, "REPRESENTATION_MAP")?;
        let origin = self.d.reference(map, m.params.first())?;
        let mapped = self.d.reference(map, m.params.get(1))?;
        let mapped_scale = self.d.length_scale_mm(mapped)?;
        let key = match self.rep_owner.get(&mapped) {
            Some(&pd) => Key::Product(pd),
            None => Key::Representation(mapped),
        };
        let part = self.part(key)?;
        let t = self
            .placement(target, scale)?
            .compose(&self.placement(origin, mapped_scale)?.inverse());
        Ok(Instance { part, transform: t })
    }

    /// AXIS2_PLACEMENT_3D or a non-scaling CARTESIAN_TRANSFORMATION_OPERATOR_3D.
    fn placement(&self, id: u64, scale: f64) -> Result<Transform, Error> {
        let e = self.d.entity(id)?;
        if e.records.iter().any(|r| r.name == "AXIS2_PLACEMENT_3D") {
            return Ok(Transform::from(&self.d.frame(id, scale)?));
        }
        let Some(op) = e
            .records
            .iter()
            .find(|r| r.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        else {
            return Err(Error::Unsupported {
                entity: e.records[0].name.clone(),
                id,
            });
        };
        if let Some(Value::Real(s)) = op.params.get(3) {
            if (s - 1.0).abs() > 1e-12 {
                return Err(Error::Unsupported {
                    entity: "scaling CARTESIAN_TRANSFORMATION_OPERATOR".into(),
                    id,
                });
            }
        }
        let dir = |v: Option<&Value>| -> Result<Option<[f64; 3]>, Error> {
            self.d
                .optional_reference(id, v)?
                .map(|a| self.d.direction(a))
                .transpose()
        };
        let x = dir(op.params.first())?.unwrap_or([1.0, 0.0, 0.0]);
        let y = dir(op.params.get(1))?;
        let origin = self
            .d
            .point(self.d.reference(id, op.params.get(2))?, scale)?;
        let z = match e
            .records
            .iter()
            .find(|r| r.name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D")
        {
            Some(r3) => dir(r3.params.first())?,
            None => None,
        };
        let z = z
            .or_else(|| y.map(|y| super::cross(x, y)))
            .unwrap_or([0.0, 0.0, 1.0]);
        Ok(Transform::from_frame(origin, z, x))
    }

    fn product_roots(&mut self) -> Result<Vec<Instance>, Error> {
        let mut children: HashSet<u64> = HashSet::new();
        let mut uses = Vec::new();
        for nauo in self.d.all("NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            let r = self.d.first(nauo)?;
            let parent = self.d.reference(nauo, r.params.get(3))?;
            let child = self.d.reference(nauo, r.params.get(4))?;
            children.insert(child);
            uses.push((nauo, parent, child));
        }
        let placements = self.occurrence_placements()?;
        for (nauo, parent, child) in uses {
            let p = self.part(Key::Product(parent))?;
            let c = self.part(Key::Product(child))?;
            let transform = match placements.get(&nauo) {
                Some(&rel) => self.relationship_transform(rel, child)?,
                // Placed by a MAPPED_ITEM of the parent instead, or not at all.
                None if self.parts[p].children.iter().any(|i| i.part == c) => continue,
                None => Transform::identity(),
            };
            self.parts[p].children.push(Instance { part: c, transform });
        }
        let mut roots = Vec::new();
        let pds: Vec<u64> = self.product_reps.keys().copied().collect();
        for pd in pds {
            if !children.contains(&pd) {
                let part = self.part(Key::Product(pd))?;
                roots.push(Instance {
                    part,
                    transform: Transform::identity(),
                });
            }
        }
        Ok(roots)
    }

    /// NAUO id → its CONTEXT_DEPENDENT_SHAPE_REPRESENTATION's representation relationship.
    fn occurrence_placements(&self) -> Result<HashMap<u64, u64>, Error> {
        let mut out = HashMap::new();
        for cdsr in self.d.all("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION") {
            let r = self.d.first(cdsr)?;
            let rel = self.d.reference(cdsr, r.params.first())?;
            let pds = self.d.reference(cdsr, r.params.get(1))?;
            let nauo = self.d.reference(pds, self.d.first(pds)?.params.get(2))?;
            out.insert(nauo, rel);
        }
        Ok(out)
    }

    /// Child-to-parent transform of a REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION.
    /// ISO 10303-42: `transform_item_1` lies in `rep_1`, `transform_item_2` in `rep_2`, and the
    /// operator maps `rep_1` into `rep_2`.
    fn relationship_transform(&self, rel: u64, child: u64) -> Result<Transform, Error> {
        let (rep1, rep2) = self.related_reps(rel)?;
        let e = self.d.entity(rel)?;
        let op = e
            .records
            .iter()
            .find(|r| r.name == "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
            .map_or_else(
                || self.d.invalid(rel, "relationship without transformation"),
                Ok,
            )?;
        let op = self.d.reference(rel, op.params.first())?;
        let t12 = if self
            .d
            .entity(op)?
            .records
            .iter()
            .any(|r| r.name == "ITEM_DEFINED_TRANSFORMATION")
        {
            let r = self.d.record(op, "ITEM_DEFINED_TRANSFORMATION")?;
            let i1 = self.d.reference(op, r.params.get(2))?;
            let i2 = self.d.reference(op, r.params.get(3))?;
            let f1 = self.placement(i1, self.d.length_scale_mm(rep1)?)?;
            let f2 = self.placement(i2, self.d.length_scale_mm(rep2)?)?;
            f2.compose(&f1.inverse())
        } else {
            self.placement(op, self.d.length_scale_mm(rep2)?)?
        };
        let child_reps = self.closure(
            self.product_reps.get(&child).map_or(&[][..], Vec::as_slice),
            Some(child),
        );
        Ok(
            if child_reps.contains(&rep1) || !child_reps.contains(&rep2) {
                t12
            } else {
                t12.inverse()
            },
        )
    }

    /// Files without product structure: every representation that holds a solid is a root.
    fn loose_roots(&mut self) -> Result<Vec<Instance>, Error> {
        let mut roots = Vec::new();
        let mut reps: Vec<u64> = self
            .d
            .x
            .entities
            .iter()
            .filter(|(_, e)| e.records[0].name.ends_with("REPRESENTATION"))
            .map(|(id, _)| *id)
            .collect();
        reps.sort_unstable();
        for rep in reps {
            let holds_solid = self
                .d
                .references(rep, self.d.first(rep)?.params.get(1))
                .is_ok_and(|items| {
                    items.iter().any(|&i| {
                        self.d.entity(i).is_ok_and(|e| {
                            e.records
                                .iter()
                                .any(|r| SOLID_ITEMS.contains(&r.name.as_str()))
                        })
                    })
                });
            if holds_solid {
                let part = self.part(Key::Representation(rep))?;
                roots.push(Instance {
                    part,
                    transform: Transform::identity(),
                });
            }
        }
        Ok(roots)
    }
}
