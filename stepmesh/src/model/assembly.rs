//! Product structure: which solids each product holds and where its sub-assemblies sit.

use std::collections::{BTreeMap, HashMap, HashSet};

use super::{Decoder, Instance, Model, Part, Transform, Units};
use crate::Error;
use crate::p21::{Exchange, Value};

/// Deeper assembly trees than this are treated as malformed.
const MAX_ASSEMBLY_DEPTH: usize = 64;
/// More placed solids than this (instances counted) are treated as malformed: a few levels of
/// shared sub-assemblies can multiply into billions.
const MAX_PLACED_SOLIDS: u64 = 100_000;

/// Representation items that carry no solid geometry and are left out on purpose.
const PASSIVE_ITEMS: [&str; 9] = [
    "AXIS2_PLACEMENT_3D",
    "AXIS2_PLACEMENT_2D",
    "AXIS1_PLACEMENT",
    "CARTESIAN_POINT",
    "DIRECTION",
    "VECTOR",
    "GEOMETRIC_CURVE_SET",
    "GEOMETRIC_SET",
    "VERTEX_POINT",
];
/// Name parts of wireframe and presentation items (annotations, PMI), also left out.
const PASSIVE_KINDS: [&str; 12] = [
    "CURVE",
    "LINE",
    "CIRCLE",
    "ELLIPSE",
    "POINT",
    "PLACEMENT",
    "ANNOTATION",
    "DRAUGHTING",
    "CALLOUT",
    "LEADER",
    "TEXT",
    "DIMENSION",
];

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
    /// Parts being read, innermost last: guards the recursion through mapped items.
    building: Vec<Key>,
    /// Assembly usages (NEXT_ASSEMBLY_USAGE_OCCURRENCE id, child part), checked for a shape.
    uses: Vec<(u64, usize)>,
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
        building: Vec::new(),
        uses: Vec::new(),
    };
    b.collect_shape_definitions()?;
    b.collect_links()?;
    let roots = if b.product_reps.is_empty() {
        b.loose_roots()?
    } else {
        b.product_roots()?
    };
    let mut model = Model {
        parts: b.parts,
        roots,
    };
    let counts = solid_counts(&model)?;
    for &(nauo, part) in &b.uses {
        if counts[part] == 0 {
            return Err(Error::Unsupported {
                entity: "assembly component without a shape".into(),
                id: nauo,
            });
        }
    }
    let solids = model
        .roots
        .iter()
        .fold(0u64, |a, r| a.saturating_add(counts[r.part]));
    if solids > MAX_PLACED_SOLIDS {
        return Err(Error::Invalid {
            id: 0,
            msg: format!("more than {MAX_PLACED_SOLIDS} placed solids"),
        });
    }
    if solids == 0 {
        return match b.surfaces.first() {
            Some(&id) => Err(Error::OpenShell { id }),
            None => Err(Error::NoGeometry),
        };
    }
    for part in &mut model.parts {
        part.children.retain(|child| counts[child.part] != 0);
    }
    model.roots.retain(|root| counts[root.part] != 0);
    Ok(model)
}

/// Placed solid count of every part with its sub-parts (instances counted, saturating); an
/// assembly that contains itself or nests too deeply is an error.
fn solid_counts(m: &Model) -> Result<Vec<u64>, Error> {
    fn visit(
        m: &Model,
        part: usize,
        depth: usize,
        state: &mut [u8],
        counts: &mut [u64],
        heights: &mut [usize],
    ) -> Result<u64, Error> {
        match state[part] {
            2 => return Ok(counts[part]),
            1 => {
                return Err(Error::Invalid {
                    id: 0,
                    msg: "assembly contains itself".into(),
                });
            }
            _ if depth > MAX_ASSEMBLY_DEPTH => {
                return Err(Error::Invalid {
                    id: 0,
                    msg: "assembly nested too deeply".into(),
                });
            }
            _ => {}
        }
        state[part] = 1;
        let mut n = m.parts[part].breps.len() as u64;
        let mut height = 0;
        for c in &m.parts[part].children {
            n = n.saturating_add(visit(m, c.part, depth + 1, state, counts, heights)?);
            height = height.max(heights[c.part] + 1);
            if height > MAX_ASSEMBLY_DEPTH {
                return Err(Error::Invalid {
                    id: 0,
                    msg: "assembly nested too deeply".into(),
                });
            }
        }
        state[part] = 2;
        counts[part] = n;
        heights[part] = height;
        Ok(n)
    }
    let mut state = vec![0u8; m.parts.len()];
    let mut counts = vec![0u64; m.parts.len()];
    let mut heights = vec![0usize; m.parts.len()];
    for p in 0..m.parts.len() {
        visit(m, p, 0, &mut state, &mut counts, &mut heights)?;
    }
    Ok(counts)
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
        if self.building.len() > MAX_ASSEMBLY_DEPTH {
            return self.d.invalid(0, "representations mapped too deeply");
        }
        let i = self.parts.len();
        self.index.insert(key, i);
        self.parts.push(Part::default());
        self.building.push(key);
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
        self.building.pop();
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
        let units = self.d.units(rep)?;
        for item in self.d.references(rep, r.params.get(1))? {
            let e = self.d.entity(item)?;
            let has = |names: &[&str]| {
                e.records
                    .iter()
                    .find(|r| names.contains(&r.name.as_str()))
                    .map(|r| r.name.clone())
            };
            let name = e.records[0].name.as_str();
            if has(&SOLID_ITEMS).is_some() {
                // A solid listed by two linked representations is one solid.
                if let Some(&(_, previous)) =
                    self.parts[part].breps.iter().find(|&&(b, _)| b == item)
                {
                    if previous != units {
                        return self
                            .d
                            .invalid(item, "linked representations disagree on solid units");
                    }
                } else {
                    self.parts[part].breps.push((item, units));
                }
            } else if has(&SURFACE_ITEMS).is_some() {
                self.surfaces.push(item);
            } else if let Some(name) = has(&UNSUPPORTED_SOLIDS) {
                return Err(Error::Unsupported {
                    entity: name,
                    id: item,
                });
            } else if name == "MAPPED_ITEM" {
                let child = self.mapped_item(item, units)?;
                self.parts[part].children.push(child);
            } else if !(PASSIVE_ITEMS.contains(&name)
                || PASSIVE_KINDS.iter().any(|k| name.contains(k)))
            {
                return Err(Error::Unsupported {
                    entity: name.to_string(),
                    id: item,
                });
            }
        }
        Ok(())
    }

    /// MAPPED_ITEM: the mapped representation placed by `target ∘ origin⁻¹`.
    fn mapped_item(&mut self, item: u64, units: Units) -> Result<Instance, Error> {
        let r = self.d.first(item)?;
        let map = self.d.reference(item, r.params.get(1))?;
        let target = self.d.reference(item, r.params.get(2))?;
        let m = self.d.record(map, "REPRESENTATION_MAP")?;
        let origin = self.d.reference(map, m.params.first())?;
        let mapped = self.d.reference(map, m.params.get(1))?;
        let mapped_scale = self.d.units(mapped)?.length_mm;
        let key = match self.rep_owner.get(&mapped) {
            Some(&pd) => Key::Product(pd),
            None => Key::Representation(mapped),
        };
        let part = self.part(key)?;
        let t = self
            .placement(target, units.length_mm)?
            .compose(&self.placement(origin, mapped_scale)?.inverse());
        Ok(Instance { part, transform: t })
    }

    /// AXIS2_PLACEMENT_3D or a non-scaling CARTESIAN_TRANSFORMATION_OPERATOR_3D.
    fn placement(&self, id: u64, scale: f64) -> Result<Transform, Error> {
        let e = self.d.entity(id)?;
        if e.records.iter().any(|r| r.name == "AXIS2_PLACEMENT_3D") {
            return Ok(Transform::from(&self.d.frame(id, scale)?));
        }
        // CARTESIAN_TRANSFORMATION_OPERATOR(_3D): (axis1, axis2, local_origin, scale[, axis3]),
        // after name and description on a simple instance.
        let (params, axis3): (Vec<&Value>, Option<&Value>) = match e
            .records
            .iter()
            .find(|r| r.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        {
            Some(op) => (
                op.params.iter().collect(),
                e.records
                    .iter()
                    .find(|r| r.name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D")
                    .and_then(|r| r.params.first()),
            ),
            None if e.records[0].name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D" => {
                let p = &e.records[0].params;
                (p.iter().skip(2).collect(), p.get(6))
            }
            None => {
                return Err(Error::Unsupported {
                    entity: e.records[0].name.clone(),
                    id,
                });
            }
        };
        if let Some(s) = params.get(3).filter(|v| !matches!(v, Value::Unset)) {
            if (self.d.real(id, Some(s))? - 1.0).abs() > 1e-12 {
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
        let origin = self
            .d
            .point(self.d.reference(id, params.get(2).copied())?, scale)?;
        // ISO 10303-42 base_axis: z from axis3, x from axis1 made orthogonal to z, y from axis2
        // made orthogonal to both (or z × x).
        let z = dir(axis3)?.unwrap_or([0.0, 0.0, 1.0]);
        let x0 = dir(params.first().copied())?.unwrap_or([1.0, 0.0, 0.0]);
        let Some(x) = super::unit(super::sub(x0, super::mul(z, super::dot(x0, z)))) else {
            return self.d.invalid(id, "axis1 parallel to axis3");
        };
        let natural = super::cross(z, x);
        if let Some(y0) = dir(params.get(1).copied())? {
            let y = super::sub(
                super::sub(y0, super::mul(x, super::dot(y0, x))),
                super::mul(z, super::dot(y0, z)),
            );
            let Some(y) = super::unit(y) else {
                return self.d.invalid(id, "axis2 parallel to the other axes");
            };
            if super::dot(y, natural) < 0.0 {
                return Err(Error::Unsupported {
                    entity: "mirroring CARTESIAN_TRANSFORMATION_OPERATOR".into(),
                    id,
                });
            }
        }
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
        let mut remaining_mapped: HashMap<usize, Vec<Instance>> = HashMap::new();
        for (nauo, parent, child) in uses {
            let p = self.part(Key::Product(parent))?;
            let c = self.part(Key::Product(child))?;
            let transform = placements
                .get(&nauo)
                .map(|&rel| self.relationship_transform(rel, child))
                .transpose()?;
            self.uses.push((nauo, c));
            let mapped = remaining_mapped
                .entry(p)
                .or_insert_with(|| self.parts[p].children.clone());
            // A mapped placement can account for only one product occurrence.
            let same = |i: &Instance| {
                i.part == c
                    && transform.as_ref().is_none_or(|transform| {
                        i.transform
                            .m
                            .iter()
                            .flatten()
                            .zip(transform.m.iter().flatten())
                            .all(|(a, b)| (a - b).abs() < 1e-9)
                    })
            };
            if let Some(i) = mapped.iter().position(same) {
                mapped.swap_remove(i);
            } else {
                self.parts[p].children.push(Instance {
                    part: c,
                    transform: transform.unwrap_or_else(Transform::identity),
                });
            }
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
            let f1 = self.placement(i1, self.d.units(rep1)?.length_mm)?;
            let f2 = self.placement(i2, self.d.units(rep2)?.length_mm)?;
            f2.compose(&f1.inverse())
        } else {
            self.placement(op, self.d.units(rep2)?.length_mm)?
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
        // A representation placed by a MAPPED_ITEM appears where it is mapped, not on its own.
        let mut mapped = HashSet::new();
        for map in self.d.all("REPRESENTATION_MAP") {
            mapped.insert(self.d.reference(map, self.d.first(map)?.params.get(1))?);
        }
        let mut reps: Vec<u64> = self
            .d
            .x
            .entities
            .iter()
            .filter(|(_, e)| e.records[0].name.ends_with("REPRESENTATION"))
            .map(|(id, _)| *id)
            .collect();
        reps.sort_unstable();
        for rep in reps.into_iter().filter(|r| !mapped.contains(r)) {
            if matches!(self.d.first(rep)?.params.get(1), Some(Value::List(_))) {
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
