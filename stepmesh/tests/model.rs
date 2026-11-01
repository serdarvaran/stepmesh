use stepmesh::Error;
use stepmesh::internals::{Model, Transform, decode_model};
use stepmesh::p21::parse;

fn fixture(name: &str) -> Vec<u8> {
    std::fs::read(format!(
        "{}/tests/fixtures/{name}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn model(name: &str) -> Result<Model, Error> {
    decode_model(&parse(&fixture(name)).unwrap())
}

/// Every B-rep with its world transform and length scale.
fn placed(m: &Model) -> Vec<(u64, f64, Transform)> {
    fn walk(m: &Model, part: usize, t: &Transform, out: &mut Vec<(u64, f64, Transform)>) {
        let p = &m.parts[part];
        out.extend(
            p.breps
                .iter()
                .map(|&(id, units)| (id, units.length_mm, t.clone())),
        );
        for c in &p.children {
            walk(m, c.part, &t.compose(&c.transform), out);
        }
    }
    let mut out = Vec::new();
    for r in &m.roots {
        walk(m, r.part, &r.transform, &mut out);
    }
    out
}

#[test]
fn box_has_one_root_and_one_brep() {
    let m = model("box.step").unwrap();
    assert_eq!(m.roots.len(), 1);
    let b = placed(&m);
    assert_eq!(b.len(), 1);
    assert_eq!(b[0].1, 1.0);
}

#[test]
fn inch_unit_scale_is_25_4() {
    let b = placed(&model("box_inch.step").unwrap());
    assert_eq!(b.len(), 1);
    assert!((b[0].1 - 25.4).abs() < 1e-12, "{}", b[0].1);
}

#[test]
fn assembly_two_boxes_has_two_instances_50mm_apart() {
    let b = placed(&model("assembly_two_boxes.step").unwrap());
    assert_eq!(b.len(), 2);
    let (p, q) = (b[0].2.apply([0.0; 3]), b[1].2.apply([0.0; 3]));
    let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let d = if d[0] < 0.0 { d.map(|c| -c) } else { d };
    assert!(
        (d[0] - 50.0).abs() < 1e-9 && d[1].abs() < 1e-9 && d[2].abs() < 1e-9,
        "{d:?}"
    );
    assert!(
        p[0].abs() < 1e-9 || q[0].abs() < 1e-9,
        "one box stays at the origin: {p:?} {q:?}"
    );
}

#[test]
fn mixed_unit_assembly_scales() {
    let b = placed(&model("mixed_units.step").unwrap());
    assert_eq!(b.len(), 1);
    assert!(
        (b[0].1 - 25.4).abs() < 1e-12,
        "part is in inches: {}",
        b[0].1
    );
    let o = b[0].2.apply([0.0; 3]);
    assert!(
        (o[0] - 30.0).abs() < 1e-9 && o[1].abs() < 1e-9,
        "placement stays in millimetres: {o:?}"
    );
}

#[test]
fn external_reference_is_no_geometry() {
    assert_eq!(
        model("external_reference.step").err(),
        Some(Error::NoGeometry)
    );
}

const CYCLE: &str = "ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
#1=PRODUCT('a','a','',());#2=PRODUCT_DEFINITION_FORMATION('','',#1);#3=PRODUCT_DEFINITION('','',#2,$);
#4=PRODUCT('b','b','',());#5=PRODUCT_DEFINITION_FORMATION('','',#4);#6=PRODUCT_DEFINITION('','',#5,$);
#7=PRODUCT('c','c','',());#8=PRODUCT_DEFINITION_FORMATION('','',#7);#9=PRODUCT_DEFINITION('','',#8,$);
#10=NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','','',#3,#6,$);
#11=NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','','',#6,#9,$);
#12=NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','','',#9,#6,$);
#20=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#21))REPRESENTATION_CONTEXT('',''));
#21=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#22=AXIS2_PLACEMENT_3D('',#23,$,$);#23=CARTESIAN_POINT('',(0.,0.,0.));
#30=SHAPE_REPRESENTATION('',(#22),#20);#31=PRODUCT_DEFINITION_SHAPE('','',#3);#32=SHAPE_DEFINITION_REPRESENTATION(#31,#30);
#33=SHAPE_REPRESENTATION('',(#22),#20);#34=PRODUCT_DEFINITION_SHAPE('','',#6);#35=SHAPE_DEFINITION_REPRESENTATION(#34,#33);
#36=SHAPE_REPRESENTATION('',(#22),#20);#37=PRODUCT_DEFINITION_SHAPE('','',#9);#38=SHAPE_DEFINITION_REPRESENTATION(#37,#36);
ENDSEC;END-ISO-10303-21;";

#[test]
fn cyclic_reference_is_an_error() {
    let r = decode_model(&parse(CYCLE.as_bytes()).unwrap());
    assert!(matches!(r, Err(Error::Invalid { .. })), "{r:?}");
}

#[test]
fn transform_compose_matches_sequential_apply() {
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    let mut rnd = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    for _ in 0..100 {
        let mut t = || {
            let z = [rnd(), rnd(), rnd() + 1.5];
            let x = [rnd() + 1.5, rnd(), rnd()];
            Transform::from_frame([rnd() * 50.0, rnd() * 50.0, rnd() * 50.0], z, x)
        };
        let (a, b) = (t(), t());
        let p = [rnd() * 10.0, rnd() * 10.0, rnd() * 10.0];
        let (u, v) = (a.compose(&b).apply(p), a.apply(b.apply(p)));
        assert!((0..3).all(|i| (u[i] - v[i]).abs() < 1e-9), "{u:?} {v:?}");
        let w = a.inverse().apply(a.apply(p));
        assert!((0..3).all(|i| (w[i] - p[i]).abs() < 1e-9), "{w:?} {p:?}");
    }
}

mod topology {
    use super::*;
    use stepmesh::internals::{Brep, Loop, Units, decode_brep};
    use stepmesh::p21::Exchange;

    fn first_brep(x: &Exchange) -> Result<Brep, Error> {
        let m = decode_model(x)?;
        let (id, units) = m
            .parts
            .iter()
            .flat_map(|p| &p.breps)
            .next()
            .copied()
            .unwrap();
        decode_brep(x, id, units)
    }

    fn brep(name: &str) -> Brep {
        first_brep(&parse(&fixture(name)).unwrap()).unwrap()
    }

    #[test]
    fn box_topology_counts() {
        let b = brep("box.step");
        assert_eq!(b.shells.len(), 1);
        assert_eq!(b.shells[0].faces.len(), 6);
        assert_eq!(b.edges.len(), 12);
        let mut vertices: Vec<u64> = b
            .edges
            .values()
            .flat_map(|e| [e.start_vertex, e.end_vertex])
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        assert_eq!(vertices.len(), 8);
        // OpenCascade writes FACE_BOUND even for the outer boundary, so `outer` stays None here.
        for f in &b.shells[0].faces {
            assert_eq!(f.loops.len(), 1);
        }
    }

    #[test]
    fn cylinder_seam_edge_used_twice_in_one_loop() {
        let b = brep("cylinder.step");
        let twice = b.shells[0].faces.iter().any(|f| {
            f.loops.iter().any(|(l, _)| match l {
                Loop::Edges(es) => es
                    .iter()
                    .any(|(e, _)| es.iter().filter(|(o, _)| o == e).count() == 2),
                Loop::Vertex(_) => false,
            })
        });
        assert!(twice);
    }

    #[test]
    fn box_with_hole_has_inner_loop() {
        let b = brep("box_with_hole.step");
        assert!(b.shells[0].faces.iter().any(|f| f.loops.len() == 2));
    }

    #[test]
    fn shelled_box_has_void_shell() {
        let b = brep("shelled_box.step");
        assert_eq!(b.shells.len(), 2);
        assert!(!b.shells[0].void);
        assert!(b.shells[1].void);
        assert_eq!(b.shells[1].faces.len(), 6);
    }

    #[test]
    fn edge_points_are_scaled() {
        let x = parse(&fixture("box.step")).unwrap();
        let m = decode_model(&x).unwrap();
        let id = m.parts.iter().flat_map(|p| &p.breps).next().unwrap().0;
        let b = decode_brep(
            &x,
            id,
            Units {
                length_mm: 2.0,
                angle_rad: 1.0,
            },
        )
        .unwrap();
        let max = b
            .edges
            .values()
            .flat_map(|e| [e.start, e.end])
            .map(|p| p[2])
            .fold(f64::MIN, f64::max);
        assert_eq!(max, 30.0, "box is 30 high, ±15 doubled");
    }

    /// box.step with the surface of its first face wrapped in an OFFSET_SURFACE.
    fn box_with_offset_surface() -> String {
        let text = String::from_utf8(fixture("box.step")).unwrap();
        let line = text.lines().find(|l| l.contains("ADVANCED_FACE(")).unwrap();
        let (head, tail) = line.rsplit_once(",#").unwrap();
        let (surface, rest) = tail.split_once(',').unwrap();
        let new_line = format!("{head},#9999,{rest}");
        text.replacen(line, &new_line, 1).replacen(
            "ENDSEC;\nEND-ISO",
            &format!("#9999 = OFFSET_SURFACE('',#{surface},1.,.F.);\nENDSEC;\nEND-ISO"),
            1,
        )
    }

    #[test]
    fn offset_surface_is_unsupported() {
        let x = parse(box_with_offset_surface().as_bytes()).unwrap();
        assert_eq!(
            first_brep(&x).err(),
            Some(Error::Unsupported {
                entity: "OFFSET_SURFACE".into(),
                id: 9999
            })
        );
    }

    #[test]
    fn open_shell_is_rejected() {
        let text = String::from_utf8(fixture("box.step")).unwrap();
        let surface_model = text.replace(
            "MANIFOLD_SOLID_BREP('',#16)",
            "SHELL_BASED_SURFACE_MODEL('',(#16))",
        );
        assert_ne!(surface_model, text);
        assert!(matches!(
            model_of(&surface_model),
            Err(Error::OpenShell { .. })
        ));
        let open = text.replace("= CLOSED_SHELL(", "= OPEN_SHELL(");
        assert_ne!(open, text);
        let x = parse(open.as_bytes()).unwrap();
        assert!(matches!(first_brep(&x), Err(Error::OpenShell { .. })));
    }

    fn model_of(text: &str) -> Result<Model, Error> {
        decode_model(&parse(text.as_bytes()).unwrap())
    }
}
