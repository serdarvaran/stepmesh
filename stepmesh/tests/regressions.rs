use stepmesh::{Error, Options, SolidMesh, StepMesh, read_step};

const BOX: &str = include_str!("fixtures/box.step");

#[test]
fn discontinuous_spline_knots_are_rejected() {
    let text = swap(
        BOX,
        "#27 = LINE('',#28,#29);",
        "#27 = B_SPLINE_CURVE_WITH_KNOTS('',1,(#28,#28,#28,#28),.UNSPECIFIED.,.F.,.F.,(2,2,2),(0.,0.5,1.),.UNSPECIFIED.);",
    );
    let x = stepmesh::p21::parse(text.as_bytes()).unwrap();
    let got = stepmesh::internals::decode_curve(
        &x,
        27,
        stepmesh::internals::Units {
            length_mm: 1.,
            angle_rad: 1.,
        },
    );
    assert!(matches!(got, Err(Error::Invalid { .. })));
}

fn swap(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "missing {from}");
    text.replacen(from, to, 1)
}

fn add(text: &str, data: &str) -> String {
    swap(
        text,
        "ENDSEC;\nEND-ISO",
        &format!("{data}\nENDSEC;\nEND-ISO"),
    )
}

fn loose_box() -> String {
    swap(
        BOX,
        "#3 = SHAPE_DEFINITION_REPRESENTATION(#4,#10);",
        "#3 = CARTESIAN_POINT('',(0.,0.,0.));",
    )
}

fn read(text: &str) -> Result<StepMesh, Error> {
    read_step(text.as_bytes(), &Options::default())
}

#[test]
fn zero_degree_bezier_curves_are_errors() {
    let text = swap(
        BOX,
        "#27 = LINE('',#28,#29);",
        "#27 = BEZIER_CURVE('',0,(#28),.UNSPECIFIED.,.F.,.F.);",
    );
    assert!(matches!(read(&text), Err(Error::Invalid { .. })));
}

#[test]
fn zero_degree_bezier_surfaces_are_errors() {
    for (du, dv) in [(0, 0), (0, 1), (1, 0)] {
        let text = swap(
            BOX,
            "#32 = PLANE('',#33);",
            &format!(
                "#32 = BEZIER_SURFACE('',{du},{dv},((#28,#28),(#28,#28)),.UNSPECIFIED.,.F.,.F.,.F.);"
            ),
        );
        assert!(matches!(read(&text), Err(Error::Invalid { .. })));
    }
}

#[test]
fn mesh_volume_is_translation_invariant() {
    let mesh = read(BOX).unwrap();
    let solid = &mesh.solids[0];
    for shift in [0., 1e8, 1e9, 1e10, 1e12] {
        let moved = SolidMesh {
            name: None,
            positions: solid
                .positions
                .iter()
                .map(|p| p.map(|x| x + shift))
                .collect(),
            triangles: solid.triangles.clone(),
        };
        assert!(
            (moved.volume() - 6000.).abs() < 1e-9,
            "shift {shift}: {}",
            moved.volume()
        );
    }
}

fn inch_rep() -> &'static str {
    "#9000=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#15),#9001);\
     #9001=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#9002))REPRESENTATION_CONTEXT('',''));\
     #9002=(LENGTH_UNIT()NAMED_UNIT(*)CONVERSION_BASED_UNIT('inch',#9003));\
     #9003=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#346);"
}

#[test]
fn shared_brep_respects_each_representation_unit() {
    let text = add(&loose_box(), inch_rep());
    let mesh = read(&text).unwrap();
    let mut volumes: Vec<_> = mesh.solids.iter().map(SolidMesh::volume).collect();
    volumes.sort_by(f64::total_cmp);
    assert_eq!(volumes.len(), 2);
    assert!((volumes[0] - 6000.).abs() < 1e-7);
    assert!((volumes[1] - 98_322_384.).abs() < 1e-5, "{volumes:?}");
}

#[test]
fn coincident_assembly_occurrences_remain_distinct() {
    let text = swap(
        include_str!("fixtures/assembly_two_boxes.step"),
        "#380 = ITEM_DEFINED_TRANSFORMATION('','',#11,#19);",
        "#380 = ITEM_DEFINED_TRANSFORMATION('','',#11,#15);",
    );
    let mesh = read(&text).unwrap();
    assert_eq!(mesh.solids.len(), 2);
    assert_eq!(
        mesh.solids.iter().map(SolidMesh::volume).sum::<f64>(),
        12_000.
    );
}

#[test]
fn mapped_placements_match_product_occurrences_one_to_one() {
    let text = swap(
        include_str!("fixtures/assembly_two_boxes.step"),
        "(#11,#15,#19),#23)",
        "(#11,#15,#19,#9000),#23)",
    );
    let text = add(
        &text,
        "#9000=MAPPED_ITEM('',#9001,#15);#9001=REPRESENTATION_MAP(#11,#36);",
    );
    assert_eq!(read(&text).unwrap().solids.len(), 2);
    let coincident = swap(
        &text,
        "#380 = ITEM_DEFINED_TRANSFORMATION('','',#11,#19);",
        "#380 = ITEM_DEFINED_TRANSFORMATION('','',#11,#15);",
    );
    assert_eq!(read(&coincident).unwrap().solids.len(), 2);
    let unplaced = swap(
        &coincident,
        "#372 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#373,#375);",
        "#372=CARTESIAN_POINT('',(0.,0.,0.));",
    );
    let unplaced = swap(
        &unplaced,
        "#378 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#379,#381);",
        "#378=CARTESIAN_POINT('',(0.,0.,0.));",
    );
    assert_eq!(read(&unplaced).unwrap().solids.len(), 2);
}

#[test]
fn unsupported_loose_solid_is_not_silently_skipped() {
    let text = add(
        &loose_box(),
        "#9000=SHAPE_REPRESENTATION('',(#9001),#345);\
         #9001=CSG_SOLID('',#9002);#9002=BLOCK('',#11,10.,20.,30.);",
    );
    assert!(matches!(
        read(&text),
        Err(Error::Unsupported { id: 9001, .. })
    ));
}

fn mapped_box(shift: &str, nested: bool) -> String {
    let text = add(
        &loose_box(),
        &format!(
            "#9001=SHAPE_REPRESENTATION('',(#9002),#345);\
         #9002=MAPPED_ITEM('',#9003,#9004);#9003=REPRESENTATION_MAP(#11,#10);\
         #9004=AXIS2_PLACEMENT_3D('',#9005,$,$);#9005=CARTESIAN_POINT('',({shift},{shift},{shift}));"
        ),
    );
    if nested {
        add(
            &text,
            "#9011=SHAPE_REPRESENTATION('',(#9012),#345);#9012=MAPPED_ITEM('',#9013,#9004);#9013=REPRESENTATION_MAP(#11,#9001);",
        )
    } else {
        text
    }
}

#[test]
fn placed_box_volume_is_translation_invariant() {
    let mesh = read(&mapped_box("1.E9", false)).unwrap();
    assert!((mesh.solids[0].volume() - 6000.).abs() < 1e-9);
}

#[test]
fn overflowing_placement_is_an_error() {
    assert!(matches!(
        read(&mapped_box("1.E308", true)),
        Err(Error::Invalid { .. })
    ));
}

#[test]
fn placement_that_collapses_a_solid_is_an_error() {
    assert!(matches!(
        read(&mapped_box("1.E308", false)),
        Err(Error::Invalid { .. })
    ));
}

#[test]
fn linked_representations_cannot_disagree_on_solid_units() {
    let text = add(
        &add(BOX, inch_rep()),
        "#9010=SHAPE_REPRESENTATION_RELATIONSHIP('','',#10,#9000);",
    );
    assert!(matches!(read(&text), Err(Error::Invalid { .. })));
}

fn assembly_chain(n: usize) -> String {
    let mut data = String::new();
    for i in 0..n {
        let p = 100_000 + 10 * i;
        let rep = if i == 0 {
            "#10".to_owned()
        } else {
            format!("#{}", p + 5)
        };
        data += &format!(
            "#{p}=PRODUCT('p{i}','p{i}','',());\
            #{}=PRODUCT_DEFINITION_FORMATION('','',#{p});\
            #{}=PRODUCT_DEFINITION('','',#{},$);\
            #{}=PRODUCT_DEFINITION_SHAPE('','',#{});\
            #{}=SHAPE_DEFINITION_REPRESENTATION(#{},{rep});\n",
            p + 1,
            p + 2,
            p + 1,
            p + 3,
            p + 2,
            p + 4,
            p + 3
        );
        if i > 0 {
            data += &format!(
                "#{}=SHAPE_REPRESENTATION('',(#11),#345);\
                #{}=NEXT_ASSEMBLY_USAGE_OCCURRENCE('{i}','','',#{},#{},$);\n",
                p + 5,
                p + 6,
                p + 2,
                p - 8
            );
        }
    }
    add(&loose_box(), &data)
}

#[test]
fn assembly_depth_is_independent_of_entity_order() {
    for n in [128, 1000] {
        let text = assembly_chain(n);
        assert!(matches!(read(&text), Err(Error::Invalid { .. })), "{n}");
    }
}

#[test]
fn shallow_assembly_chain_still_reads() {
    let mesh = read(&assembly_chain(64)).unwrap();
    assert_eq!(mesh.solids.len(), 1);
    assert_eq!(mesh.solids[0].volume(), 6000.);
}

fn empty_shared_graph(levels: usize) -> String {
    let text = swap(BOX, "(#11,#15),#345)", "(#11,#15,#20001),#345)");
    let mut data = String::new();
    for i in 0..levels {
        let r = 20_000 + 4 * i;
        let items = if i + 1 == levels {
            "#11".to_owned()
        } else {
            format!("#{},#{}", r + 5, r + 6)
        };
        data += &format!(
            "#{r}=SHAPE_REPRESENTATION('',({items}),#345);\
            #{}=MAPPED_ITEM('',#{},#11);#{}=MAPPED_ITEM('',#{},#11);\
            #{}=REPRESENTATION_MAP(#11,#{r});\n",
            r + 1,
            r + 3,
            r + 2,
            r + 3,
            r + 3
        );
    }
    add(&text, &data)
}

#[test]
fn empty_shared_subassemblies_finish() {
    const CHILD: &str = "STEPMESH_EMPTY_GRAPH_CHILD";
    if std::env::var_os(CHILD).is_some() {
        let mesh = read(&empty_shared_graph(40)).unwrap();
        assert_eq!(mesh.solids.len(), 1);
        assert_eq!(mesh.solids[0].volume(), 6000.);
        return;
    }
    let mut child = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "empty_shared_subassemblies_finish"])
        .env(CHILD, "1")
        .stdout(std::process::Stdio::null())
        .spawn()
        .unwrap();
    let start = std::time::Instant::now();
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            assert!(status.success(), "{status}");
            break;
        }
        if start.elapsed() > std::time::Duration::from_secs(5) {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("empty subassemblies did not finish within 5 s");
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
}

#[test]
fn boundary_outside_requested_surface_tolerance_is_an_error() {
    let text = swap(
        BOX,
        "#34 = CARTESIAN_POINT('',(-5.,-10.,-15.));",
        "#34 = CARTESIAN_POINT('',(-4.95,-10.,-15.));",
    );
    assert!(matches!(read(&text), Err(Error::Projection { .. })));
}

#[test]
fn boundary_within_requested_surface_tolerance_still_reads() {
    let text = swap(
        BOX,
        "#34 = CARTESIAN_POINT('',(-5.,-10.,-15.));",
        "#34 = CARTESIAN_POINT('',(-4.998,-10.,-15.));",
    );
    assert_eq!(read(&text).unwrap().solids[0].volume(), 6000.);
}

#[test]
fn spline_edge_cannot_hide_between_sample_points() {
    use stepmesh::internals::{Curve, Edge, Nurbs1, discretize_edge};
    let ys = [
        0.,
        0.,
        -1. / 224.,
        13. / 1120.,
        -13. / 1120.,
        1. / 224.,
        0.,
        0.,
    ];
    let curve = Curve::Nurbs(Nurbs1 {
        degree: 7,
        knots: [vec![0.; 8], vec![1.; 8]].concat(),
        ctrl: ys
            .iter()
            .enumerate()
            .map(|(i, y)| [-5., -10. + 1000. * y, -15. + 30. * i as f64 / 7., 1.])
            .collect(),
    });
    let edge = Edge {
        start: curve.point(0.),
        end: curve.point(1.),
        curve,
        start_vertex: 1,
        end_vertex: 2,
        same_sense: true,
    };
    let options = Options::default();
    let points = discretize_edge(&edge, &options).unwrap();
    for i in 0..=2000 {
        let p = edge.curve.point(i as f64 / 2000.);
        let distance = points
            .windows(2)
            .map(|w| {
                let d: [f64; 3] = std::array::from_fn(|k| w[1][k] - w[0][k]);
                let q: [f64; 3] = std::array::from_fn(|k| p[k] - w[0][k]);
                let t = (q.iter().zip(d).map(|(x, y)| x * y).sum::<f64>()
                    / d.iter().map(|x| x * x).sum::<f64>())
                .clamp(0., 1.);
                (0..3)
                    .map(|k| (q[k] - t * d[k]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(f64::INFINITY, f64::min);
        assert!(
            distance <= options.linear_deflection,
            "sample {i}: {distance}"
        );
    }
}

#[test]
fn spline_surface_cannot_hide_between_sample_points() {
    let bu = [
        0.,
        0.,
        1. / 3645.,
        -7. / 9720.,
        11. / 9720.,
        -53. / 40824.,
        11. / 9720.,
        -7. / 9720.,
        1. / 3645.,
        0.,
        0.,
    ];
    let bv = [0., 0., 1. / 6., 0., 0.];
    let mut data = String::new();
    let mut rows = Vec::new();
    for (i, u) in bu.iter().enumerate() {
        let mut row = Vec::new();
        for (j, v) in bv.iter().enumerate() {
            let id = 10000 + 10 * i + j;
            data += &format!(
                "#{id}=CARTESIAN_POINT('',({:.17e},{:.17e},{:.17e}));",
                -5. + 1e6 * u * v,
                -10. + 2. * i as f64,
                -15. + 7.5 * j as f64
            );
            row.push(format!("#{id}"));
        }
        rows.push(format!("({})", row.join(",")));
    }
    let text = swap(
        BOX,
        "#32 = PLANE('',#33);",
        &format!(
            "#32 = BEZIER_SURFACE('',10,4,({}),.UNSPECIFIED.,.F.,.F.,.F.);",
            rows.join(",")
        ),
    );
    let mut text = add(&text, &data);
    for line in BOX.lines().filter(|s| s.contains("SURFACE_CURVE('',#")) {
        let (id, rest) = line.split_once(" = SURFACE_CURVE('',#").unwrap();
        let basis = rest.split(',').next().unwrap();
        text = text.replace(&format!(",{id},.T.);"), &format!(",#{basis},.T.);"));
        text = text.replace(&format!(",{id},.F.);"), &format!(",#{basis},.F.);"));
    }
    let mesh = read(&text).unwrap();
    let exact = 109_766_000. / 18_711.;
    assert!(
        (mesh.solids[0].volume() - exact).abs() < 6.,
        "{}",
        mesh.solids[0].volume()
    );
}
