//! Broken, hostile or unusual input: an explicit error (or the right result), never a panic,
//! hang, abort or silently wrong mesh.

use std::sync::mpsc;
use std::time::Duration;

use stepmesh::internals::{Units, decode_curve, decode_model};
use stepmesh::p21::parse;
use stepmesh::{Error, Options, read_step};

const HEAD: &str = "ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;\n";
const TAIL: &str = "\nENDSEC;END-ISO-10303-21;";
const CONTEXT: &str = "#4=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#5))\
    REPRESENTATION_CONTEXT('',''));#5=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\
    #6=CARTESIAN_POINT('',(0.,0.,0.));";
const MM: Units = Units {
    length_mm: 1.0,
    angle_rad: 1.0,
};

fn fixture(name: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/tests/fixtures/{name}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn add(text: &str, lines: &str) -> String {
    assert!(text.contains("ENDSEC;\nEND-ISO"));
    text.replacen("ENDSEC;\nEND-ISO", &format!("{lines}\nENDSEC;\nEND-ISO"), 1)
}

fn swap(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from}");
    text.replacen(from, to, 1)
}

/// Runs `f` on a thread with the default stack and fails the test if it does not finish.
fn timed<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(f());
    });
    rx.recv_timeout(Duration::from_secs(20))
        .expect("did not finish within 20 s")
}

/// Solid count and the centre x of the first solid.
fn read_with(text: String, o: Options) -> Result<(usize, f64), Error> {
    timed(move || {
        read_step(text.as_bytes(), &o).map(|m| {
            let x = m.solids.first().map_or(f64::NAN, |s| {
                let [lo, hi] = s.bounds();
                (lo[0] + hi[0]) / 2.0
            });
            (m.solids.len(), x)
        })
    })
}

fn read(text: String) -> Result<(usize, f64), Error> {
    read_with(text, Options::default())
}

fn curve(data: &str, id: u64, units: Units) -> Result<(), Error> {
    let text = format!("{HEAD}{data}{TAIL}");
    timed(move || decode_curve(&parse(text.as_bytes()).unwrap(), id, units).map(|_| ()))
}

fn model(text: String) -> Result<(), Error> {
    timed(move || decode_model(&parse(text.as_bytes()).unwrap()).map(|_| ()))
}

#[test]
fn shell_that_contains_itself() {
    let text = format!(
        "{HEAD}#1=MANIFOLD_SOLID_BREP('',#2);#2=ORIENTED_CLOSED_SHELL('',*,#2,.T.);\
         #3=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#1),#4);{CONTEXT}{TAIL}"
    );
    assert!(matches!(read(text), Err(Error::Invalid { .. })));
}

#[test]
fn curve_that_trims_itself() {
    let r = curve("#9=TRIMMED_CURVE('',#9,(),(),.T.,.PARAMETER.);", 9, MM);
    assert!(matches!(r, Err(Error::Invalid { .. })), "{r:?}");
}

/// Product #22 shown by representation #100; `body` adds the representations from there.
fn product_file(body: &str) -> String {
    format!(
        "{HEAD}{CONTEXT}#7=MANIFOLD_SOLID_BREP('',#8);#8=CLOSED_SHELL('',());\
         #20=PRODUCT('p','p','',());#21=PRODUCT_DEFINITION_FORMATION('','',#20);\
         #22=PRODUCT_DEFINITION('','',#21,$);#23=PRODUCT_DEFINITION_SHAPE('','',#22);\
         #24=SHAPE_DEFINITION_REPRESENTATION(#23,#100);\n{body}{TAIL}"
    )
}

#[test]
fn deep_chain_of_mapped_items() {
    let n = 20_000;
    let mut body = String::new();
    for i in 0..n {
        let r = 100 + 4 * i;
        body += &format!(
            "#{r}=SHAPE_REPRESENTATION('',(#{}),#4);#{}=MAPPED_ITEM('',#{},#{});#{}=REPRESENTATION_MAP(#{},#{});\
             #{}=AXIS2_PLACEMENT_3D('',#6,$,$);\n",
            r + 1,
            r + 1,
            r + 2,
            r + 3,
            r + 2,
            r + 3,
            r + 4,
            r + 3
        );
    }
    body += &format!("#{}=SHAPE_REPRESENTATION('',(#7),#4);", 100 + 4 * n);
    assert!(matches!(
        model(product_file(&body)),
        Err(Error::Invalid { .. })
    ));
}

#[test]
fn mapped_items_doubling_at_every_level() {
    let levels = 40;
    let mut body = String::new();
    for i in 0..levels {
        let r = 100 + 6 * i;
        body += &format!(
            "#{r}=SHAPE_REPRESENTATION('',(#{},#{}),#4);#{}=MAPPED_ITEM('',#{},#{});#{}=MAPPED_ITEM('',#{},#{});\
             #{}=REPRESENTATION_MAP(#{},#{});#{}=AXIS2_PLACEMENT_3D('',#6,$,$);\n",
            r + 1,
            r + 2,
            r + 1,
            r + 3,
            r + 4,
            r + 2,
            r + 3,
            r + 4,
            r + 3,
            r + 4,
            r + 6,
            r + 4
        );
    }
    body += &format!("#{}=SHAPE_REPRESENTATION('',(#7),#4);", 100 + 6 * levels);
    assert!(matches!(
        model(product_file(&body)),
        Err(Error::Invalid { .. })
    ));
}

#[test]
fn absurd_knot_multiplicity() {
    let pts = "#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(1.,1.,0.));\
               #3=CARTESIAN_POINT('',(2.,1.,0.));#4=CARTESIAN_POINT('',(3.,0.,0.));";
    let r = curve(
        &format!(
            "{pts}#10=B_SPLINE_CURVE_WITH_KNOTS('',3,(#1,#2,#3,#4),.UNSPECIFIED.,.F.,.F.,\
             (4000000000000000,4),(0.,1.),.UNSPECIFIED.);"
        ),
        10,
        MM,
    );
    assert!(matches!(r, Err(Error::Invalid { .. })), "{r:?}");
}

#[test]
fn unread_representation_item_is_an_error() {
    let text = add(
        &swap(
            &fixture("box.step"),
            "(#11,#15),#345)",
            "(#11,#15,#9001),#345)",
        ),
        "#9001=TESSELLATED_SOLID('',(),$);",
    );
    let r = read(text);
    assert!(
        matches!(&r, Err(Error::Unsupported { entity, .. }) if entity == "TESSELLATED_SOLID"),
        "{r:?}"
    );
}

#[test]
fn assembly_component_without_a_shape_is_an_error() {
    let text = add(
        &fixture("assembly_two_boxes.step"),
        "#9001=PRODUCT('c','c','',());#9002=PRODUCT_DEFINITION_FORMATION('','',#9001);\
         #9003=PRODUCT_DEFINITION('design','',#9002,#35);\
         #9004=NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','c','',#5,#9003,$);",
    );
    assert!(read(text).is_err());
}

/// box.step without product structure, plus a representation that maps the box representation
/// with `target` as its placement.
fn box_mapped_by(target: &str) -> String {
    let text = swap(
        &fixture("box.step"),
        "#3 = SHAPE_DEFINITION_REPRESENTATION(#4,#10);",
        "#3 = CARTESIAN_POINT('',(0.,0.,0.));",
    );
    add(
        &text,
        &format!(
            "#9001=SHAPE_REPRESENTATION('',(#9002),#345);#9002=MAPPED_ITEM('',#9003,#9004);\
             #9003=REPRESENTATION_MAP(#11,#10);#9005=CARTESIAN_POINT('',(100.,0.,0.));\
             #9006=DIRECTION('',(1.,0.,0.));#9007=DIRECTION('',(0.,-1.,0.));#9008=DIRECTION('',(0.,0.,1.));\n{target}"
        ),
    )
}

#[test]
fn mapped_representation_is_not_also_a_root() {
    let r = read(box_mapped_by("#9004=AXIS2_PLACEMENT_3D('',#9005,$,$);"));
    let (n, x) = r.unwrap();
    assert_eq!(n, 1);
    assert!((x - 100.0).abs() < 1e-9, "{x}");
}

#[test]
fn solid_listed_by_two_linked_representations_counts_once() {
    let text = add(
        &fixture("box.step"),
        "#9001=SHAPE_REPRESENTATION('',(#15),#345);#9002=SHAPE_REPRESENTATION_RELATIONSHIP('','',#10,#9001);",
    );
    assert_eq!(read(text).unwrap().0, 1);
}

#[test]
fn cartesian_transformation_operators() {
    let op = |args: &str| {
        box_mapped_by(&format!(
            "#9004=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',{args});"
        ))
    };
    // translation only
    let (n, x) = read(op("$,$,#9005,1.,$")).unwrap();
    assert_eq!(n, 1);
    assert!((x - 100.0).abs() < 1e-9, "{x}");
    // a scale written as an integer is still a scale
    assert!(matches!(
        read(op("$,$,#9005,2,$")),
        Err(Error::Unsupported { .. })
    ));
    // axis2 against axis3 x axis1: a mirror
    assert!(matches!(
        read(op("#9006,#9007,#9005,1.,#9008")),
        Err(Error::Unsupported { .. })
    ));
    // axis1 parallel to axis2 without axis3: no frame
    assert!(matches!(
        read(op("#9006,#9006,#9005,1.,$")),
        Err(Error::Invalid { .. })
    ));
}

const PTS: &str = "#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(1.,1.,0.));\
                   #3=CARTESIAN_POINT('',(2.,1.,0.));#4=CARTESIAN_POINT('',(3.,0.,0.));";

#[test]
fn invalid_bsplines_are_errors() {
    for (what, data) in [
        (
            "zero weights",
            "#10=(BOUNDED_CURVE()B_SPLINE_CURVE(3,(#1,#2,#3,#4),.UNSPECIFIED.,.F.,.F.)\
            B_SPLINE_CURVE_WITH_KNOTS((4,4),(0.,1.),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()\
            RATIONAL_B_SPLINE_CURVE((0.,0.,0.,0.))REPRESENTATION_ITEM(''));",
        ),
        (
            "equal knots",
            "#10=B_SPLINE_CURVE_WITH_KNOTS('',3,(#1,#2,#3,#4),.UNSPECIFIED.,.F.,.F.,(8),(1.),.UNSPECIFIED.);",
        ),
        (
            "no control points",
            "#10=B_SPLINE_CURVE_WITH_KNOTS('',3,(),.UNSPECIFIED.,.F.,.F.,(),(),.UNSPECIFIED.);",
        ),
        (
            "too few points",
            "#10=B_SPLINE_CURVE_WITH_KNOTS('',3,(#1,#2),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);",
        ),
        (
            "decreasing knots",
            "#10=B_SPLINE_CURVE_WITH_KNOTS('',3,(#1,#2,#3,#4),.UNSPECIFIED.,.F.,.F.,(4,4),(1.,0.),.UNSPECIFIED.);",
        ),
    ] {
        let r = curve(&format!("{PTS}{data}"), 10, MM);
        assert!(matches!(r, Err(Error::Invalid { .. })), "{what}: {r:?}");
    }
    let huge = "#5=AXIS2_PLACEMENT_3D('',#1,$,$);#10=CIRCLE('',#5,1.E306);";
    let r = curve(
        &format!("{PTS}{huge}"),
        10,
        Units {
            length_mm: 1000.0,
            angle_rad: 1.0,
        },
    );
    assert!(matches!(r, Err(Error::Invalid { .. })), "{r:?}");
}

#[test]
fn invalid_options_are_errors() {
    for o in [
        Options {
            linear_deflection: 0.0,
            ..Options::default()
        },
        Options {
            linear_deflection: f64::NAN,
            ..Options::default()
        },
        Options {
            angular_deflection: 0.0,
            ..Options::default()
        },
        Options {
            angular_deflection: f64::INFINITY,
            ..Options::default()
        },
    ] {
        assert!(matches!(
            read_with(fixture("box.step"), o),
            Err(Error::InvalidOptions(_))
        ));
    }
}

#[test]
fn triangle_limit_holds_inside_a_face() {
    let o = Options {
        linear_deflection: 1e-7,
        max_triangles: Some(1000),
        ..Options::default()
    };
    assert!(matches!(
        read_with(fixture("sphere.step"), o),
        Err(Error::TriangleLimit { .. })
    ));
}

#[test]
fn annotation_items_are_left_out() {
    let text = add(
        &swap(
            &fixture("box.step"),
            "(#11,#15),#345)",
            "(#11,#15,#9001,#9002),#345)",
        ),
        "#9001=TEXT_LITERAL('','A',#11,'',.LEFT.);#9002=LEADER_DIRECTED_CALLOUT('',());",
    );
    assert_eq!(read(text).unwrap().0, 1);
}
