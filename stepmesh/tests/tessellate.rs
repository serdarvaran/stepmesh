use std::f64::consts::{FRAC_1_SQRT_2, TAU};

use stepmesh::Options;
use stepmesh::internals::{Curve, Edge, Frame, Nurbs1, discretize_edge, expand_knots};

const DEFLECTION: f64 = 0.003;

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn xy_frame() -> Frame {
    Frame::new([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}

fn edge(curve: Curve, start: [f64; 3], end: [f64; 3], same_sense: bool, closed: bool) -> Edge {
    Edge {
        curve,
        start,
        end,
        start_vertex: 1,
        end_vertex: if closed { 1 } else { 2 },
        same_sense,
    }
}

fn on_circle(r: f64, t: f64) -> [f64; 3] {
    [r * t.cos(), r * t.sin(), 0.0]
}

/// Angle swept by a polyline on a circle around the origin, signed by direction.
fn swept(points: &[[f64; 3]]) -> f64 {
    points
        .windows(2)
        .map(|w| {
            let (a, b) = (w[0], w[1]);
            (a[0] * b[1] - a[1] * b[0]).atan2(a[0] * b[0] + a[1] * b[1])
        })
        .sum()
}

#[test]
fn line_edge_has_two_points() {
    let e = edge(
        Curve::Line {
            origin: [0.0; 3],
            dir: [1.0, 0.0, 0.0],
        },
        [1.0, 0.0, 0.0],
        [7.0, 0.0, 0.0],
        true,
        false,
    );
    assert_eq!(
        discretize_edge(&e, &Options::default()).unwrap(),
        [[1.0, 0.0, 0.0], [7.0, 0.0, 0.0]]
    );
}

#[test]
fn circle_edge_respects_sagitta_and_angle() {
    let r = 5.0;
    let e = edge(
        Curve::Circle {
            frame: xy_frame(),
            r,
        },
        on_circle(r, 0.0),
        on_circle(r, 2.0),
        true,
        false,
    );
    let pts = discretize_edge(&e, &Options::default()).unwrap();
    assert!((swept(&pts) - 2.0).abs() < 1e-9);
    for w in pts.windows(2) {
        let mid = [(w[0][0] + w[1][0]) / 2.0, (w[0][1] + w[1][1]) / 2.0, 0.0];
        assert!(r - dist(mid, [0.0; 3]) <= DEFLECTION + 1e-12);
        assert!(swept(w) <= 0.5 + 1e-12);
    }
}

#[test]
fn closed_circle_edge_is_a_full_turn() {
    let r = 2.0;
    let e = edge(
        Curve::Circle {
            frame: xy_frame(),
            r,
        },
        on_circle(r, 1.0),
        on_circle(r, 1.0),
        true,
        true,
    );
    let pts = discretize_edge(&e, &Options::default()).unwrap();
    assert_eq!(pts.first(), pts.last());
    assert!((swept(&pts) - TAU).abs() < 1e-9, "{}", swept(&pts));
}

#[test]
fn reversed_sense_goes_the_other_way_round() {
    let r = 3.0;
    let e = edge(
        Curve::Circle {
            frame: xy_frame(),
            r,
        },
        on_circle(r, 0.0),
        on_circle(r, 1.0),
        false,
        false,
    );
    let pts = discretize_edge(&e, &Options::default()).unwrap();
    assert!((swept(&pts) + (TAU - 1.0)).abs() < 1e-9, "{}", swept(&pts));
}

#[test]
fn nurbs_edge_adaptive_deviation() {
    let r = 50.0;
    let w = FRAC_1_SQRT_2;
    let pts = [
        (1., 0., 1.),
        (1., 1., w),
        (0., 1., 1.),
        (-1., 1., w),
        (-1., 0., 1.),
    ];
    let curve = Curve::Nurbs(Nurbs1 {
        degree: 2,
        knots: expand_knots(&[3, 2, 3], &[0.0, 1.0, 2.0]),
        ctrl: pts
            .iter()
            .map(|&(x, y, w)| [x * r * w, y * r * w, 0.0, w])
            .collect(),
    });
    let e = edge(curve, [r, 0.0, 0.0], [-r, 0.0, 0.0], true, false);
    let pts = discretize_edge(&e, &Options::default()).unwrap();
    assert!((swept(&pts) - std::f64::consts::PI).abs() < 1e-9);
    for w in pts.windows(2) {
        let mid = [(w[0][0] + w[1][0]) / 2.0, (w[0][1] + w[1][1]) / 2.0, 0.0];
        assert!(
            r - dist(mid, [0.0; 3]) <= DEFLECTION + 1e-9,
            "{}",
            r - dist(mid, [0.0; 3])
        );
    }
    assert!(pts.len() < 400, "not oversampled: {}", pts.len());
}

#[test]
fn endpoints_are_exact_vertices() {
    let r = 4.0;
    // vertex coordinates a little off the curve, as in real files
    let (s, t) = ([r + 1e-8, 0.0, 0.0], [0.0, r - 1e-8, 0.0]);
    let e = edge(
        Curve::Circle {
            frame: xy_frame(),
            r,
        },
        s,
        t,
        true,
        false,
    );
    let pts = discretize_edge(&e, &Options::default()).unwrap();
    assert_eq!(pts[0], s);
    assert_eq!(*pts.last().unwrap(), t);
}

#[test]
fn fixture_edges_discretize() {
    use stepmesh::internals::{decode_brep, decode_model};
    let dir = format!("{}/tests/fixtures", env!("CARGO_MANIFEST_DIR"));
    let o = Options::default();
    for f in std::fs::read_dir(dir).unwrap() {
        let path = f.unwrap().path();
        if path.extension().is_none_or(|e| e != "step") || path.ends_with("external_reference.step")
        {
            continue;
        }
        let x = stepmesh::p21::parse(&std::fs::read(&path).unwrap()).unwrap();
        let m = decode_model(&x).unwrap();
        for &(id, units) in m.parts.iter().flat_map(|p| &p.breps) {
            let b = decode_brep(&x, id, units).unwrap();
            for (eid, e) in &b.edges {
                let pts = discretize_edge(e, &o).unwrap();
                assert!(pts.len() >= 2, "{path:?} #{eid}");
                for w in pts.windows(2) {
                    assert!(
                        dist(w[0], w[1]) > 0.0 || e.start_vertex == e.end_vertex,
                        "{path:?} #{eid} repeats a point"
                    );
                }
            }
        }
    }
}
