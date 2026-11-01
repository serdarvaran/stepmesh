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

mod faces {
    use super::*;
    use std::collections::HashMap;
    use stepmesh::internals::{
        Brep, FaceMesh, Surface, VRef, decode_brep, decode_model, tessellate_face,
    };

    struct Meshed {
        brep: Brep,
        faces: Vec<(usize, usize, FaceMesh)>,
        polys: HashMap<u64, Vec<[f64; 3]>>,
    }

    fn fixture_area(name: &str) -> f64 {
        let text = std::fs::read_to_string(format!(
            "{}/tests/fixtures/fixtures.json",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let start = text.find(&format!("\"{name}\"")).unwrap();
        let area = &text[start..][text[start..].find("\"area\": ").unwrap() + 8..];
        area[..area.find(',').unwrap()].trim().parse().unwrap()
    }

    fn mesh(name: &str) -> Vec<Meshed> {
        let x = stepmesh::p21::parse(
            &std::fs::read(format!(
                "{}/tests/fixtures/{name}",
                env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap(),
        )
        .unwrap();
        let m = decode_model(&x).unwrap();
        let o = Options::default();
        m.parts
            .iter()
            .flat_map(|p| &p.breps)
            .map(|&(id, units)| {
                let brep = decode_brep(&x, id, units).unwrap();
                let polys: HashMap<u64, Vec<[f64; 3]>> = brep
                    .edges
                    .iter()
                    .map(|(&k, e)| (k, discretize_edge(e, &o).unwrap()))
                    .collect();
                let mut faces = Vec::new();
                for (si, s) in brep.shells.iter().enumerate() {
                    for (fi, f) in s.faces.iter().enumerate() {
                        let fm = tessellate_face(f, &polys, &o)
                            .unwrap_or_else(|e| panic!("{name} face #{}: {e}", f.id));
                        faces.push((si, fi, fm));
                    }
                }
                Meshed { brep, faces, polys }
            })
            .collect()
    }

    fn point(m: &Meshed, fm: &FaceMesh, r: VRef) -> [f64; 3] {
        match r {
            VRef::Edge { edge, i } => m.polys[&edge][i],
            VRef::Inner(i) => fm.inner[i],
        }
    }

    fn tri_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let (u, v) = (
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
            [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
        );
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        0.5 * dist(n, [0.0; 3])
    }

    fn total_area(ms: &[Meshed]) -> f64 {
        ms.iter()
            .flat_map(|m| m.faces.iter().map(move |(_, _, fm)| (m, fm)))
            .flat_map(|(m, fm)| {
                fm.tris.iter().map(move |t| {
                    tri_area(point(m, fm, t[0]), point(m, fm, t[1]), point(m, fm, t[2]))
                })
            })
            .sum()
    }

    fn assert_area(name: &str, rel: f64) {
        let (got, want) = (total_area(&mesh(name)), fixture_area(name));
        assert!(
            (got - want).abs() <= rel * want,
            "{name}: area {got} vs {want}"
        );
    }

    #[test]
    fn planar_face_with_hole_area_exact() {
        // The planar faces are triangulated from their boundary polylines alone, so their area is
        // that of the polygons; with the cylindrical hole wall the total is within deflection.
        assert_area("box_with_hole.step", 1e-4);
        assert_area("box.step", 1e-12);
    }

    #[test]
    fn cylinder_lateral_area_within_deflection() {
        assert_area("cylinder.step", 1e-3);
    }

    #[test]
    fn cylinder_with_hole_across_seam() {
        assert_area("cylinder_hole_seam.step", 1e-3);
        assert_oriented("cylinder_hole_seam.step");
    }

    #[test]
    fn full_sphere() {
        assert_area("sphere.step", 1e-3);
        assert_oriented("sphere.step");
    }

    #[test]
    fn cone_to_apex() {
        assert_area("cone.step", 1e-3);
        assert_oriented("cone.step");
    }

    #[test]
    fn torus_face() {
        assert_area("torus.step", 1e-3);
        assert_oriented("torus.step");
    }

    #[test]
    fn loft_bspline_face_area() {
        assert_area("loft.step", 1e-3);
        assert_area("revolve_spline.step", 1e-3);
        assert_area("extrude_spline.step", 1e-3);
        assert_area("fillet_box.step", 1e-3);
    }

    /// Every triangle's normal points the way of its face (surface normal, flipped when the face
    /// does not follow its surface).
    fn assert_oriented(name: &str) {
        for m in mesh(name) {
            for (si, fi, fm) in &m.faces {
                let f = &m.brep.shells[*si].faces[*fi];
                for t in &fm.tris {
                    let [a, b, c] = t.map(|r| point(&m, fm, r));
                    let (u, v) = (
                        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                    );
                    let n = [
                        u[1] * v[2] - u[2] * v[1],
                        u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0],
                    ];
                    if dist(n, [0.0; 3]) <= 1e-9 * dist(u, [0.0; 3]) * dist(v, [0.0; 3]) {
                        continue;
                    }
                    let centroid = [
                        (a[0] + b[0] + c[0]) / 3.0,
                        (a[1] + b[1] + c[1]) / 3.0,
                        (a[2] + b[2] + c[2]) / 3.0,
                    ];
                    let (pu, pv) = f.surface.project(centroid, None).unwrap();
                    let s = f.surface.normal(pu, pv);
                    let d = (n[0] * s[0] + n[1] * s[1] + n[2] * s[2])
                        * if f.same_sense { 1.0 } else { -1.0 };
                    // Slivers standing on the surface next to a pole have a normal at right angles to it;
                    // a flipped triangle points the other way.
                    assert!(
                        d > -0.01 * dist(n, [0.0; 3]),
                        "{name}: face #{} has a triangle against its normal",
                        f.id
                    );
                }
            }
        }
    }

    #[test]
    fn every_triangle_within_deflection() {
        let dir = format!("{}/tests/fixtures", env!("CARGO_MANIFEST_DIR"));
        for f in std::fs::read_dir(dir).unwrap() {
            let path = f.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_string();
            if !name.ends_with(".step") || name == "external_reference.step" {
                continue;
            }
            for m in mesh(&name) {
                for (si, fi, fm) in &m.faces {
                    let face = &m.brep.shells[*si].faces[*fi];
                    if matches!(face.surface, Surface::Plane(_)) {
                        continue;
                    }
                    for t in &fm.tris {
                        let [a, b, c] = t.map(|r| point(&m, fm, r));
                        let centroid = [
                            (a[0] + b[0] + c[0]) / 3.0,
                            (a[1] + b[1] + c[1]) / 3.0,
                            (a[2] + b[2] + c[2]) / 3.0,
                        ];
                        let (u, v) = face.surface.project(centroid, None).unwrap();
                        let d = dist(face.surface.point(u, v), centroid);
                        assert!(
                            d <= DEFLECTION * 1.5 + 1e-9,
                            "{name} face #{}: centroid {d} from the surface",
                            face.id
                        );
                    }
                }
            }
        }
    }
}
