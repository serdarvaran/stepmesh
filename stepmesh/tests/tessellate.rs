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
                    assert!(dist(w[0], w[1]) > 1e-9, "{path:?} #{eid} repeats a point");
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
                        let (u, v) = face.surface.project(centroid, None).unwrap_or_else(|| {
                            panic!("{name} face #{}: no projection for {centroid:?}", face.id)
                        });
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

mod seamless {
    use super::*;
    use std::collections::HashMap;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    use stepmesh::internals::{Face, FaceMesh, Loop, Surface, VRef, tessellate_face};

    fn frame() -> Frame {
        Frame::new([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
    }

    /// A circle at height `h` on a sphere of radius `r`, split into three arcs that run
    /// counter-clockwise (seen from +z) when `ccw`.
    type Polylines = HashMap<u64, Vec<[f64; 3]>>;

    fn ring(r: f64, h: f64, ccw: bool, first_id: u64) -> (Polylines, Vec<(u64, bool)>) {
        let rho = (r * r - h * h).sqrt();
        let circle = Frame::new(
            [0.0, 0.0, h],
            [0.0, 0.0, if ccw { 1.0 } else { -1.0 }],
            [1.0, 0.0, 0.0],
        )
        .unwrap();
        let at = |a: f64| circle.to_world([rho * a.cos(), rho * a.sin(), 0.0]);
        let mut polys = HashMap::new();
        let mut list = Vec::new();
        for k in 0..3 {
            let (a, b) = (TAU * k as f64 / 3.0, TAU * (k + 1) as f64 / 3.0);
            let e = Edge {
                curve: Curve::Circle {
                    frame: circle.clone(),
                    r: rho,
                },
                start: at(a),
                end: at(b),
                start_vertex: 100 + k,
                end_vertex: 100 + (k + 1) % 3,
                same_sense: true,
            };
            polys.insert(
                first_id + k,
                discretize_edge(&e, &Options::default()).unwrap(),
            );
            list.push((first_id + k, true));
        }
        (polys, list)
    }

    fn area(fm: &FaceMesh, polys: &HashMap<u64, Vec<[f64; 3]>>) -> f64 {
        let pt = |r: VRef| match r {
            VRef::Edge { edge, i } => polys[&edge][i],
            VRef::Inner(i) => fm.inner[i],
        };
        fm.tris
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(pt);
                let (u, v) = (
                    [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                    [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                );
                0.5 * dist(
                    [
                        u[1] * v[2] - u[2] * v[1],
                        u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0],
                    ],
                    [0.0; 3],
                )
            })
            .sum()
    }

    #[test]
    fn spherical_cap_without_a_seam() {
        let (r, h) = (10.0, 6.0);
        // North cap: the ring runs counter-clockwise seen from outside (+z), the cap on its left.
        let (polys, list) = ring(r, h, true, 1);
        let face = Face {
            id: 1,
            surface: Surface::Sphere { frame: frame(), r },
            same_sense: true,
            outer: Some(0),
            loops: vec![(Loop::Edges(list.clone()), true)],
        };
        let fm = tessellate_face(&face, &polys, &Options::default()).unwrap();
        let cap = TAU * r * (r - h);
        assert!(
            (area(&fm, &polys) - cap).abs() < 1e-3 * cap,
            "north cap {} vs {cap}",
            area(&fm, &polys)
        );
        // The same ring bounding the rest of the sphere: the loop reversed.
        let rest_face = Face {
            loops: vec![(Loop::Edges(list), false)],
            ..face
        };
        let fm = tessellate_face(&rest_face, &polys, &Options::default()).unwrap();
        let rest = 4.0 * PI * r * r - cap;
        assert!(
            (area(&fm, &polys) - rest).abs() < 1e-3 * rest,
            "rest {} vs {rest}",
            area(&fm, &polys)
        );
        let _ = FRAC_PI_2;
    }

    #[test]
    fn cylinder_band_without_a_seam() {
        let (r, h) = (5.0, 8.0);
        let cyl = Surface::Cylinder { frame: frame(), r };
        let mut polys = HashMap::new();
        let mut loops = Vec::new();
        for (z, ccw, id) in [(0.0, true, 1u64), (h, false, 2)] {
            let circle = Frame::new(
                [0.0, 0.0, z],
                [0.0, 0.0, if ccw { 1.0 } else { -1.0 }],
                [1.0, 0.0, 0.0],
            )
            .unwrap();
            let start = circle.to_world([r, 0.0, 0.0]);
            let e = Edge {
                curve: Curve::Circle { frame: circle, r },
                start,
                end: start,
                start_vertex: id,
                end_vertex: id,
                same_sense: true,
            };
            polys.insert(id, discretize_edge(&e, &Options::default()).unwrap());
            loops.push((Loop::Edges(vec![(id, true)]), true));
        }
        let face = Face {
            id: 1,
            surface: cyl,
            same_sense: true,
            outer: None,
            loops,
        };
        let fm = tessellate_face(&face, &polys, &Options::default()).unwrap();
        let band = TAU * r * h;
        assert!(
            (area(&fm, &polys) - band).abs() < 1e-3 * band,
            "band {} vs {band}",
            area(&fm, &polys)
        );
    }

    #[test]
    fn sphere_with_a_hole_off_its_poles() {
        // A circle round the +x axis, 30° from it: the loop does not go round the sphere, so
        // the face is either the cap inside it or the rest of the sphere outside it.
        let (r, half) = (10.0, PI / 6.0);
        let circle =
            Frame::new([r * half.cos(), 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let rho = r * half.sin();
        let start = circle.to_world([rho, 0.0, 0.0]);
        let e = Edge {
            curve: Curve::Circle {
                frame: circle,
                r: rho,
            },
            start,
            end: start,
            start_vertex: 1,
            end_vertex: 1,
            same_sense: true,
        };
        let polys = HashMap::from([(1u64, discretize_edge(&e, &Options::default()).unwrap())]);
        let cap = TAU * r * r * (1.0 - half.cos());
        // Counter-clockwise seen from outside (+x): the cap lies on its left.
        for (ccw, want) in [(true, cap), (false, 4.0 * PI * r * r - cap)] {
            let face = Face {
                id: 1,
                surface: Surface::Sphere { frame: frame(), r },
                same_sense: true,
                outer: Some(0),
                loops: vec![(Loop::Edges(vec![(1, true)]), ccw)],
            };
            let fm = tessellate_face(&face, &polys, &Options::default()).unwrap();
            let got = area(&fm, &polys);
            assert!(
                (got - want).abs() < 1e-3 * want,
                "ccw {ccw}: {got} vs {want}"
            );
        }
    }

    #[test]
    fn torus_band_lies_on_the_side_its_loops_give() {
        // Between the outer equator (v = 0) and the top circle (v = π/2) a torus has two bands:
        // the quarter and the other three quarters. The loops' direction picks one.
        let (major, minor) = (40.0, 10.0);
        let mut polys = HashMap::new();
        let mut circles = Vec::new();
        for (z, rho, id) in [(0.0, major + minor, 1u64), (minor, major, 2)] {
            let start = [rho, 0.0, z];
            let e = Edge {
                curve: Curve::Circle {
                    frame: Frame::new([0.0, 0.0, z], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
                    r: rho,
                },
                start,
                end: start,
                start_vertex: id,
                end_vertex: id,
                same_sense: true,
            };
            polys.insert(id, discretize_edge(&e, &Options::default()).unwrap());
            circles.push(id);
        }
        let quarter = TAU * minor * (major * FRAC_PI_2 + minor);
        let whole = TAU * TAU * major * minor;
        // Counter-clockwise (seen from +z) along the equator, clockwise along the top circle:
        // the face lies on the +v side of the first, between the two.
        for (equator_ccw, want) in [(true, quarter), (false, whole - quarter)] {
            let face = Face {
                id: 1,
                surface: Surface::Torus {
                    frame: frame(),
                    major,
                    minor,
                },
                same_sense: true,
                outer: Some(1),
                loops: vec![
                    (Loop::Edges(vec![(circles[0], true)]), equator_ccw),
                    (Loop::Edges(vec![(circles[1], true)]), !equator_ccw),
                ],
            };
            let fm = tessellate_face(&face, &polys, &Options::default()).unwrap();
            let got = area(&fm, &polys);
            assert!(
                (got - want).abs() < 1e-3 * want,
                "equator ccw {equator_ccw}: {got} vs {want}"
            );
        }
    }

    #[test]
    fn band_seam_follows_a_wavy_tube() {
        // A circle swept along a wavy axis of four quadratic humps (tangent-continuous at their
        // joints): the seam (u = 0) sampled only at 1/4, 1/2 and 3/4 of its length lies on its
        // chord there, with parallel end tangents, though it is far from straight.
        use stepmesh::internals::{Nurbs2, expand_knots};
        let (r, a, h) = (5.0, 10.0, 10.0);
        let w = FRAC_1_SQRT_2;
        let circle = [
            ([1.0, 0.0], 1.0),
            ([1.0, 1.0], w),
            ([0.0, 1.0], 1.0),
            ([-1.0, 1.0], w),
            ([-1.0, 0.0], 1.0),
            ([-1.0, -1.0], w),
            ([0.0, -1.0], 1.0),
            ([1.0, -1.0], w),
            ([1.0, 0.0], 1.0),
        ];
        let offsets = [0.0, a, 0.0, -a, 0.0, a, 0.0, -a, 0.0];
        let ctrl = circle
            .iter()
            .map(|&([x, y], wt)| {
                offsets
                    .iter()
                    .enumerate()
                    .map(|(j, &dx)| {
                        let p = [r * x + dx, r * y, h * j as f64];
                        [p[0] * wt, p[1] * wt, p[2] * wt, wt]
                    })
                    .collect()
            })
            .collect();
        let knots = expand_knots(&[3, 2, 2, 2, 3], &[0.0, 0.25, 0.5, 0.75, 1.0]);
        let tube = Surface::Nurbs(Nurbs2 {
            du: 2,
            dv: 2,
            ku: knots.clone(),
            kv: knots,
            ctrl,
            closed_u: true,
            closed_v: false,
        });
        let mut polys = HashMap::new();
        let mut loops = Vec::new();
        for (z, ccw, id) in [(0.0, true, 1u64), (8.0 * h, false, 2)] {
            let frame = Frame::new(
                [0.0, 0.0, z],
                [0.0, 0.0, if ccw { 1.0 } else { -1.0 }],
                [1.0, 0.0, 0.0],
            )
            .unwrap();
            let start = frame.to_world([r, 0.0, 0.0]);
            let e = Edge {
                curve: Curve::Circle { frame, r },
                start,
                end: start,
                start_vertex: id,
                end_vertex: id,
                same_sense: true,
            };
            polys.insert(id, discretize_edge(&e, &Options::default()).unwrap());
            loops.push((Loop::Edges(vec![(id, true)]), true));
        }
        let face = Face {
            id: 1,
            surface: tube.clone(),
            same_sense: true,
            outer: None,
            loops,
        };
        let fm = tessellate_face(&face, &polys, &Options::default()).unwrap();
        let n = 2000;
        let mut want = 0.0;
        for i in 0..n {
            for j in 0..n {
                let (u, v) = ((i as f64 + 0.5) / n as f64, (j as f64 + 0.5) / n as f64);
                let (_, su, sv) = tube.derivs(u, v);
                let c = [
                    su[1] * sv[2] - su[2] * sv[1],
                    su[2] * sv[0] - su[0] * sv[2],
                    su[0] * sv[1] - su[1] * sv[0],
                ];
                want += dist(c, [0.0; 3]) / (n * n) as f64;
            }
        }
        let got = area(&fm, &polys);
        assert!((got - want).abs() < 1e-3 * want, "tube {got} vs {want}");
    }
}

#[test]
fn quarter_sphere_through_a_pole() {
    use std::collections::HashMap;
    use std::f64::consts::{FRAC_PI_2, PI};
    use stepmesh::internals::{Face, Loop, Surface, VRef, tessellate_face};
    // Sphere r = 10 about z. Face: the lune between the meridians u = 0 and u = pi below the
    // equator: half the equator (u from pi down to 0) then down the meridian u = 0 to the south
    // pole and up the meridian u = pi. Its loop passes through the pole instead of round it.
    let r = 10.0;
    let f = Frame::new([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let s = Surface::Sphere {
        frame: f.clone(),
        r,
    };
    let o = Options::default();
    let south = [0.0, 0.0, -r];
    let equator = Curve::Circle {
        frame: f.clone(),
        r,
    };
    // Meridian circles: in the xz plane (through u = 0 and u = pi).
    let meridian = Curve::Circle {
        frame: Frame::new([0.0; 3], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]).unwrap(),
        r,
    };
    let edges = [
        (
            1,
            Edge {
                curve: equator,
                start: [-r, 0.0, 0.0],
                end: [r, 0.0, 0.0],
                start_vertex: 10,
                end_vertex: 11,
                same_sense: false,
            },
        ),
        (
            2,
            Edge {
                curve: meridian.clone(),
                start: [r, 0.0, 0.0],
                end: south,
                start_vertex: 11,
                end_vertex: 12,
                same_sense: false,
            },
        ),
        (
            3,
            Edge {
                curve: meridian,
                start: south,
                end: [-r, 0.0, 0.0],
                start_vertex: 12,
                end_vertex: 10,
                same_sense: false,
            },
        ),
    ];
    let polys: HashMap<u64, Vec<[f64; 3]>> = edges
        .iter()
        .map(|(id, e)| (*id, discretize_edge(e, &o).unwrap()))
        .collect();
    for (same_sense, list) in [
        (true, vec![(1, true), (2, true), (3, true)]),
        (false, vec![(3, false), (2, false), (1, false)]),
    ] {
        let face = Face {
            id: 1,
            surface: s.clone(),
            same_sense,
            outer: Some(0),
            loops: vec![(Loop::Edges(list), true)],
        };
        let fm = tessellate_face(&face, &polys, &o).unwrap();
        let pt = |q: VRef| match q {
            VRef::Edge { edge, i } => polys[&edge][i],
            VRef::Inner(i) => fm.inner[i],
        };
        let area: f64 = fm
            .tris
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(pt);
                let (u, v) = (
                    [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                    [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                );
                0.5 * dist(
                    [
                        u[1] * v[2] - u[2] * v[1],
                        u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0],
                    ],
                    [0.0; 3],
                )
            })
            .sum();
        let quarter = PI * r * r;
        assert!(
            (area - quarter).abs() < 1e-3 * quarter,
            "same_sense {same_sense}: area {area} vs {quarter}"
        );
        let _ = FRAC_PI_2;
    }
}

#[test]
fn quarter_sphere_through_a_pole_vertex_slightly_off() {
    use std::collections::HashMap;
    use std::f64::consts::{FRAC_PI_2, PI};
    use stepmesh::internals::{Face, Loop, Surface, VRef, tessellate_face};
    // Sphere r = 10 about z. Face: the lune between the meridians u = 0 and u = pi below the
    // equator: half the equator (u from pi down to 0) then down the meridian u = 0 to the south
    // pole and up the meridian u = pi. Its loop passes through the pole instead of round it.
    let r = 10.0;
    let f = Frame::new([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let s = Surface::Sphere {
        frame: f.clone(),
        r,
    };
    let o = Options::default();
    // The pole vertex a little off the pole, as file tolerances leave it.
    for off in [
        [0.0, 1e-7, 0.0],
        [1e-7, 0.0, 0.0],
        [-1e-7, -1e-7, 0.0],
        [3e-6, -2e-6, 0.0],
        [0.0, -1e-6, 0.0],
    ] {
        let south = [off[0], off[1], -r];
        let equator = Curve::Circle {
            frame: f.clone(),
            r,
        };
        // Meridian circles: in the xz plane (through u = 0 and u = pi).
        let meridian = Curve::Circle {
            frame: Frame::new([0.0; 3], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]).unwrap(),
            r,
        };
        let edges = [
            (
                1,
                Edge {
                    curve: equator,
                    start: [-r, 0.0, 0.0],
                    end: [r, 0.0, 0.0],
                    start_vertex: 10,
                    end_vertex: 11,
                    same_sense: false,
                },
            ),
            (
                2,
                Edge {
                    curve: meridian.clone(),
                    start: [r, 0.0, 0.0],
                    end: south,
                    start_vertex: 11,
                    end_vertex: 12,
                    same_sense: false,
                },
            ),
            (
                3,
                Edge {
                    curve: meridian,
                    start: south,
                    end: [-r, 0.0, 0.0],
                    start_vertex: 12,
                    end_vertex: 10,
                    same_sense: false,
                },
            ),
        ];
        let polys: HashMap<u64, Vec<[f64; 3]>> = edges
            .iter()
            .map(|(id, e)| (*id, discretize_edge(e, &o).unwrap()))
            .collect();
        for (same_sense, list) in [
            (true, vec![(1, true), (2, true), (3, true)]),
            (false, vec![(3, false), (2, false), (1, false)]),
        ] {
            let face = Face {
                id: 1,
                surface: s.clone(),
                same_sense,
                outer: Some(0),
                loops: vec![(Loop::Edges(list), true)],
            };
            let fm = tessellate_face(&face, &polys, &o).unwrap();
            let pt = |q: VRef| match q {
                VRef::Edge { edge, i } => polys[&edge][i],
                VRef::Inner(i) => fm.inner[i],
            };
            let area: f64 = fm
                .tris
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(pt);
                    let (u, v) = (
                        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                    );
                    0.5 * dist(
                        [
                            u[1] * v[2] - u[2] * v[1],
                            u[2] * v[0] - u[0] * v[2],
                            u[0] * v[1] - u[1] * v[0],
                        ],
                        [0.0; 3],
                    )
                })
                .sum();
            let quarter = PI * r * r;
            assert!(
                (area - quarter).abs() < 1e-3 * quarter,
                "same_sense {same_sense}: area {area} vs {quarter}"
            );
            let _ = FRAC_PI_2;
        }
    }
}
