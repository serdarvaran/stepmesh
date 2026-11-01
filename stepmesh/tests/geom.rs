use std::f64::consts::{FRAC_PI_2, PI, TAU};

use stepmesh::internals::{
    Curve, Frame, Nurbs1, Nurbs2, Surface, Units, decode_curve, decode_surface, expand_knots,
};
use stepmesh::p21::parse;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn frame() -> Frame {
    Frame::new([1.0, -2.0, 3.0], [0.3, 0.4, 1.0], [1.0, 0.2, 0.0]).unwrap()
}

/// Rational quadratic full circle of radius `r` in the XY plane (9 control points).
fn nurbs_circle(r: f64) -> Nurbs1 {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let pts = [
        (1., 0., 1.),
        (1., 1., w),
        (0., 1., 1.),
        (-1., 1., w),
        (-1., 0., 1.),
        (-1., -1., w),
        (0., -1., 1.),
        (1., -1., w),
        (1., 0., 1.),
    ];
    Nurbs1 {
        degree: 2,
        knots: expand_knots(&[3, 2, 2, 2, 3], &[0.0, 1.0, 2.0, 3.0, 4.0]),
        ctrl: pts
            .iter()
            .map(|&(x, y, w)| [x * r * w, y * r * w, 0.0, w])
            .collect(),
    }
}

fn nurbs_patch(rng: &mut Rng) -> Nurbs2 {
    let ctrl = (0..4)
        .map(|i| {
            (0..5)
                .map(|j| {
                    let w = rng.range(0.6, 1.4);
                    [
                        i as f64 * 3.0 * w,
                        j as f64 * 2.0 * w,
                        rng.range(-1.0, 1.0) * w,
                        w,
                    ]
                })
                .collect()
        })
        .collect();
    Nurbs2 {
        du: 3,
        dv: 2,
        ku: expand_knots(&[4, 4], &[0.0, 1.0]),
        kv: expand_knots(&[3, 1, 1, 3], &[0.0, 0.3, 0.7, 1.0]),
        ctrl,
        closed_u: false,
        closed_v: false,
    }
}

fn surfaces(rng: &mut Rng) -> Vec<(&'static str, Surface)> {
    let f = frame();
    vec![
        ("plane", Surface::Plane(f.clone())),
        (
            "cylinder",
            Surface::Cylinder {
                frame: f.clone(),
                r: 4.0,
            },
        ),
        (
            "cone",
            Surface::Cone {
                frame: f.clone(),
                r: 3.0,
                half_angle: 0.4,
            },
        ),
        (
            "sphere",
            Surface::Sphere {
                frame: f.clone(),
                r: 5.0,
            },
        ),
        (
            "torus",
            Surface::Torus {
                frame: f.clone(),
                major: 10.0,
                minor: 3.0,
            },
        ),
        ("nurbs", Surface::Nurbs(nurbs_patch(rng))),
        (
            "extrusion",
            Surface::Extrusion {
                curve: Curve::Nurbs(nurbs_circle(2.0)),
                dir: [0.0, 0.6, 0.8],
            },
        ),
        (
            "revolution",
            Surface::Revolution {
                curve: Curve::Line {
                    origin: [3.0, 0.0, 0.0],
                    dir: [0.287_347_885_566_345_2, 0.0, 0.957_826_285_221_150_7],
                },
                axis_origin: [0.0, 0.0, 0.0],
                axis_dir: [0.0, 0.0, 1.0],
            },
        ),
    ]
}

#[test]
fn rational_quadratic_circle_is_exact() {
    let c = Curve::Nurbs(nurbs_circle(3.0));
    for i in 0..=1000 {
        let p = c.point(4.0 * i as f64 / 1000.0);
        assert!(
            ((p[0].powi(2) + p[1].powi(2)).sqrt() - 3.0).abs() < 1e-12,
            "{p:?}"
        );
    }
}

#[test]
fn knot_multiplicities_expand() {
    assert_eq!(
        expand_knots(&[4, 4], &[0.0, 1.0]),
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    );
}

#[test]
fn curve_derivative_matches_finite_difference() {
    let curves = [
        Curve::Nurbs(nurbs_circle(3.0)),
        Curve::Circle {
            frame: frame(),
            r: 2.0,
        },
        Curve::Ellipse {
            frame: frame(),
            a: 4.0,
            b: 1.5,
        },
    ];
    for c in &curves {
        for t in [0.3, 1.1, 2.5] {
            let h = 1e-6;
            let (a, b) = (c.point(t - h), c.point(t + h));
            let fd = [
                (b[0] - a[0]) / (2.0 * h),
                (b[1] - a[1]) / (2.0 * h),
                (b[2] - a[2]) / (2.0 * h),
            ];
            let d = c.derivative(t);
            assert!(
                dist(fd, d) < 1e-5 * (1.0 + dist(d, [0.0; 3])),
                "{c:?} {fd:?} {d:?}"
            );
        }
    }
}

#[test]
fn curve_project_round_trip() {
    let curves = [
        Curve::Line {
            origin: [1.0, 2.0, 3.0],
            dir: [0.6, 0.8, 0.0],
        },
        Curve::Nurbs(nurbs_circle(3.0)),
        Curve::Circle {
            frame: frame(),
            r: 2.0,
        },
        Curve::Ellipse {
            frame: frame(),
            a: 4.0,
            b: 1.5,
        },
    ];
    let mut rng = Rng(7);
    for c in &curves {
        let (t0, t1) = c.param_range();
        for _ in 0..200 {
            let t = rng.range(t0.max(-50.0), t1.min(50.0));
            let p = c.point(t);
            let q = c.point(c.project(p));
            assert!(dist(p, q) < 1e-9, "{c:?} t={t} {p:?} {q:?}");
        }
    }
}

#[test]
fn surface_project_round_trip() {
    let mut rng = Rng(11);
    for (name, s) in surfaces(&mut rng) {
        let [(u0, u1), (v0, v1)] = s.param_bounds();
        for _ in 0..200 {
            let u = rng.range(u0.max(-20.0), u1.min(20.0));
            let v = rng.range(v0.max(-20.0), v1.min(20.0));
            let p = s.point(u, v);
            let (pu, pv) = s
                .project(p, None)
                .unwrap_or_else(|| panic!("{name}: no projection for ({u}, {v})"));
            let q = s.point(pu, pv);
            assert!(
                dist(p, q) < 1e-9,
                "{name}: ({u}, {v}) -> ({pu}, {pv}), {}",
                dist(p, q)
            );
        }
    }
}

#[test]
fn project_off_surface_finds_closest() {
    let s = Surface::Cylinder {
        frame: frame(),
        r: 4.0,
    };
    let f = frame();
    let on = s.point(1.0, 2.0);
    let n = s.normal(1.0, 2.0);
    let off = [on[0] + n[0], on[1] + n[1], on[2] + n[2]];
    let (u, v) = s.project(off, None).unwrap();
    assert!(dist(s.point(u, v), on) < 1e-9);
    assert!((dist(off, s.point(u, v)) - 1.0).abs() < 1e-9);
    let _ = f;
}

#[test]
fn surface_normal_is_unit_and_orthogonal() {
    let mut rng = Rng(5);
    for (name, s) in surfaces(&mut rng) {
        let [(u0, u1), (v0, v1)] = s.param_bounds();
        for _ in 0..50 {
            let u = rng.range(u0.max(-20.0), u1.min(20.0));
            let v = rng
                .range(v0.max(-20.0), v1.min(20.0))
                .clamp(-FRAC_PI_2 + 0.01, 20.0);
            let v = if name == "sphere" {
                v.min(FRAC_PI_2 - 0.01)
            } else {
                v
            };
            let n = s.normal(u, v);
            assert!((dist(n, [0.0; 3]) - 1.0).abs() < 1e-9, "{name}");
            let h = 1e-6;
            for (a, b) in [
                (s.point(u - h, v), s.point(u + h, v)),
                (s.point(u, v - h), s.point(u, v + h)),
            ] {
                let t = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let len = dist(t, [0.0; 3]);
                if len > 1e-9 {
                    let c = (t[0] * n[0] + t[1] * n[1] + t[2] * n[2]) / len;
                    assert!(c.abs() < 1e-5, "{name}: normal not orthogonal ({c})");
                }
            }
        }
    }
}

#[test]
fn periods() {
    let f = frame();
    assert_eq!(
        Surface::Cylinder {
            frame: f.clone(),
            r: 1.0
        }
        .period(),
        (Some(TAU), None)
    );
    assert_eq!(
        Surface::Torus {
            frame: f.clone(),
            major: 3.0,
            minor: 1.0
        }
        .period(),
        (Some(TAU), Some(TAU))
    );
    assert_eq!(Surface::Plane(f).period(), (None, None));
    let _ = PI;
}

const RATIONAL: &str = "ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
#1=CARTESIAN_POINT('',(3.,0.,0.));#2=CARTESIAN_POINT('',(3.,3.,0.));#3=CARTESIAN_POINT('',(0.,3.,0.));
#4=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#1,#2,#3),.UNSPECIFIED.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)
CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,0.7071067811865476,1.))REPRESENTATION_ITEM(''));
#5=CARTESIAN_POINT('',(0.,0.,0.));#6=DIRECTION('',(0.,0.,1.));#7=DIRECTION('',(1.,0.,0.));
#8=AXIS2_PLACEMENT_3D('',#5,#6,#7);#9=CONICAL_SURFACE('',#8,2.,30.);
#10=OFFSET_SURFACE('',#9,1.,.F.);
ENDSEC;END-ISO-10303-21;";

#[test]
fn rational_bspline_complex_entity_decodes() {
    let x = parse(RATIONAL.as_bytes()).unwrap();
    let c = decode_curve(
        &x,
        4,
        Units {
            length_mm: 2.0,
            angle_rad: 1.0,
        },
    )
    .unwrap();
    let (t0, t1) = c.param_range();
    for i in 0..=20 {
        let p = c.point(t0 + (t1 - t0) * i as f64 / 20.0);
        assert!(
            ((p[0].powi(2) + p[1].powi(2)).sqrt() - 6.0).abs() < 1e-12,
            "quarter circle of radius 3 in units of 2 mm: {p:?}"
        );
    }
}

#[test]
fn cone_angle_uses_the_angle_unit() {
    let x = parse(RATIONAL.as_bytes()).unwrap();
    let s = decode_surface(
        &x,
        9,
        Units {
            length_mm: 1.0,
            angle_rad: PI / 180.0,
        },
    )
    .unwrap();
    match s {
        Surface::Cone { r, half_angle, .. } => {
            assert_eq!(r, 2.0);
            assert!((half_angle - PI / 6.0).abs() < 1e-15);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn unknown_surface_is_unsupported() {
    let x = parse(RATIONAL.as_bytes()).unwrap();
    let r = decode_surface(
        &x,
        10,
        Units {
            length_mm: 1.0,
            angle_rad: 1.0,
        },
    );
    assert_eq!(
        r.err(),
        Some(stepmesh::Error::Unsupported {
            entity: "OFFSET_SURFACE".into(),
            id: 10
        })
    );
}
