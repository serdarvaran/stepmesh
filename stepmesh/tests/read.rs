use std::collections::HashMap;

use stepmesh::{Error, Options, StepMesh, read_step};

fn fixture(name: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/tests/fixtures/{name}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

/// (volume, solid count) per fixture from fixtures.json; `None` for an expected error.
fn expected() -> Vec<(String, Option<(f64, usize)>)> {
    let text = fixture("fixtures.json");
    let mut out = Vec::new();
    for block in text.split("\n  \"").skip(1) {
        let name = block[..block.find('"').unwrap()].to_string();
        let field = |key: &str| {
            block.find(&format!("\"{key}\": ")).map(|i| {
                let rest = &block[i + key.len() + 4..];
                rest[..rest.find([',', '\n']).unwrap()].trim().to_string()
            })
        };
        let entry = match (field("volume"), field("solids")) {
            (Some(v), Some(s)) => Some((v.parse().unwrap(), s.parse().unwrap())),
            _ => None,
        };
        out.push((name, entry));
    }
    out
}

fn read(text: &str) -> Result<StepMesh, Error> {
    read_step(text.as_bytes(), &Options::default())
}

fn volume(m: &StepMesh) -> f64 {
    m.solids.iter().map(|s| s.volume()).sum()
}

#[test]
fn fixtures_match_analytic_volume() {
    let all = expected();
    assert!(all.len() >= 16);
    for (name, want) in all {
        let got = read(&fixture(&name));
        match want {
            Some((v, solids)) => {
                let m = got.unwrap_or_else(|e| panic!("{name}: {e}"));
                assert_eq!(m.solids.len(), solids, "{name}");
                assert!(
                    (volume(&m) - v).abs() <= 1e-3 * v,
                    "{name}: volume {} vs {v}",
                    volume(&m)
                );
            }
            None => assert_eq!(got.err(), Some(Error::NoGeometry), "{name}"),
        }
    }
}

#[test]
fn fixtures_are_closed_oriented_manifolds() {
    for (name, want) in expected() {
        if want.is_none() {
            continue;
        }
        for (i, s) in read(&fixture(&name)).unwrap().solids.iter().enumerate() {
            let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
            for t in &s.triangles {
                for k in 0..3 {
                    *directed.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
                }
            }
            for (&(a, b), &n) in &directed {
                assert_eq!(n, 1, "{name} solid {i}: edge {a}->{b} used {n} times");
                assert_eq!(
                    directed.get(&(b, a)),
                    Some(&1),
                    "{name} solid {i}: edge {a}->{b} has no twin"
                );
            }
        }
    }
}

#[test]
fn mixed_unit_assembly_volumes() {
    let m = read(&fixture("mixed_units.step")).unwrap();
    assert!(
        (volume(&m) - 16_387.064).abs() <= 1e-3 * 16_387.064,
        "{}",
        volume(&m)
    );
    let [lo, hi] = m.solids[0].bounds();
    assert!(((lo[0] + hi[0]) / 2.0 - 30.0).abs() < 1e-6, "{lo:?} {hi:?}");
}

/// box.step with edge #21 (a line from #22 up to #24) written the other way round: a line that
/// runs downwards from #24's point, used with same_sense .F.
fn box_with_reversed_edge() -> String {
    let text = fixture("box.step");
    let edge = "#21 = EDGE_CURVE('',#22,#24,#26,.T.);";
    assert!(text.contains(edge));
    text.replacen(edge, "#21 = EDGE_CURVE('',#22,#24,#9003,.F.);", 1)
        .replacen(
            "ENDSEC;\nEND-ISO",
            "#9001 = DIRECTION('',(0.,0.,-1.));\n#9002 = VECTOR('',#9001,1.);\n\
         #9003 = LINE('',#25,#9002);\nENDSEC;\nEND-ISO",
            1,
        )
}

#[test]
fn reversed_edge_curve_is_shared() {
    let m = read(&box_with_reversed_edge()).unwrap();
    assert!((volume(&m) - 6000.0).abs() < 1e-9, "{}", volume(&m));
}

#[test]
fn zero_length_edge_joins_its_vertices() {
    // Face #17 reaches the corner of vertex #22 through a second vertex at the same point; a
    // sliver face holds the zero-length edge between the two, as some exporters write it.
    let text = fixture("box.step");
    let face_loop = "#19 = EDGE_LOOP('',(#20,#55,#83,#111));";
    let shell = "#16 = CLOSED_SHELL('',(#17,";
    assert!(text.contains(face_loop) && text.contains(shell));
    let text = text
        .replacen(face_loop, "#19 = EDGE_LOOP('',(#9020,#9022,#83,#111));", 1)
        .replacen(shell, "#16 = CLOSED_SHELL('',(#9030,#17,", 1)
        .replacen(
            "ENDSEC;\nEND-ISO",
            "#9010 = VERTEX_POINT('',#9014);\n\
             #9014 = CARTESIAN_POINT('',(-5.,-10.,-15.));\n\
             #9011 = EDGE_CURVE('',#9010,#24,#26,.T.);\n\
             #9012 = EDGE_CURVE('',#9010,#22,#26,.T.);\n\
             #9013 = EDGE_CURVE('',#9010,#57,#59,.T.);\n\
             #9020 = ORIENTED_EDGE('',*,*,#9011,.F.);\n\
             #9021 = ORIENTED_EDGE('',*,*,#9012,.T.);\n\
             #9022 = ORIENTED_EDGE('',*,*,#9013,.T.);\n\
             #9023 = ORIENTED_EDGE('',*,*,#9012,.F.);\n\
             #9030 = ADVANCED_FACE('',(#9031),#32,.F.);\n\
             #9031 = FACE_OUTER_BOUND('',#9032,.T.);\n\
             #9032 = EDGE_LOOP('',(#9021,#9023));\nENDSEC;\nEND-ISO",
            1,
        );
    let m = read(&text).unwrap();
    assert!((volume(&m) - 6000.0).abs() < 1e-9, "{}", volume(&m));
}

#[test]
fn missing_face_is_not_watertight() {
    let text = fixture("box.step");
    let shell = text.lines().find(|l| l.contains("CLOSED_SHELL(")).unwrap();
    let cut = shell.replacen(",#137", "", 1);
    assert_ne!(cut, shell);
    assert_eq!(
        read(&text.replacen(shell, &cut, 1)).err(),
        Some(Error::NotWatertight { solid: 0 })
    );
}

#[test]
fn triangle_limit() {
    let o = Options {
        max_triangles: Some(10),
        ..Options::default()
    };
    let r = read_step(fixture("box.step").as_bytes(), &o);
    assert!(
        matches!(r, Err(Error::TriangleLimit { limit: 10, .. })),
        "{r:?}"
    );
}

#[test]
fn instances_are_placed() {
    let m = read(&fixture("assembly_two_boxes.step")).unwrap();
    let mut x: Vec<f64> = m
        .solids
        .iter()
        .map(|s| {
            let [lo, hi] = s.bounds();
            (lo[0] + hi[0]) / 2.0
        })
        .collect();
    x.sort_by(f64::total_cmp);
    assert!(
        (x[0] - 0.0).abs() < 1e-9 && (x[1] - 50.0).abs() < 1e-9,
        "{x:?}"
    );
}
