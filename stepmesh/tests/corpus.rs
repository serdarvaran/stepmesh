//! Corpus gate against OpenCascade reference measurements (reference/occt_reference.py).
//!
//!     STEPMESH_CORPUS=list.txt STEPMESH_REFERENCE=ref.jsonl \
//!         cargo test --release --test corpus -- --ignored --nocapture
//!
//! Optional: STEPMESH_TIMEOUT_S (default 120), STEPMESH_OUT (one JSON line per file),
//! STEPMESH_DISPUTED (lines `path<TAB>reason`: files whose OpenCascade reference was shown to be
//! wrong; reported as "reference disputed", never as ok).
//! A file passes when stepmesh returns an error, or returns solids whose count matches, whose
//! bounds are within 0.1 % of the diagonal (or the deflections of both meshes on small parts),
//! and whose volume is within 0.1 % of OpenCascade's own mesh at the same tolerances or within
//! max(0.1 %, deflection x area) of its exact B-rep volume (a mesh within the deflection of the
//! surface cannot be closer than that). Anything else read
//! without an error is WRONG; the gate allows none, and no panic or timeout either. The share of
//! files within 0.1 % of the exact volume is reported separately.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::sync::mpsc;
use std::time::{Duration, Instant};

use stepmesh::{Error, Options, StepMesh, read_step};

/// Floor of the deflection of the fine OpenCascade mesh the reference takes its bounds from
/// (reference/occt_reference.py), in millimetres.
const REFERENCE_BOUNDS_DEFLECTION: f64 = 0.01;

struct Reference {
    solids: usize,
    volume: f64,
    mesh_volume: Option<f64>,
    area: f64,
    min: [f64; 3],
    max: [f64; 3],
}

/// Minimal reader for the reference lines: numbers, arrays of numbers, strings and booleans.
fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let i = line.find(&format!("\"{key}\": "))? + key.len() + 4;
    let rest = &line[i..];
    let end = if rest.starts_with('[') {
        rest.find(']')? + 1
    } else if let Some(stripped) = rest.strip_prefix('"') {
        stripped.find('"')? + 2
    } else {
        rest.find([',', '}'])?
    };
    Some(&rest[..end])
}

fn triple(s: &str) -> Option<[f64; 3]> {
    let v: Vec<f64> = s
        .trim_matches(['[', ']'])
        .split(',')
        .filter_map(|x| x.trim().parse().ok())
        .collect();
    (v.len() == 3).then(|| [v[0], v[1], v[2]])
}

fn references(path: &str) -> HashMap<String, Option<Reference>> {
    let mut out = HashMap::new();
    for line in std::fs::read_to_string(path)
        .expect("reference file")
        .lines()
    {
        let Some(file) = field(line, "file") else {
            continue;
        };
        let file = file.trim_matches('"').to_string();
        let ok = field(line, "ok") == Some("true") && field(line, "valid") == Some("true");
        let r = ok.then(|| Reference {
            solids: field(line, "solids")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0),
            volume: field(line, "volume")
                .and_then(|s| s.parse().ok())
                .unwrap_or(f64::NAN),
            mesh_volume: field(line, "mesh_volume").and_then(|s| s.parse().ok()),
            area: field(line, "area")
                .and_then(|s| s.parse().ok())
                .unwrap_or(f64::NAN),
            min: field(line, "min").and_then(triple).unwrap_or([f64::NAN; 3]),
            max: field(line, "max").and_then(triple).unwrap_or([f64::NAN; 3]),
        });
        out.insert(file, r);
    }
    out
}

fn judge(m: &StepMesh, r: &Reference) -> Result<(), String> {
    let volume: f64 = m.solids.iter().map(|s| s.volume()).sum();
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for s in &m.solids {
        let [a, b] = s.bounds();
        for i in 0..3 {
            lo[i] = lo[i].min(a[i]);
            hi[i] = hi[i].max(b[i]);
        }
    }
    let diag = (0..3)
        .map(|i| (r.max[i] - r.min[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    let bounds = (0..3)
        .map(|i| (lo[i] - r.min[i]).abs().max((hi[i] - r.max[i]).abs()))
        .fold(0.0, f64::max);
    let dv = (volume - r.volume).abs() / r.volume.abs().max(1e-12);
    if m.solids.len() != r.solids {
        return Err(format!("{} solids, reference {}", m.solids.len(), r.solids));
    }
    let allowed =
        1e-3_f64.max(Options::default().linear_deflection * r.area / r.volume.abs().max(1e-12));
    let dm = r
        .mesh_volume
        .map_or(f64::INFINITY, |m| (volume - m).abs() / m.abs().max(1e-12));
    if dv > allowed && dm > 1e-3 {
        return Err(format!(
            "volume {volume} vs {} ({dv:.2e}), mesh {:?} ({dm:.2e})",
            r.volume, r.mesh_volume
        ));
    }
    // The reference bounds come from a mesh at no finer than REFERENCE_BOUNDS_DEFLECTION: on a
    // small part they cannot settle 0.1 % of the diagonal.
    let bounds_allowed =
        (1e-3 * diag).max(REFERENCE_BOUNDS_DEFLECTION + Options::default().linear_deflection);
    if bounds > bounds_allowed {
        return Err(format!("bounds off by {bounds} (diagonal {diag})"));
    }
    Ok(())
}

fn error_class(e: &Error) -> String {
    match e {
        Error::Unsupported { entity, .. } => format!("Unsupported {entity}"),
        Error::Syntax { .. } => "Syntax".into(),
        Error::Invalid { msg, .. } => format!("Invalid ({msg})"),
        Error::Projection { reason, .. } => format!("Projection: {reason}"),
        other => format!("{other:?}")
            .split([' ', '{', '('])
            .next()
            .unwrap_or("?")
            .to_string(),
    }
}

#[test]
#[ignore = "needs STEPMESH_CORPUS and STEPMESH_REFERENCE"]
fn corpus_gate() {
    let list = std::env::var("STEPMESH_CORPUS").expect("STEPMESH_CORPUS");
    let reference = references(&std::env::var("STEPMESH_REFERENCE").expect("STEPMESH_REFERENCE"));
    let timeout = Duration::from_secs(
        std::env::var("STEPMESH_TIMEOUT_S")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(120),
    );
    let mut out = std::env::var("STEPMESH_OUT")
        .ok()
        .map(|p| std::fs::File::create(p).expect("STEPMESH_OUT"));
    let files: Vec<String> = std::fs::read_to_string(list)
        .expect("list")
        .lines()
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect();
    let disputed: HashMap<String, String> = std::env::var("STEPMESH_DISPUTED")
        .ok()
        .map(|p| std::fs::read_to_string(p).expect("STEPMESH_DISPUTED"))
        .unwrap_or_default()
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .collect();
    std::panic::set_hook(Box::new(|_| {}));

    let mut classes: BTreeMap<String, usize> = BTreeMap::new();
    let (mut wrong, mut panics, mut timeouts, mut ok, mut read_by_reference) = (0, 0, 0, 0, 0);
    let mut volume_errors = Vec::new();
    let (mut ours_s, start_all) = (0.0, Instant::now());
    for path in &files {
        let r = reference.get(path).and_then(Option::as_ref);
        if r.is_some() {
            read_by_reference += 1;
        }
        let bytes = std::fs::read(path).expect("corpus file");
        let (tx, rx) = mpsc::channel();
        let start = Instant::now();
        std::thread::spawn(move || {
            let res = std::panic::catch_unwind(|| read_step(&bytes, &Options::default()));
            let _ = tx.send(res);
        });
        let (class, detail) = match rx.recv_timeout(timeout) {
            Err(_) => {
                timeouts += 1;
                ("TIMEOUT".to_string(), String::new())
            }
            Ok(Err(_)) => {
                panics += 1;
                ("PANIC".to_string(), String::new())
            }
            Ok(Ok(Err(e))) => (format!("error: {}", error_class(&e)), e.to_string()),
            Ok(Ok(Ok(_))) if disputed.contains_key(path) => {
                ("reference disputed".to_string(), disputed[path].clone())
            }
            Ok(Ok(Ok(m))) => match r {
                None => ("ok (no reference)".to_string(), String::new()),
                Some(r) => match judge(&m, r) {
                    Ok(()) => {
                        ok += 1;
                        let v: f64 = m.solids.iter().map(|s| s.volume()).sum();
                        volume_errors.push((v - r.volume).abs() / r.volume.abs());
                        ("ok".to_string(), String::new())
                    }
                    Err(why) => {
                        wrong += 1;
                        ("WRONG".to_string(), why)
                    }
                },
            },
        };
        let secs = start.elapsed().as_secs_f64();
        ours_s += secs;
        *classes.entry(class.clone()).or_default() += 1;
        if let Some(f) = out.as_mut() {
            let esc = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"");
            writeln!(
                f,
                "{{\"file\": \"{}\", \"class\": \"{}\", \"detail\": \"{}\", \"s\": {secs:.3}}}",
                esc(path),
                esc(&class),
                esc(&detail)
            )
            .ok();
        }
        if class == "WRONG" || class == "PANIC" || class == "TIMEOUT" {
            eprintln!("{class}: {path} {detail}");
        }
    }
    let elapsed = start_all.elapsed().as_secs_f64();
    volume_errors.sort_by(f64::total_cmp);
    let q = |f: f64| {
        volume_errors
            .get(((volume_errors.len() as f64 - 1.0) * f).round() as usize)
            .copied()
            .unwrap_or(f64::NAN)
    };
    println!(
        "files {}  read by OpenCascade {}  ok {}  ({:.1}% of those)",
        files.len(),
        read_by_reference,
        ok,
        100.0 * ok as f64 / read_by_reference.max(1) as f64
    );
    println!(
        "volume error of ok files: median {:.2e}  p90 {:.2e}  max {:.2e}  within 0.1 %: {}",
        q(0.5),
        q(0.9),
        q(1.0),
        volume_errors.iter().filter(|&&e| e <= 1e-3).count()
    );
    println!("time {ours_s:.1} s ({elapsed:.1} s wall)");
    for (c, n) in &classes {
        println!("{n:6}  {c}");
    }
    assert_eq!(
        (wrong, panics, timeouts),
        (0, 0, 0),
        "WRONG, PANIC, TIMEOUT"
    );
}
