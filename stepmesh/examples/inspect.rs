//! Prints each solid's volume, area and bounds: `cargo run --release --example inspect -- file.step`.
fn main() {
    let path = std::env::args().nth(1).expect("path");
    let bytes = std::fs::read(&path).expect("read");
    match stepmesh::read_step(&bytes, &stepmesh::Options::default()) {
        Ok(m) => {
            for (i, s) in m.solids.iter().enumerate() {
                let [lo, hi] = s.bounds();
                println!(
                    "solid {i} {:?}: volume {:.4} area {:.4} triangles {} bounds {lo:.3?} {hi:.3?}",
                    s.name,
                    s.volume(),
                    s.area(),
                    s.triangles.len()
                );
            }
        }
        Err(e) => println!("error: {e}"),
    }
}
