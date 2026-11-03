use std::collections::HashMap;

use crate::model::{P3, cross, dot, sub};

/// Every directed edge is used once and its reverse once: a closed, consistently oriented
/// 2-manifold surface (possibly touching itself at vertices).
pub fn closed_oriented(tris: &[[u32; 3]]) -> bool {
    let mut directed: HashMap<(u32, u32), u32> = HashMap::with_capacity(tris.len() * 3);
    for t in tris {
        for k in 0..3 {
            let n = directed.entry((t[k], t[(k + 1) % 3])).or_default();
            *n += 1;
            if *n > 1 {
                return false;
            }
        }
    }
    !tris.is_empty()
        && directed
            .keys()
            .all(|&(a, b)| directed.contains_key(&(b, a)))
}

/// Volume enclosed by counter-clockwise (outward facing) triangles.
pub fn signed_volume(positions: &[P3], tris: &[[u32; 3]]) -> f64 {
    let Some(first) = tris.first() else {
        return 0.0;
    };
    let origin = positions[first[0] as usize];
    tris.iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| sub(positions[i as usize], origin));
            dot(a, cross(b, c))
        })
        .sum::<f64>()
        / 6.0
}
