//! Tessellation of B-rep solids.

mod edges;
mod face;

pub use edges::discretize_edge;
pub use face::{FaceMesh, VRef, tessellate_face};
