//! Tessellation of B-rep solids.

mod check;
mod edges;
mod face;
mod param;
mod seam;
mod solid;

pub(crate) use check::signed_volume;
pub use edges::discretize_edge;
pub use face::{FaceMesh, VRef, tessellate_face};
pub use solid::edge_polylines;
pub(crate) use solid::{LocalMesh, mesh_brep, triangle_area};
