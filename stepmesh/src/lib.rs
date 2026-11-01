//! Reads STEP (ISO 10303-21) files and turns their B-rep solids into watertight triangle meshes.
#![forbid(unsafe_code)]

mod error;
mod fxhash;
mod geom;
mod mesh;
mod model;
pub mod p21;
mod read;

pub use error::Error;
pub use read::{SolidMesh, StepMesh, read_step};

/// Tessellation tolerances.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Largest distance between the mesh and the exact surface, in millimetres.
    pub linear_deflection: f64,
    /// Largest angle between neighbouring mesh edges along a curve, in radians.
    pub angular_deflection: f64,
    /// Fail with [`Error::TriangleLimit`] above this many triangles.
    pub max_triangles: Option<usize>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            linear_deflection: 0.003,
            angular_deflection: 0.5,
            max_triangles: None,
        }
    }
}

/// Intermediate stages, exposed for the crate's own tests; not a stable interface.
#[doc(hidden)]
pub mod internals {
    pub use crate::geom::{Curve, Nurbs1, Nurbs2, Surface, expand_knots};
    pub use crate::mesh::{FaceMesh, VRef, discretize_edge, tessellate_face};
    pub use crate::model::{
        Brep, Edge, Face, Frame, Instance, Loop, Model, P3, Part, Shell, Transform, Units,
        decode_brep, decode_curve, decode_model, decode_surface,
    };
}
