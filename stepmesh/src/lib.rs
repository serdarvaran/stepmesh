//! Reads STEP (ISO 10303-21) files and turns their B-rep solids into watertight triangle meshes.
#![forbid(unsafe_code)]

mod error;
mod geom;
mod model;
pub mod p21;

pub use error::Error;

/// Intermediate stages, exposed for the crate's own tests; not a stable interface.
#[doc(hidden)]
pub mod internals {
    pub use crate::geom::{Curve, Nurbs1, Nurbs2, Surface, expand_knots};
    pub use crate::model::{
        Frame, Instance, Model, P3, Part, Transform, Units, decode_curve, decode_model,
        decode_surface,
    };
}
