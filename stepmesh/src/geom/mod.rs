//! Curve and surface evaluation.

mod curve;
mod nurbs;
mod surface;

pub use curve::Curve;
pub use nurbs::{Nurbs1, Nurbs2, expand_knots};
pub use surface::Surface;
