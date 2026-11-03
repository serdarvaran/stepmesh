//! Curve and surface evaluation.

pub(crate) mod bounds;
mod curve;
mod nurbs;
mod surface;

pub use curve::Curve;
pub use nurbs::{MAX_DEGREE, Nurbs1, Nurbs2, expand_knots};
pub use surface::Surface;
