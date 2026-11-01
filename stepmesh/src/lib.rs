//! Reads STEP (ISO 10303-21) files and turns their B-rep solids into watertight triangle meshes.
#![forbid(unsafe_code)]

mod error;
pub mod p21;

pub use error::Error;
