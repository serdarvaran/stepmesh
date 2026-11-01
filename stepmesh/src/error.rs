use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// The input is not an ISO 10303-21 exchange structure.
    NotStep,
    Syntax {
        line: usize,
        msg: String,
    },
    /// The file uses an entity this crate does not read.
    Unsupported {
        entity: String,
        id: u64,
    },
    /// An instance does not have the structure its type requires.
    Invalid {
        id: u64,
        msg: String,
    },
    /// The file has no solid geometry (for example only references to external documents).
    NoGeometry,
    /// A shell is open, so the shape has no volume.
    OpenShell {
        id: u64,
    },
    /// The mesh of a solid is not a closed, consistently oriented surface.
    NotWatertight {
        solid: usize,
    },
    TriangleLimit {
        count: usize,
        limit: usize,
    },
    /// A point of a face boundary could not be mapped onto the face's surface.
    Projection {
        face: u64,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NotStep => write!(f, "not a STEP (ISO 10303-21) file"),
            Error::Syntax { line, msg } => write!(f, "syntax error on line {line}: {msg}"),
            Error::Unsupported { entity, id } => write!(f, "unsupported entity {entity} (#{id})"),
            Error::Invalid { id, msg } => write!(f, "invalid instance #{id}: {msg}"),
            Error::NoGeometry => write!(f, "no solid geometry"),
            Error::OpenShell { id } => write!(f, "shell #{id} is open"),
            Error::NotWatertight { solid } => write!(f, "mesh of solid {solid} is not watertight"),
            Error::TriangleLimit { count, limit } => {
                write!(f, "{count} triangles exceed the limit of {limit}")
            }
            Error::Projection { face } => {
                write!(f, "boundary of face #{face} does not lie on its surface")
            }
        }
    }
}

impl std::error::Error for Error {}
