use std::collections::HashMap;

/// One parameter value of an entity instance.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Int(i64),
    Real(f64),
    Str(String),
    /// Enumeration or logical (`.T.`, `.F.`, `.U.`), without the dots.
    Enum(String),
    /// Entity instance reference `#n`.
    Ref(u64),
    List(Vec<Value>),
    /// Typed parameter such as `LENGTH_MEASURE(25.4)`.
    Typed(String, Box<Value>),
    /// `$`
    Unset,
    /// `*`
    Derived,
    /// Binary literal `"…"`, as written (hex digits after the unused-bits digit).
    Binary(String),
}

/// One `NAME(params)` record: a simple instance has one, a complex instance several.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// Upper case.
    pub name: String,
    pub params: Vec<Value>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Entity {
    pub records: Vec<Record>,
}

#[derive(Debug, Clone, Default)]
pub struct Exchange {
    /// `FILE_SCHEMA` names.
    pub schemas: Vec<String>,
    /// Instances of every DATA section, by `#id`.
    pub entities: HashMap<u64, Entity>,
}
