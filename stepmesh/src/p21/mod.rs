//! ISO 10303-21 ("Part 21") exchange structure.

mod lexer;
mod value;

pub use value::{Entity, Exchange, Record, Value};

use crate::Error;
use lexer::Cursor;
use std::collections::hash_map::Entry;

/// Parses the header (for `FILE_SCHEMA`) and every DATA section.
pub fn parse(bytes: &[u8]) -> Result<Exchange, Error> {
    let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    let mut c = Cursor::new(bytes);
    if c.keyword().ok().as_deref() != Some("ISO-10303-21") || !c.eat(b';').unwrap_or(false) {
        return Err(Error::NotStep);
    }
    let mut x = Exchange::default();
    let mut has_data = false;
    loop {
        if c.at_end()? {
            // A file cut after its last section is still read; one cut inside a section is not.
            return if has_data {
                Ok(x)
            } else {
                c.err("missing DATA section")
            };
        }
        let section = c.keyword()?;
        match section.as_str() {
            "END-ISO-10303-21" => {
                c.expect(b';')?;
                return Ok(x);
            }
            "HEADER" => {
                c.expect(b';')?;
                header(&mut c, &mut x)?;
            }
            "DATA" => {
                if c.eat(b'(')? {
                    c.params(1)?;
                }
                c.expect(b';')?;
                data(&mut c, &mut x)?;
                has_data = true;
            }
            "ANCHOR" | "REFERENCE" | "SIGNATURE" => {
                return c.err(format!("unsupported section {section}"));
            }
            _ => return c.err(format!("unexpected '{section}'")),
        }
    }
}

fn header(c: &mut Cursor, x: &mut Exchange) -> Result<(), Error> {
    loop {
        let name = c.keyword()?;
        if name == "ENDSEC" {
            return c.expect(b';');
        }
        c.expect(b'(')?;
        let params = c.params(1)?;
        c.expect(b';')?;
        if name == "FILE_SCHEMA" {
            if let Some(Value::List(names)) = params.first() {
                x.schemas = names
                    .iter()
                    .filter_map(|v| match v {
                        Value::Str(s) => Some(s.clone()),
                        _ => None,
                    })
                    .collect();
            }
        }
    }
}

fn data(c: &mut Cursor, x: &mut Exchange) -> Result<(), Error> {
    loop {
        if c.peek()? != Some(b'#') {
            let end = c.keyword()?;
            if end != "ENDSEC" {
                return c.err(format!("expected an instance or ENDSEC, found '{end}'"));
            }
            return c.expect(b';');
        }
        let line = c.line;
        let id = c.instance_id()?;
        c.expect(b'=')?;
        let mut records = Vec::new();
        if c.eat(b'(')? {
            while !c.eat(b')')? {
                records.push(record(c)?);
            }
            if records.is_empty() {
                return c.err("empty complex instance");
            }
        } else {
            records.push(record(c)?);
        }
        c.expect(b';')?;
        match x.entities.entry(id) {
            Entry::Occupied(_) => {
                return Err(Error::Syntax {
                    line,
                    msg: format!("#{id} defined twice"),
                });
            }
            Entry::Vacant(v) => {
                v.insert(Entity { records });
            }
        }
    }
}

fn record(c: &mut Cursor) -> Result<Record, Error> {
    let name = c.keyword()?;
    c.expect(b'(')?;
    let params = c.params(1)?;
    Ok(Record { name, params })
}
