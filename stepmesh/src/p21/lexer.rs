use super::value::Value;
use crate::Error;

const MAX_DEPTH: usize = 64;

/// Byte cursor over a Part 21 file; skips whitespace and comments between tokens.
pub(super) struct Cursor<'a> {
    b: &'a [u8],
    pub pos: usize,
    pub line: usize,
}

impl<'a> Cursor<'a> {
    pub fn new(b: &'a [u8]) -> Self {
        Cursor { b, pos: 0, line: 1 }
    }

    pub fn err<T>(&self, msg: impl Into<String>) -> Result<T, Error> {
        Err(Error::Syntax {
            line: self.line,
            msg: msg.into(),
        })
    }

    fn bump(&mut self) -> Option<u8> {
        let c = *self.b.get(self.pos)?;
        self.pos += 1;
        if c == b'\n' {
            self.line += 1;
        }
        Some(c)
    }

    pub fn skip_ws(&mut self) -> Result<(), Error> {
        loop {
            match self.b.get(self.pos) {
                Some(c) if c.is_ascii_whitespace() => {
                    self.bump();
                }
                Some(b'/') if self.b.get(self.pos + 1) == Some(&b'*') => {
                    self.pos += 2;
                    loop {
                        match self.bump() {
                            None => return self.err("unterminated comment"),
                            Some(b'*') if self.b.get(self.pos) == Some(&b'/') => {
                                self.pos += 1;
                                break;
                            }
                            _ => {}
                        }
                    }
                }
                _ => return Ok(()),
            }
        }
    }

    pub fn peek(&mut self) -> Result<Option<u8>, Error> {
        self.skip_ws()?;
        Ok(self.b.get(self.pos).copied())
    }

    pub fn at_end(&mut self) -> Result<bool, Error> {
        Ok(self.peek()?.is_none())
    }

    pub fn expect(&mut self, c: u8) -> Result<(), Error> {
        match self.peek()? {
            Some(x) if x == c => {
                self.bump();
                Ok(())
            }
            Some(x) => self.err(format!("expected '{}', found '{}'", c as char, x as char)),
            None => self.err(format!("expected '{}', found end of file", c as char)),
        }
    }

    /// Accepts `c` if it is the next token.
    pub fn eat(&mut self, c: u8) -> Result<bool, Error> {
        if self.peek()? == Some(c) {
            self.bump();
            return Ok(true);
        }
        Ok(false)
    }

    /// Standard or user-defined (`!NAME`) keyword; `-` is allowed for `ISO-10303-21` and its end mark.
    pub fn keyword(&mut self) -> Result<String, Error> {
        self.skip_ws()?;
        let start = self.pos;
        if self.b.get(self.pos) == Some(&b'!') {
            self.pos += 1;
        }
        while let Some(&c) = self.b.get(self.pos) {
            if c.is_ascii_alphanumeric() || c == b'_' || c == b'-' {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return self.err("expected a keyword");
        }
        Ok(String::from_utf8_lossy(&self.b[start..self.pos]).to_ascii_uppercase())
    }

    fn integer(&mut self) -> Result<u64, Error> {
        let start = self.pos;
        while self.b.get(self.pos).is_some_and(u8::is_ascii_digit) {
            self.pos += 1;
        }
        let text = std::str::from_utf8(&self.b[start..self.pos]).unwrap_or("");
        text.parse()
            .or_else(|_| self.err(format!("bad integer '{text}'")))
    }

    /// `#n` after `#` has been seen.
    pub fn instance_id(&mut self) -> Result<u64, Error> {
        self.expect(b'#')?;
        self.integer()
    }

    fn number(&mut self) -> Result<Value, Error> {
        let start = self.pos;
        let mut real = false;
        if matches!(self.b.get(self.pos), Some(b'+' | b'-')) {
            self.pos += 1;
        }
        while let Some(&c) = self.b.get(self.pos) {
            match c {
                b'0'..=b'9' => self.pos += 1,
                b'.' => {
                    real = true;
                    self.pos += 1;
                }
                b'E' | b'e' => {
                    real = true;
                    self.pos += 1;
                    if matches!(self.b.get(self.pos), Some(b'+' | b'-')) {
                        self.pos += 1;
                    }
                }
                _ => break,
            }
        }
        let text = std::str::from_utf8(&self.b[start..self.pos]).unwrap_or("");
        if real {
            let t = text.trim_end_matches(['E', 'e']);
            t.parse()
                .map(Value::Real)
                .or_else(|_| self.err(format!("bad real '{text}'")))
        } else {
            text.parse()
                .map(Value::Int)
                .or_else(|_| self.err(format!("bad integer '{text}'")))
        }
    }

    fn enumeration(&mut self) -> Result<Value, Error> {
        self.bump(); // '.'
        let start = self.pos;
        while self
            .b
            .get(self.pos)
            .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
        {
            self.pos += 1;
        }
        let name = String::from_utf8_lossy(&self.b[start..self.pos]).to_ascii_uppercase();
        if self.bump() != Some(b'.') || name.is_empty() {
            return self.err("bad enumeration");
        }
        Ok(Value::Enum(name))
    }

    fn binary(&mut self) -> Result<Value, Error> {
        self.bump(); // '"'
        let start = self.pos;
        while self.b.get(self.pos).is_some_and(u8::is_ascii_hexdigit) {
            self.pos += 1;
        }
        let text = String::from_utf8_lossy(&self.b[start..self.pos]).into_owned();
        if self.bump() != Some(b'"') {
            return self.err("unterminated binary");
        }
        Ok(Value::Binary(text))
    }

    fn hex(&mut self, digits: usize) -> Result<u32, Error> {
        let s = self
            .b
            .get(self.pos..self.pos + digits)
            .ok_or(())
            .ok()
            .and_then(|s| std::str::from_utf8(s).ok())
            .and_then(|s| u32::from_str_radix(s, 16).ok());
        match s {
            Some(v) => {
                self.pos += digits;
                Ok(v)
            }
            None => self.err("bad hex escape"),
        }
    }

    fn starts_with(&self, s: &[u8]) -> bool {
        self.b[self.pos..].starts_with(s)
    }

    /// String literal with its `''` and backslash directives decoded.
    fn string(&mut self) -> Result<Value, Error> {
        let start_line = self.line;
        self.bump(); // '\''
        let mut out: Vec<u8> = Vec::new();
        let push = |out: &mut Vec<u8>, ch: u32| {
            let c = char::from_u32(ch).unwrap_or(char::REPLACEMENT_CHARACTER);
            let mut buf = [0; 4];
            out.extend_from_slice(c.encode_utf8(&mut buf).as_bytes());
        };
        loop {
            let Some(c) = self.bump() else {
                return Err(Error::Syntax {
                    line: start_line,
                    msg: "unterminated string".into(),
                });
            };
            match c {
                b'\'' if self.b.get(self.pos) == Some(&b'\'') => {
                    self.pos += 1;
                    out.push(b'\'');
                }
                b'\'' => break,
                b'\n' | b'\r' => {} // line breaks inside strings are not part of the value
                b'\\' if self.starts_with(b"\\") => {
                    self.pos += 1;
                    out.push(b'\\');
                }
                b'\\' if self.starts_with(b"X\\") => {
                    self.pos += 2;
                    let v = self.hex(2)?;
                    push(&mut out, v);
                }
                b'\\' if self.starts_with(b"X2\\") || self.starts_with(b"X4\\") => {
                    let digits = if self.b[self.pos + 1] == b'2' { 4 } else { 8 };
                    self.pos += 3;
                    let mut units = Vec::new();
                    while !self.starts_with(b"\\X0\\") {
                        if self.pos >= self.b.len() {
                            return self.err("unterminated \\X2\\ or \\X4\\ escape");
                        }
                        units.push(self.hex(digits)?);
                    }
                    self.pos += 4;
                    if digits == 4 {
                        let u16s: Vec<u16> = units.iter().map(|&u| u as u16).collect();
                        for ch in char::decode_utf16(u16s) {
                            push(&mut out, ch.map_or(0xFFFD, u32::from));
                        }
                    } else {
                        for u in units {
                            push(&mut out, u);
                        }
                    }
                }
                b'\\' if self.starts_with(b"S\\") => {
                    self.pos += 2;
                    let Some(ch) = self.bump() else {
                        return self.err("unterminated string");
                    };
                    push(&mut out, u32::from(ch) + 128);
                }
                b'\\'
                    if self.b.get(self.pos) == Some(&b'P')
                        && self.b.get(self.pos + 2) == Some(&b'\\') =>
                {
                    self.pos += 3; // code page switch: \S\ is always read as ISO 8859-1
                }
                b'\\' if self.starts_with(b"N\\") => {
                    self.pos += 2;
                    out.push(b'\n');
                }
                _ => out.push(c),
            }
        }
        Ok(Value::Str(String::from_utf8_lossy(&out).into_owned()))
    }

    /// Comma-separated values up to `)`, after the `(` has been read.
    pub fn params(&mut self, depth: usize) -> Result<Vec<Value>, Error> {
        if depth > MAX_DEPTH {
            return self.err("lists nested too deeply");
        }
        let mut out = Vec::new();
        if self.eat(b')')? {
            return Ok(out);
        }
        loop {
            out.push(self.value(depth)?);
            if self.eat(b')')? {
                return Ok(out);
            }
            self.expect(b',')?;
        }
    }

    fn value(&mut self, depth: usize) -> Result<Value, Error> {
        match self.peek()? {
            None => self.err("unexpected end of file"),
            Some(b'$') => {
                self.bump();
                Ok(Value::Unset)
            }
            Some(b'*') => {
                self.bump();
                Ok(Value::Derived)
            }
            Some(b'#') => Ok(Value::Ref(self.instance_id()?)),
            Some(b'\'') => self.string(),
            Some(b'"') => self.binary(),
            Some(b'.') if self.b.get(self.pos + 1).is_some_and(u8::is_ascii_digit) => self.number(),
            Some(b'.') => self.enumeration(),
            Some(b'(') => {
                self.bump();
                Ok(Value::List(self.params(depth + 1)?))
            }
            Some(c) if c.is_ascii_digit() || c == b'+' || c == b'-' => self.number(),
            Some(c) if c.is_ascii_alphabetic() || c == b'!' => {
                let name = self.keyword()?;
                self.expect(b'(')?;
                let mut inner = self.params(depth + 1)?;
                if inner.len() != 1 {
                    return self.err(format!("typed parameter {name} needs one value"));
                }
                Ok(Value::Typed(
                    name,
                    Box::new(inner.pop().unwrap_or(Value::Unset)),
                ))
            }
            Some(c) => self.err(format!("unexpected '{}'", c as char)),
        }
    }
}
