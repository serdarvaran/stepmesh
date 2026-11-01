use stepmesh::Error;
use stepmesh::p21::{Value, parse};

fn file(data: &str) -> String {
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME('t','',(''),(''),'','','');\n\
         FILE_SCHEMA(('AUTOMOTIVE_DESIGN {{ 1 0 10303 214 1 1 1 1 }}'));\nENDSEC;\nDATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;\n"
    )
}

fn single(data: &str) -> Vec<Value> {
    let x = parse(file(data).as_bytes()).unwrap();
    x.entities[&1].records[0].params.clone()
}

#[test]
fn minimal_file_parses() {
    let x = parse(file("#1=CARTESIAN_POINT('',(0.,1.5,-2.E3));").as_bytes()).unwrap();
    assert_eq!(x.schemas, ["AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }"]);
    let r = &x.entities[&1].records[0];
    assert_eq!(r.name, "CARTESIAN_POINT");
    assert_eq!(
        r.params,
        [
            Value::Str(String::new()),
            Value::List(vec![
                Value::Real(0.0),
                Value::Real(1.5),
                Value::Real(-2000.0)
            ])
        ]
    );
}

#[test]
fn string_escapes() {
    let p = single(r"#1=A('it''s','\X2\00E700FC\X0\','\X\E7','\S\a');");
    assert_eq!(
        p,
        [
            Value::Str("it's".into()),
            Value::Str("çü".into()),
            Value::Str("ç".into()),
            Value::Str("á".into()),
        ]
    );
}

#[test]
fn complex_instance() {
    let x =
        parse(file("#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));").as_bytes()).unwrap();
    let e = &x.entities[&1];
    assert_eq!(e.records.len(), 3);
    let si = e.records.iter().find(|r| r.name == "SI_UNIT").unwrap();
    assert_eq!(
        si.params,
        [Value::Enum("MILLI".into()), Value::Enum("METRE".into())]
    );
    assert_eq!(e.records[1].params, [Value::Derived]);
}

#[test]
fn typed_and_special_values() {
    let p = single("#1=A(LENGTH_MEASURE(25.4),$,*,.T.,#12,12,\"0FF\");");
    assert_eq!(
        p,
        [
            Value::Typed("LENGTH_MEASURE".into(), Box::new(Value::Real(25.4))),
            Value::Unset,
            Value::Derived,
            Value::Enum("T".into()),
            Value::Ref(12),
            Value::Int(12),
            Value::Binary("0FF".into()),
        ]
    );
}

#[test]
fn comments_and_whitespace_ignored() {
    let x =
        parse(file("/* a */ #1 = A ( /* b */ 1 , \n 2 ) ; /* c */\n#2=B();").as_bytes()).unwrap();
    assert_eq!(
        x.entities[&1].records[0].params,
        [Value::Int(1), Value::Int(2)]
    );
    assert!(x.entities[&2].records[0].params.is_empty());
}

#[test]
fn not_step_is_rejected() {
    assert!(matches!(parse(b"solid x"), Err(Error::NotStep)));
    assert!(matches!(parse(b""), Err(Error::NotStep)));
}

#[test]
fn duplicate_id_is_an_error() {
    assert!(matches!(
        parse(file("#1=A();\n#1=B();").as_bytes()),
        Err(Error::Syntax { .. })
    ));
}

#[test]
fn unterminated_string_reports_line() {
    // `DATA;` is line 7 of `file`, so #3 with the open string is on line 10.
    let r = parse(file("#1=A();\n#2=B();\n#3=C('open);").as_bytes());
    assert!(matches!(r, Err(Error::Syntax { line: 10, .. })), "{r:?}");
}

#[test]
fn nesting_deeper_than_64_is_an_error() {
    let deep = format!("#1=A({}1{});", "(".repeat(65), ")".repeat(65));
    assert!(matches!(
        parse(file(&deep).as_bytes()),
        Err(Error::Syntax { .. })
    ));
    let ok = format!("#1=A({}1{});", "(".repeat(63), ")".repeat(63));
    assert!(parse(file(&ok).as_bytes()).is_ok());
}

#[test]
fn integer_overflow_is_an_error() {
    assert!(matches!(
        parse(file("#1=A(99999999999999999999);").as_bytes()),
        Err(Error::Syntax { .. })
    ));
}

#[test]
fn every_prefix_is_rejected_or_parsed() {
    let text = file(
        "#1=CARTESIAN_POINT('p\\X2\\00E7\\X0\\',(0.,1.5,-2.E3));\n\
         #2=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n\
         #3=A(LENGTH_MEASURE(1.),$,*,#2,(#1,#2),\"0F\");",
    );
    for n in 0..=text.len() {
        let _ = parse(&text.as_bytes()[..n]);
    }
    assert!(parse(text.as_bytes()).is_ok());
}
