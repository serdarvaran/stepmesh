"""Generates tests/fixtures/*.step and fixtures.json with OpenCascade (CadQuery + OCP).

Volumes and areas in fixtures.json are OpenCascade's exact B-rep properties (mm, mm²);
where a closed form exists the script checks it first.

    python reference/make_fixtures.py   # with CadQuery installed
"""

import json
import math
import re
from pathlib import Path

import cadquery as cq
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer

OUT = Path(__file__).resolve().parent.parent / "stepmesh" / "tests" / "fixtures"


def props(shape):
    v, a = GProp_GProps(), GProp_GProps()
    BRepGProp.VolumeProperties_s(shape.wrapped, v)
    BRepGProp.SurfaceProperties_s(shape.wrapped, a)
    return v.Mass(), a.Mass()


def write(shape, name):
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    w = STEPControl_Writer()
    w.Transfer(shape.wrapped, STEPControl_AsIs)
    assert w.Write(str(OUT / name)) == IFSelect_RetDone, name
    normalize(OUT / name)


def normalize(path):
    """Fixed header timestamp and author so regenerating gives the same bytes."""
    text = path.read_text()
    text = re.sub(r"FILE_NAME\('[^']*','[^']*'", "FILE_NAME('fixture','2026-01-01T00:00:00'", text, count=1)
    path.write_text(text)


def spline_profile():
    pts = [(0, 0), (8, 1), (12, 6), (9, 12), (3, 10), (0, 6)]
    return cq.Workplane("XY").spline(pts, includeCurrent=False, periodic=True).close()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fixtures = {}

    def add(name, shape, analytic=None, solids=1):
        if isinstance(shape, cq.Workplane):
            shape = shape.val()
        volume, area = props(shape)
        if analytic is not None:
            assert abs(volume - analytic) <= 1e-6 * analytic, (name, volume, analytic)
        write(shape, name)
        fixtures[name] = {"volume": volume, "area": area, "solids": solids}

    add("box.step", cq.Workplane().box(10, 20, 30), 6000)
    add("box_inch.step", cq.Workplane().box(1, 1, 1), 1)
    to_inch(OUT / "box_inch.step", whole_file=True)
    fixtures["box_inch.step"] = {"volume": 25.4**3, "area": 6 * 25.4**2, "solids": 1}
    add("box_with_hole.step", cq.Workplane().box(20, 20, 10).faces(">Z").workplane().hole(5),
        4000 - math.pi * 2.5**2 * 10)
    outer = cq.Workplane().box(20, 20, 20).val()
    inner = cq.Workplane().box(10, 10, 10).val()
    add("shelled_box.step", outer.cut(inner), 7000)
    add("cylinder.step", cq.Workplane().cylinder(10, 5), math.pi * 25 * 10)
    big = cq.Workplane().cylinder(20, 10).val()
    drill = cq.Workplane("YZ").cylinder(40, 2).val()  # along X: crosses the seam at u = 0
    add("cylinder_hole_seam.step", big.cut(drill))
    add("cone.step", cq.Solid.makeCone(5, 0, 10), math.pi * 25 * 10 / 3)
    add("sphere.step", cq.Solid.makeSphere(7, angleDegrees1=-90, angleDegrees2=90), 4 / 3 * math.pi * 7**3)
    add("torus.step", cq.Solid.makeTorus(10, 3), 2 * math.pi**2 * 10 * 9)
    add("revolve_spline.step", cq.Workplane("XZ").spline([(2, 0), (6, 3), (4, 8), (3, 12)], includeCurrent=False)
        .lineTo(0, 12).lineTo(0, 0).close().revolve(360, (0, 0, 0), (0, 1, 0)))
    add("extrude_spline.step", spline_profile().extrude(5))
    add("loft.step", cq.Workplane("XY").rect(10, 10).workplane(offset=10).circle(4).loft())
    add("fillet_box.step", cq.Workplane().box(20, 20, 20).edges().fillet(3))

    part = cq.Workplane().box(10, 20, 30)
    asm = cq.Assembly(name="pair")
    asm.add(part, name="a", loc=cq.Location((0, 0, 0)))
    asm.add(part, name="b", loc=cq.Location((50, 0, 0)))
    asm.save(str(OUT / "assembly_two_boxes.step"))
    normalize(OUT / "assembly_two_boxes.step")
    fixtures["assembly_two_boxes.step"] = {"volume": 12000.0, "area": 2 * 2200.0, "solids": 2}

    cube = cq.Workplane().box(1, 1, 1)
    asm = cq.Assembly(name="mixed")
    asm.add(cube, name="inch_part", loc=cq.Location((30, 0, 0)))
    asm.save(str(OUT / "mixed_units.step"))
    to_inch(OUT / "mixed_units.step", whole_file=False)
    normalize(OUT / "mixed_units.step")
    fixtures["mixed_units.step"] = {"volume": 25.4**3, "area": 6 * 25.4**2, "solids": 1}

    (OUT / "external_reference.step").write_text(EXTERNAL_REFERENCE)
    fixtures["external_reference.step"] = {"error": "NoGeometry"}

    (OUT / "fixtures.json").write_text(json.dumps(fixtures, indent=2, sort_keys=True) + "\n")


def to_inch(path, whole_file):
    """Makes the length unit of the solid's representation context the inch.

    With `whole_file` the context itself is changed (every representation using it is in inches);
    otherwise only the solid's representation gets an inch copy of it, so an assembly around it
    keeps millimetres.
    """
    text = path.read_text()
    ents = dict(re.findall(r"^#(\d+)\s*=\s*(.*?);\s*$", text, re.M | re.S))
    brep_rep = next(i for i, e in ents.items()
                    if e.startswith("ADVANCED_BREP_SHAPE_REPRESENTATION")
                    and any(ents[r].startswith("MANIFOLD_SOLID_BREP") for r in re.findall(r"#(\d+)", e)))
    ctx = re.findall(r"#(\d+)", ents[brep_rep])[-1]
    units = re.search(r"GLOBAL_UNIT_ASSIGNED_CONTEXT\s*\(\s*\(([^)]*)\)", ents[ctx]).group(1)
    mm_unit = next(u for u in re.findall(r"#(\d+)", units) if "LENGTH_UNIT" in ents[u])
    n = max(map(int, ents)) + 1
    inch, measure, dims, new_ctx = n, n + 1, n + 2, n + 3
    ents[str(dims)] = "DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)"
    ents[str(inch)] = f"( CONVERSION_BASED_UNIT('INCH',#{measure}) LENGTH_UNIT() NAMED_UNIT(#{dims}) )"
    ents[str(measure)] = f"LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#{mm_unit})"
    inch_ctx = ents[ctx].replace(units, re.sub(rf"#{mm_unit}\b", f"#{inch}", units))
    if whole_file:
        ents[ctx] = inch_ctx
    else:
        ents[str(new_ctx)] = inch_ctx
        ents[brep_rep] = re.sub(rf"#{ctx}(\s*\)\s*)$", rf"#{new_ctx}\1", ents[brep_rep])
    body = "".join(f"#{i} = {ents[i]};\n" for i in sorted(ents, key=int))
    text = re.sub(r"DATA;.*ENDSEC;", lambda _: "DATA;\n" + body + "ENDSEC;", text, count=1, flags=re.S)
    path.write_text(text)


EXTERNAL_REFERENCE = """ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('fixture','2026-01-01T00:00:00',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
ENDSEC;
DATA;
#1=APPLICATION_CONTEXT('automotive design');
#2=PRODUCT_CONTEXT('',#1,'mechanical');
#3=PRODUCT('part','part','',(#2));
#4=PRODUCT_DEFINITION_FORMATION('','',#3);
#5=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');
#6=PRODUCT_DEFINITION('design','',#4,#5);
#7=DOCUMENT_TYPE('');
#8=DOCUMENT_FILE('part.stp','','',#7,'',$);
#9=APPLIED_DOCUMENT_REFERENCE(#8,'',(#6));
ENDSEC;
END-ISO-10303-21;
"""

if __name__ == "__main__":
    main()
