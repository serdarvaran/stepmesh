"""Reference measurements with OpenCascade (OCP): one JSON line per STEP file.

    python reference/occt_reference.py <list-of-paths.txt | file.step ...> <out.jsonl>

Each line: {file, ok, solids, volume, mesh_volume, area, min, max, valid} in millimetres, with the
exact B-rep volume and area (BRepGProp), the volume of OpenCascade's own mesh at stepmesh's default
tolerances (0.003 mm, 0.5 rad), or {file, ok: false, error}. Solids are counted after the reader has
applied assembly placements, so an instanced part counts once per instance. The bounds come from
a fine OpenCascade mesh (deflection 1e-4 of the diagonal, at least 0.01 mm): Bnd_Box, even the
"optimal" one, can reach past trimmed B-spline faces.
"""

import json
import sys
from pathlib import Path

from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRep import BRep_Tool
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS


def measure(path):
    Interface_Static.SetCVal_s("xstep.cascade.unit", "MM")
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone:
        return {"ok": False, "error": "read failed"}
    reader.TransferRoots()
    shape = reader.OneShape()
    solids = []
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        solids.append(explorer.Current())
        explorer.Next()
    if not solids:
        return {"ok": False, "error": "no solids"}
    volume = area = 0.0
    valid = True
    box = Bnd_Box()
    for s in solids:
        v, a = GProp_GProps(), GProp_GProps()
        BRepGProp.VolumeProperties_s(s, v)
        BRepGProp.SurfaceProperties_s(s, a)
        volume += v.Mass()
        area += a.Mass()
        valid = valid and BRepCheck_Analyzer(s).IsValid()
        BRepBndLib.AddOptimal_s(s, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    lo, hi = mesh_bounds(shape, max(0.01, 1e-4 * ((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2) ** 0.5))
    if lo is None:
        lo, hi = [xmin, ymin, zmin], [xmax, ymax, zmax]
    return {"ok": True, "solids": len(solids), "volume": volume, "mesh_volume": mesh_volume(shape),
            "area": area, "min": lo, "max": hi, "valid": valid}


def mesh_volume(shape):
    """Signed volume of OpenCascade's mesh at 0.003 mm / 0.5 rad, the tolerances stepmesh uses."""
    BRepMesh_IncrementalMesh(shape, 0.003, False, 0.5, True)
    total = 0.0
    faces = TopExp_Explorer(shape, TopAbs_FACE)
    while faces.More():
        face = TopoDS.Face_s(faces.Current())
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            return None
        trsf = loc.Transformation()
        for i in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(i).Get()
            if face.Orientation() == TopAbs_REVERSED:
                b, c = c, b
            (ax, ay, az), (bx, by, bz), (cx, cy, cz) = (
                (q.X(), q.Y(), q.Z()) for q in (tri.Node(k).Transformed(trsf) for k in (a, b, c)))
            total += (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6.0
        faces.Next()
    return total


def mesh_bounds(shape, deflection):
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.5, True)
    lo, hi = [float("inf")] * 3, [float("-inf")] * 3
    faces = TopExp_Explorer(shape, TopAbs_FACE)
    while faces.More():
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(TopoDS.Face_s(faces.Current()), loc)
        if tri is not None:
            trsf = loc.Transformation()
            for i in range(1, tri.NbNodes() + 1):
                p = tri.Node(i).Transformed(trsf)
                for k, c in enumerate((p.X(), p.Y(), p.Z())):
                    lo[k], hi[k] = min(lo[k], c), max(hi[k], c)
        faces.Next()
    return (lo, hi) if lo[0] != float("inf") else (None, None)


def main():
    *inputs, out = sys.argv[1:]
    files = []
    for i in inputs:
        p = Path(i)
        files += p.read_text().split("\n") if p.suffix == ".txt" else [i]
    with open(out, "w") as sink:
        for f in filter(None, files):
            try:
                row = measure(f)
            except Exception as e:  # the reference must not stop the run on one bad file
                row = {"ok": False, "error": f"{type(e).__name__}: {e}"}
            sink.write(json.dumps({"file": f, **row}) + "\n")
            sink.flush()


if __name__ == "__main__":
    main()
