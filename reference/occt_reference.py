"""Reference measurements with OpenCascade (OCP): one JSON line per STEP file.

    python reference/occt_reference.py <list-of-paths.txt | file.step ...> <out.jsonl>

Each line: {file, ok, solids, volume, area, min, max, valid} in millimetres, with the exact B-rep
volume and area (BRepGProp), or {file, ok: false, error}. Solids are counted after the reader has
applied assembly placements, so an instanced part counts once per instance.
"""

import json
import sys
from pathlib import Path

from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer


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
    return {"ok": True, "solids": len(solids), "volume": volume, "area": area,
            "min": [xmin, ymin, zmin], "max": [xmax, ymax, zmax], "valid": valid}


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
