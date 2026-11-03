# Changelog

## 0.1.0 (unreleased)

First release: STEP (ISO 10303-21) files to closed, oriented triangle meshes in millimetres.

- Part 21 parser (simple and complex instances, string escapes, all DATA sections).
- Units per representation context, assemblies (`NEXT_ASSEMBLY_USAGE_OCCURRENCE`, `MAPPED_ITEM`).
- Analytic, B-spline (explicit, uniform, quasi-uniform and Bézier knots; rational), extrusion,
  revolution and offset surfaces; lines, conics and B-spline curves.
- Tessellation in each face's parameter space with shared edge samples, seams, poles, apices and
  seamless periodic faces (seams winding between holes); refinement to the chord and angle
  tolerances; edges sampled finer where their polylines cross in a face.
- Edges running past the end of a B-spline or swept surface within the file's tolerance, and
  faces written against their surface's normal, are read as their neighbours meet them.
- Triangles of B-spline, swept and offset faces are certified within the deflection by bounds
  from the control net, and edges of spline faces get half of it; a boundary farther from its face
  than the deflection, or a face that cannot be certified, is an error.
- Every solid is checked to be closed and consistently oriented; anything else is an error.
