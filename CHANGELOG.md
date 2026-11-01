# Changelog

## 0.1.0 (unreleased)

First release: STEP (ISO 10303-21) files to closed, oriented triangle meshes in millimetres.

- Part 21 parser (simple and complex instances, string escapes, all DATA sections).
- Units per representation context, assemblies (`NEXT_ASSEMBLY_USAGE_OCCURRENCE`, `MAPPED_ITEM`).
- Analytic, B-spline (explicit, uniform, quasi-uniform and Bézier knots; rational), extrusion and
  revolution surfaces; lines, conics and B-spline curves.
- Tessellation in each face's parameter space with shared edge samples, seams, poles, apices and
  seamless periodic faces; refinement to the chord and angle tolerances.
- Every solid is checked to be closed and consistently oriented; anything else is an error.
