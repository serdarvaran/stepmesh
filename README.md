# stepmesh

Pure Rust reader for STEP files (ISO 10303-21, the geometry of AP203, AP214 and AP242) that turns
their B-rep solids into closed, consistently oriented triangle meshes in millimetres. No C or C++
code: the only dependency is [`spade`](https://crates.io/crates/spade) for constrained Delaunay
triangulation.

```toml
[dependencies]
stepmesh = "0.1"
```

```rust,no_run
let bytes = std::fs::read("part.step")?;
let mesh = stepmesh::read_step(&bytes, &stepmesh::Options::default())?;
for solid in &mesh.solids {
    println!("{:?}: {:.1} mm³, {} triangles", solid.name, solid.volume(), solid.triangles.len());
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Guarantees

- Every returned solid is a closed, consistently oriented triangle surface (each directed edge used
  once, its reverse once); outer shells enclose positive volume, voids negative.
- Coordinates are millimetres whatever the file's units, per representation (an inch part in a
  millimetre assembly is converted on its own), with assembly placements applied: one solid per
  placed instance.
- Triangles stay within `Options::linear_deflection` (default 0.003 mm) of their surfaces, and the
  normals of neighbouring vertices within `angular_deflection` (default 0.5 rad).
- A file this crate cannot mesh faithfully fails with an [`Error`] naming the reason (unsupported
  entity, open shell, a face that cannot be laid out in its parameter space, a mesh that would not
  close). It never returns an approximate or partial result.

## Supported entities

| Area | Entities |
|---|---|
| File | ISO 10303-21 exchange structure, simple and complex instances, all DATA sections, string escapes (`\X\`, `\X2\`, `\X4\`, `\S\`) |
| Units | `SI_UNIT` with prefixes, `CONVERSION_BASED_UNIT` (inch, foot, degree, …), per representation context |
| Product structure | `SHAPE_DEFINITION_REPRESENTATION`, `SHAPE_REPRESENTATION_RELATIONSHIP`, `NEXT_ASSEMBLY_USAGE_OCCURRENCE` with `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION` and `ITEM_DEFINED_TRANSFORMATION`, `MAPPED_ITEM`; files without product structure |
| Solids | `MANIFOLD_SOLID_BREP`, `BREP_WITH_VOIDS`, `CLOSED_SHELL`, `ORIENTED_CLOSED_SHELL`, `ADVANCED_FACE`, `FACE_SURFACE`, `FACE_BOUND`, `FACE_OUTER_BOUND`, `EDGE_LOOP`, `VERTEX_LOOP`, `ORIENTED_EDGE`, `EDGE_CURVE`, `VERTEX_POINT` |
| Surfaces | `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`, B-spline surfaces (explicit knots, uniform, quasi-uniform, Bézier; rational), `SURFACE_OF_LINEAR_EXTRUSION`, `SURFACE_OF_REVOLUTION`, `OFFSET_SURFACE` |
| Curves | `LINE`, `CIRCLE`, `ELLIPSE`, B-spline curves (same forms), `SURFACE_CURVE`, `SEAM_CURVE`, `INTERSECTION_CURVE`, `TRIMMED_CURVE` |

Not read yet, each an explicit error: `RECTANGULAR_TRIMMED_SURFACE`, AP242
tessellated shapes, faceted B-reps and CSG solids, any other representation item that is not
wireframe or annotation, open shells and surface models (no volume), assembly components without
a shape (references to external documents), mirroring or scaling transformations, B-splines above
degree 16. Malformed input (cyclic references, inconsistent B-splines, out-of-range values, more
than 100 000 placed solids) is an `Error` too, never a panic or a hang.

## How it meshes

Each edge is sampled once, to the chord and angle tolerances, and both faces using it share those
samples. Each face is laid out in its surface's (u, v) parameter space: boundary points are
projected onto the surface, kept continuous across periodic seams, poles and apices are spread
along their parameter row, and loops that go round a closed surface without a seam are closed
through the pole or joined into a band. The region is triangulated with constraints on the
boundary and refined where the surface leaves the triangles by more than the deflection. The
result is checked: a solid whose mesh is not closed and consistently oriented is an error.

## Testing

- Fixtures generated with CadQuery / OpenCascade (`reference/make_fixtures.py`): primitives, holes,
  voids, seams, poles, apices, splines, fillets, assemblies, inch and mixed units.
- A corpus gate (`tests/corpus.rs`, ignored by default) compares every file with OpenCascade
  (`reference/occt_reference.py`): solid count, bounds within 0.1 % of the diagonal, and volume
  within 0.1 % of OpenCascade's mesh at the same tolerances or of its exact B-rep volume (allowing
  the deflection times the area). Any result outside that is a failure, as is a panic or timeout.
- On the CAx-IF test files shipped with occt-import-js 0.0.23, of the 44 files OpenCascade reads,
  43 match; the remaining one is a file OpenCascade misreads (its rim plane puts the material on
  the other side; the same part from another exporter reads as stepmesh reads it).

## License

MIT OR Apache-2.0
