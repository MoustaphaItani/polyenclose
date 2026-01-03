# polyenclose

**polyenclose** is an R package implementing a *symbolic, coordinate-free* method
for testing whether a given multiset of polygonal faces can enclose a closed,
genus-0 polyhedral surface, and for enumerating all face-type configurations
consistent with global combinatorial constraints.

The package implements the symbolic framework developed in the accompanying
preprint:

> *Symbolic Constraints in Polyhedral Enclosure and Tetrahedral Decomposition in Genus-0 Polyhedra*
Available on ArXiv via https://arxiv.org/abs/2508.18222

No geometric embedding or construction is attempted.

---

## What this package does

- Tests **symbolic enclosure feasibility** of polygon multisets using Euler,
  incidence, and flatness constraints (`assess()`) (false positives are possible).

- Uses the total flatness parameter, defined as the count of external / boundary
  triangulation segments (`S`), as a global obstruction variable.

- Computes exact **external bounds and ladders** on edges, faces, and flatness, and
  **internal bounds and ladders** on space-filling tetrahedra (`T`), internal triangulation segments (`S_i`),
  and internal gluing triangles (`N_i`), from a given vertex count
  (`V_report()`).

- Enumerates **all admissible face-type distributions** for fixed vertex count
  `V` and flatness `S`, purely combinatorially (`Pk_wizard()`) (false positives are possible).

- Summarizes the *solution space* of admissible face-type configurations
  (face-degree prevalence, co-occurrence patterns, and diversity measures)
  for fixed `V` and `S` (`Pk_summary()`).

All results represent **necessary combinatorial conditions** for enclosure.

---

## What this package does *not* do

- It does **not** construct geometric polyhedra.
- It does **not** test geometric realizability in ℝ³.
- It does **not** generate tetrahedral decompositions.
- It does **not** claim sufficiency of symbolic feasibility.

The package deliberately separates symbolic admissibility from geometry.

---

## Quick examples

```r

# Test whether a set of polygons can enclose a genus-0 surface (symbolically)
assess("6x4")          # cube
assess("12x5")         # dodecahedron
assess("4x3+5x4")      # capped-cube / house-like
assess("36x6")         # too flat to enclose (symbolically)

# External + internal symbolic bounds from vertex count
V_report(8)

# Internal decomposition ladder only (decomposition via MININAL INTERNAL EDGE (SALT+MIE))
decompose(8)

# Vertex ladder compatible with a fixed tetrahedron count T (SALT+MIE)
V_ladder_from_T(6)

# Enumerate admissible face-type configurations for fixed V and S
out <- Pk_wizard(V = 8, S = 3, show_solutions = FALSE)
Pk_summary(out)

# Restrict enumeration to face degrees ≤ k_max
out <- Pk_wizard(V = 8, S = 3, k_max = 5, show_solutions = FALSE)
Pk_summary(out)

---

## Interpreting enumeration results

For fixed vertex count `V` and flatness `S`, `Pk_wizard()` may admit many
distinct face-type configurations.  
`Pk_summary()` provides a compact summary of this solution space, reporting:

- how often each face degree appears across admissible configurations,
- common face-degree co-occurrence patterns (“footprints”),
- and a diversity measure describing how mixed face types are.

All quantities are **symbolic and combinatorial**, and summarize the set of
admissible face-type multisets rather than any single polyhedron.

```
