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

- Computes exact **external bounds** on edges, faces, and flatness, and
  **internal bounds** on space-filling tetrahedra (`T`), internal triangulation segments (`S_i`),
  and internal gluing triangles (`N_i`), from a given vertex count
  (`V_report()`).

- Enumerates **all admissible face-type distributions** for fixed vertex count
  `V` and flatness `S`, purely combinatorially (`Pk_wizard()`).

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
assess("4x3+6x4")

# External and internal symbolic bounds from vertex count
V_report(8)

# Enumerate all admissible face-type configurations for fixed V and S
Pk_wizard(V = 8, S = 6)

# Enumerate all admissible face-type configurations for fixed V and S at certain k_max (face degrees)
Pk_wizard(V = 8, S = 6, k_max=5)

```
