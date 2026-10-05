# Lattices

An infinite lattice is described by the [`Lattice`](@ref) type. It consists of
the Bravais lattice vectors, the atom positions inside the unit cell, and
optionally a type for each atom.

## Vectors

Latlib distinguishes two kinds of vectors:

- [`EuclideanVector`](@ref): a vector in Cartesian coordinates.
- [`LatticeVector`](@ref): a vector expressed in the basis of the Bravais
  lattice vectors of a given `Lattice`.

Conversions between the two are performed with [`to_lattice_basis`](@ref) and
[`to_euclidean_basis`](@ref).

```@example vectors
using Latlib

v = LatticeVector(honeycomb, [1, 1])
to_euclidean_basis(v)
```

```@docs
EuclideanVector
LatticeVector
to_lattice_basis
to_euclidean_basis
in_lattice
```

## Defining lattices

```@example lattice
using Latlib

# Kagome lattice: three atoms per unit cell, positions in the lattice basis
A = [1.0 0.0;
     0.5 sqrt(3)/2]
positions = [0.0 0.0;
             0.5 0.0;
             0.0 0.5]
kag = Lattice(A, positions)
```

```@docs
Lattice
dim
natoms
positions
get_position
lattice_vecs
```

## Predefined lattices

The following lattices are available as constants:

| Name | Dimension | Atoms per unit cell | Space group | Symmorphic |
|:-----|:---------:|:-------------------:|:------------|:----------:|
| [`square`](@ref) | 2 | 1 | p4mm (#11) | yes |
| [`triangular`](@ref) | 2 | 1 | p6mm (#17) | yes |
| [`honeycomb`](@ref) | 2 | 2 | p6mm (#17) | yes |
| [`kagome`](@ref) | 2 | 3 | p6mm (#17) | yes |
| [`shastry_sutherland`](@ref) | 2 | 4 | p4mm (#11)¹ | yes |
| [`shastry_sutherland_non_symmorphic`](@ref) | 2 | 4 | p4gm (#12) | no |
| [`lieb`](@ref) | 2 | 3 | p4mm (#11) | yes |
| [`trellis`](@ref) | 2 | 2 | c2mm (#9) | yes |
| [`maple_leaf`](@ref) | 2 | 6 | p6 (#16) | yes |
| [`hyperhoneycomb`](@ref) | 3 | 4 | Fddd (#70) | no |
| [`simple_cubic`](@ref) | 3 | 1 | Pm-3m (#221) | yes |
| [`bcc`](@ref) | 3 | 1 | Im-3m (#229) | yes |
| [`fcc`](@ref) | 3 | 1 | Fm-3m (#225) | yes |
| [`diamond`](@ref) | 3 | 2 | Fd-3m (#227) | no |
| [`pyrochlore`](@ref) | 3 | 4 | Fd-3m (#227) | no |
| [`simple_hexagonal`](@ref) | 3 | 1 | P6/mmm (#191) | yes |
| [`hcp`](@ref) | 3 | 2 | P6_3/mmc (#194) | no |

The space groups (plane groups in two dimensions) are determined by [`spacegroup`](@ref) from the
atom positions and types; couplings are not taken into account (see [Symmetries](symmetries.md)).
Irreducible representations are available for the symmorphic ones.

¹ The sites of `shastry_sutherland` form a square lattice (its unit cell of four sites is not
primitive). The dimer bonds of the Shastry–Sutherland model, which lower the symmetry to p4gm, are
not part of the lattice. In `shastry_sutherland_non_symmorphic`, the p4gm symmetry is built into
the geometry.

```@example predefined
using Latlib
honeycomb
```

```@docs
square
triangular
honeycomb
kagome
shastry_sutherland
shastry_sutherland_non_symmorphic
lieb
trellis
maple_leaf
hyperhoneycomb
simple_cubic
bcc
fcc
diamond
pyrochlore
simple_hexagonal
hcp
```
