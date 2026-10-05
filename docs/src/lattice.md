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

| Name | Dimension | Atoms per unit cell |
|:-----|:---------:|:-------------------:|
| [`square`](@ref) | 2 | 1 |
| [`triangular`](@ref) | 2 | 1 |
| [`honeycomb`](@ref) | 2 | 2 |
| [`kagome`](@ref) | 2 | 3 |
| [`shastry_sutherland`](@ref) | 2 | 4 |
| [`shastry_sutherland_non_symmorphic`](@ref) | 2 | 4 |
| [`lieb`](@ref) | 2 | 3 |
| [`trellis`](@ref) | 2 | 2 |
| [`maple_leaf`](@ref) | 2 | 6 |
| [`hyperhoneycomb`](@ref) | 3 | 4 |
| [`simple_cubic`](@ref) | 3 | 1 |
| [`bcc`](@ref) | 3 | 1 |
| [`fcc`](@ref) | 3 | 1 |
| [`diamond`](@ref) | 3 | 2 |
| [`pyrochlore`](@ref) | 3 | 4 |
| [`simple_hexagonal`](@ref) | 3 | 1 |
| [`hcp`](@ref) | 3 | 2 |

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
