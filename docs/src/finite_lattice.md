# Finite lattices

A [`FiniteLattice`](@ref) is a finite cluster cut out of a [`Lattice`](@ref).
The cluster is spanned by boundary vectors ``\mathbf{t}_i`` that are integer
combinations of the lattice vectors. Each boundary direction can be periodic
or open.

```@example flattice
using Latlib

# 6x4 cylinder of the square lattice: open along a1, periodic along a2
fl = FiniteLattice(square, [6 0; 0 4], [false, true])
```

The boundary vectors can also be given as [`LatticeVector`](@ref)s, which is
convenient for tilted clusters:

```@example flattice
t1 = LatticeVector(honeycomb, [2, -2])
t2 = LatticeVector(honeycomb, [1, 1])
fl = FiniteLattice([t1, t2], true)
atoms(fl)
```

```@docs
FiniteLattice
atoms
bravais_cells
boundary
periodic_boundary
periodicity
```

## Site ordering

The order of the sites returned by [`atoms`](@ref) determines the site indices
used in [`Op`](@ref)s and in the TOML output. It is controlled by the
`bravais_order` and `atom_order` keyword arguments of [`FiniteLattice`](@ref).

```@example order
using Latlib

# sort sites by their x coordinate first, then by y
fl = FiniteLattice(shastry_sutherland, [3 0; 0 2], [false, true]; atom_order=order_xy)
atoms(fl)
```

```@docs
order_xy
order_yx
order_xyz
```

## Vectors in the basis of boundary vectors

```@docs
FiniteLatticeVector
to_finite_lattice_basis
```
