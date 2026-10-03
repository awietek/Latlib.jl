# Symmetries

[`spacegroup`](@ref) determines the symmetry group of a [`Lattice`](@ref) or of a
periodic [`FiniteLattice`](@ref) using [spglib](https://spglib.readthedocs.io).
Atom positions and atom `types` are taken into account; couplings between the atoms
are not. Both symmorphic and non-symmorphic groups are supported.

## Infinite lattices

For an infinite lattice, [`spacegroup`](@ref) returns a [`SpaceGroup`](@ref). It holds
the Hermann–Mauguin symbol and ITA number of the group, its point group, and one
[`SymmetryOperation`](@ref) per coset of the Bravais translations. Two-dimensional
lattices are treated as a planar layer, so the plane group (wallpaper group) is found.

```@example symmetries
using Latlib
sg = spacegroup(honeycomb)
```

For symmorphic groups, `sg.origin` is a point (in the lattice basis) whose site-symmetry
group is the full point group. For the honeycomb lattice this is the center of a hexagon:

```@example symmetries
sg.origin
```

A symmetry operation ``\mathcal{X} \mapsto W\mathcal{X} + \mathbf{w}`` acts on
coordinates in the lattice basis. [`cartesian_rotation`](@ref) returns its rotational
part in Cartesian coordinates:

```@example symmetries
op = operations(sg)[2]
cartesian_rotation(op, honeycomb)
```

Non-symmorphic groups such as that of the diamond lattice are identified as well:

```@example symmetries
spacegroup(diamond).symbol
```

## Finite lattices

On a finite lattice with periodic boundaries, only those operations survive whose
rotational part maps the torus onto itself. [`spacegroup`](@ref) returns a
[`FiniteSpaceGroup`](@ref) with all operations of the cluster (rotational parts times
all Bravais translations of the cluster), the corresponding site permutations, and the
type of the cluster's symmetry group:

```@example symmetries
fl = FiniteLattice(honeycomb, [2 -2; 1 1], true)   # rectangular torus
g = spacegroup(fl)
```

The rectangular torus breaks the sixfold rotation of the honeycomb lattice, so the
cluster only has the symmetry of the plane group c2mm.

[`site_permutations`](@ref) gives the action of every operation on the sites:
`site_permutations(g)[k][i]` is the index (in [`atoms`](@ref)`(fl)`) of the image of
site `i` under the `k`-th operation.

```@example symmetries
site_permutations(g)[2]
```

Only fully periodic finite lattices are supported.

!!! note
    On very small clusters several operations can act identically on the sites, for
    example a two-site cluster has many more operations than distinct site permutations.

## API

```@docs
spacegroup
SpaceGroup
FiniteSpaceGroup
SymmetryOperation
operations
site_permutations
issymmorphic
pointgroup_operations
cartesian_rotation
```
