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

## Operations acting trivially on the sites

On small or thin clusters, operations other than the identity can fix every site. For
example, on a torus that is only two unit cells long along ``y``, the mirror
``y \mapsto -y`` maps every site onto itself. [`trivial_operations`](@ref) returns these
operations; they form a normal subgroup, and each site permutation is then realized by
several operations:

```@example symmetries
g = spacegroup(FiniteLattice(square, [4 0; 0 2], true))
operations(g)[trivial_operations(g)]
```

Irreducible representations that are not trivial on these operations vanish on such a
cluster. [`distinct_operations`](@ref) selects one operation per distinct site permutation;
only these are written to TOML files by [`toml_symmetries`](@ref).

## Writing symmetries to TOML files

Pass `symmetries=true` to [`write_toml`](@ref) to append the site permutations of the
cluster's symmetry operations as a `Symmetries` section, see [`toml_symmetries`](@ref).

## API

```@docs
spacegroup
SpaceGroup
FiniteSpaceGroup
SymmetryOperation
operations
site_permutations
trivial_operations
distinct_operations
issymmorphic
pointgroup_operations
cartesian_rotation
```
