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

## Irreducible representations

For exact diagonalization, the Hilbert space of a cluster is split into sectors labelled by
the irreducible representations (irreps) of its space group. [`symmetries`](@ref) prepares the
symmetry operations of a periodic finite lattice, [`momenta`](@ref) lists the momenta it
resolves, and [`irreps`](@ref) returns the sectors with their characters. This is implemented
for two- and three-dimensional lattices with a symmorphic space group.

```@example symmetries
cs = symmetries(FiniteLattice(maple_leaf, [1 1; 1 -2], true))
[irrep.label for irrep in irreps(cs)]
```

The irreps are constructed for one representative momentum of each star, i.e. of each orbit of
momenta under the point group of the cluster. Two momenta that are related by a symmetry of
the infinite lattice but not of the cluster belong to different stars. A sector is labelled
`"<momentum>.<little co-group>.<irrep>"`:

- **Momenta** are labelled following the Bilbao Crystallographic Server (CDML notation) for the
  Bravais lattice, e.g. `Gamma`, `K`, `M` for the hexagonal lattice and `Sigma`, `Lambda`, `T` on
  the lines between them, or `Gamma`, `H`, `N`, `P` for the body-centered cubic lattice. Generic
  momenta are labelled `GP0`, `GP1`, …; other labels are numbered only if they occur more than once.
- **Little co-groups** are named by their Schoenflies symbol, e.g. `C2v`, `C6v`, `D4h`, `Oh`.
- **Irreps** carry their Mulliken symbol, with the orientation conventions of the character
  tables of the Bilbao Crystallographic Server and `p`/`pp` for primes (`Ap`, `App`, `A1pp`).
  The complex conjugate pairs of one-dimensional irreps are labelled `a`/`b` (e.g. `E1a`, `E1b`,
  `Ega`), where `a` belongs to ``\exp(+2\pi i m/n)`` on the rotation by ``+2\pi/n`` about the
  principal axis.
- **Two-dimensional irreps** (e.g. `E1` of `C6v`, `Eg` of `D4h`) are represented by two exactly
  degenerate partners `a`/`b`, one-dimensional representations of the subgroup
  ``\ker(\det E)`` (e.g. `C6`, `C4h`).
- **Three-dimensional irreps** (`T` irreps of the cubic little co-groups) cannot be written as
  one-dimensional characters. They are skipped with a warning; TOML files get a prominent banner.

The axes that the Mulliken names refer to are chosen from the momentum, the lattice and the
cluster, so that every momentum of a star gets the same names (carried over by the symmetries of
the cluster), and the names do not depend on how the cluster is described: on the Cartesian
orientation, the basis of the lattice, the boundary vectors or the order of the sites. If the
cluster has less symmetry than the lattice, it distinguishes between equivalent axes of the
lattice (e.g. the two axes of the square lattice on a ``4 \times 2`` torus). The conventional
basis ``\mathbf{a}, \mathbf{b}, (\mathbf{c})`` is then the one in which the cluster has the
lexicographically smallest description: the Hermite normal form of the lattice of torus vectors,
then the positions of the atoms relative to the symmetry center. Repeated momentum labels are
numbered in an order fixed by the same basis. Where the little co-group alone does not fix the axes:

- `C2v` (`B1` even under ``\sigma(xz)``): the mirror containing ``\mathbf{k}``. If the twofold axis is
  along ``\mathbf{k}``, the mirror perpendicular to the principal axis of the lattice (for cubic
  lattices: the mirror whose normal is a cubic axis). For a twofold axis along the principal axis,
  the mirror whose normal is a lattice vector (hexagonal lattices) or the mirror containing ``\mathbf{a}``.
- `D2`, `D2h` (`B1`, `B2`, `B3` even under the twofold rotations about ``z``, ``y``, ``x``): ``z``
  along the principal axis of the lattice, ``x`` along ``\mathbf{k}`` or in the plane of
  ``\mathbf{k}`` and ``z``.

The character of an operation ``\mathcal{X} \mapsto W(\mathcal{X} - \mathbf{c}) + \mathbf{c} + \mathbf{t}``,
written relative to the symmetry center ``\mathbf{c}``, is ``\rho(W)\, e^{+i\mathbf{k}\cdot\mathbf{t}}``.

The symmetry center is a point with the full point-group symmetry of the lattice. If it is
unique, it is chosen automatically, like the center of a hexagon of the maple-leaf lattice
above or the site of the body-centered cubic lattice. Otherwise it has to be given, e.g. for the
square lattice (site or plaquette center) or the simple cubic, face-centered cubic and simple
hexagonal lattices:

```@example symmetries
cs = symmetries(FiniteLattice(square, [4 0; 0 4], true); origin=LatticeVector(square, [0.0, 0.0]))
[(k.label, k.littlegroup_name) for k in momenta(cs) if k.representative]
```

On small or thin clusters, irreps that are not trivial on operations acting trivially on the
sites vanish; they are left out.

Three-dimensional lattices work the same way:

```@example symmetries
cs = symmetries(FiniteLattice(bcc, 2 * [0 1 1; 1 0 1; 1 1 0], true))   # 16 sites
[(k.label, k.littlegroup_name) for k in momenta(cs) if k.representative]
```

## Dimensions of the sectors

[`sector_dimension`](@ref) gives the dimension of a sector in the Hilbert space of spin-1/2, i.e. the
size of the block an exact diagonalization code works with, optionally for a fixed number of up
spins and a parity under the global spin flip. It is computed exactly from the cycles of the site
permutations, without enumerating any states, and takes milliseconds also for large clusters. For
example, the ground-state sector of the Heisenberg model on the ``6 \times 6`` square lattice:

```@example symmetries
cs = symmetries(FiniteLattice(square, [6 0; 0 6], true); origin=LatticeVector(square, [0.0, 0.0]))
A1 = only(irrep for irrep in irreps(cs) if irrep.label == "Gamma.C4v.A1")
sector_dimension(cs, A1; nup=18, spinflip=1)
```

Whenever the irreps are computed, also by [`write_toml`](@ref), the sectors are checked against
the spin-1/2 Hilbert space. The dimension of every sector has to be an integer, which fails for
characters that do not form a representation. At every momentum, the sectors have to span the
subspace of this momentum, apart from skipped three-dimensional irreps, so that all sectors
together, weighted with the sizes of the stars, span the ``2^N`` states. A failing check raises an
error.

## Writing symmetries to TOML files

Pass `symmetries=true` to [`write_toml`](@ref) to append the site permutations of the cluster's
symmetry operations as a `Symmetries` section (see [`toml_symmetries`](@ref)). This works for
every lattice:

- **Symmorphic space groups**: the irreducible representations are written as well (see
  [`toml_irreps`](@ref)), unless `irreps=false`. If the symmetry center is ambiguous, it is passed
  with the keyword `origin`. Skipped three-dimensional irreps are listed in a banner at the top
  of the file.
- **Non-symmorphic space groups**: only the symmetry operations are written, and a warning states
  that irreducible representations are not implemented for them yet.

Instead of `true`, a precomputed [`FiniteSpaceGroup`](@ref) or [`ClusterSymmetries`](@ref) can be
passed.

## API

```@docs
symmetries
ClusterSymmetries
momenta
ClusterMomentum
irreps
Irrep
sector_dimension
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
