# Symmetries

Latlib determines the symmetry group of a [`Lattice`](@ref) or of a periodic
[`FiniteLattice`](@ref) with [spglib](https://spglib.readthedocs.io), and the irreducible
representations of the space group of a cluster, which split its Hilbert space into symmetry
sectors for exact diagonalization.

## Quick start

```@example symmetries
using Latlib
fl = FiniteLattice(triangular, [3 0; 0 3], true)   # 9 sites with the full symmetry p6mm
H = neighbor_interaction("SdotS", "J", fl)
cs = symmetries(fl)
[irrep.label for irrep in irreps(cs)]
```

Every sector is labelled by a momentum, its little co-group and an irrep of the latter (see
[Irreducible representations](@ref)). Passing `symmetries=true` to [`write_toml`](@ref) adds the
site permutations of all symmetry operations and the characters of every sector to the TOML file
(see [Format of the TOML file](@ref)):

```julia
write_toml(fl, H, "triangular-9.toml"; zero_based=true, symmetries=true)
```

An exact diagonalization code builds the block of a sector from these entries, e.g.
[XDiag](https://github.com/awietek/xdiag):

```julia
using XDiag
file = FileToml("triangular-9.toml")
H = read_opsum(file, "Interactions")
H["J"] = 1.0
block = Spinhalf(9, read_representation(file, "Gamma.C6v.A1"))
e0 = eigval0(H, block)
```

!!! warning "Symmetries of the lattice, not of the Hamiltonian"
    The symmetries are determined from the atom positions and atom `types` of the lattice only;
    the couplings of the Hamiltonian are not taken into account. If a model has less symmetry
    than its lattice, e.g. because of bond-dependent couplings as in Kitaev models, or because
    bonds of the same length carry different couplings, the TOML file contains operations that
    are not symmetries of the Hamiltonian. An exact diagonalization in these sectors then gives
    wrong results without any error. A check of the Hamiltonian is not implemented yet, so make
    sure that the lattice has the symmetry of the model:

    - Inequivalent sites can be distinguished by atom `types`, e.g. the two sublattices of the
      honeycomb lattice for a staggered field (below).
    - The sites of [`shastry_sutherland`](@ref) form a square lattice (p4mm). The dimer bonds of
      the Shastry–Sutherland model, which lower the symmetry to p4gm, are invisible to
      [`spacegroup`](@ref). Use [`shastry_sutherland_non_symmorphic`](@ref), whose geometry has the
      symmetry of the model (only the site permutations are written for it, see
      [Limitations](@ref)).

```@example symmetries
staggered = Lattice(honeycomb.A, honeycomb.positions; types=[1, 2])   # sublattices A and B
spacegroup(honeycomb).symbol, spacegroup(staggered).symbol
```

## Infinite lattices

For an infinite lattice, [`spacegroup`](@ref) returns a [`SpaceGroup`](@ref). It holds
the Hermann–Mauguin symbol and ITA number of the group, its point group, and one
[`SymmetryOperation`](@ref) per coset of the Bravais translations. Two-dimensional
lattices are treated as a planar layer, so the plane group (wallpaper group) is found. The groups
of the predefined lattices are listed in [Predefined lattices](@ref).

```@example symmetries
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

[`symmetries`](@ref) prepares the symmetry operations of a periodic finite lattice,
[`momenta`](@ref) lists the momenta it resolves, and [`irreps`](@ref) returns the sectors with
their characters. This is implemented for two- and three-dimensional lattices with a symmorphic
space group.

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

How the axes that the Mulliken names refer to are chosen, so that the names do not depend on how
the cluster is described, is explained in [Conventions in detail](@ref).

The character of an operation ``\mathcal{X} \mapsto W(\mathcal{X} - \mathbf{c}) + \mathbf{c} + \mathbf{t}``,
written relative to the symmetry center ``\mathbf{c}``, is ``\rho(W)\, e^{+i\mathbf{k}\cdot\mathbf{t}}``.

On small or thin clusters, irreps that are not trivial on operations acting trivially on the
sites vanish; they are left out.

### Symmetry center

The symmetry center is a point with the full point-group symmetry of the lattice. If it is
unique, it is chosen automatically, like the center of a hexagon of the maple-leaf lattice
above or the site of the body-centered cubic lattice. Otherwise it has to be given with the
keyword `origin`, e.g. for the square lattice (site or plaquette center) or the simple cubic,
face-centered cubic and simple hexagonal lattices:

```@example symmetries
cs = symmetries(FiniteLattice(square, [4 0; 0 4], true); origin=LatticeVector(square, [0.0, 0.0]))
[(k.label, k.littlegroup_name) for k in momenta(cs) if k.representative]
```

Without `origin`, the error message lists the candidates. As everywhere in Latlib, the type of
the vector determines its meaning: a [`LatticeVector`](@ref) refers to the basis of its lattice,
an [`EuclideanVector`](@ref) to Cartesian coordinates. For the trellis lattice, whose basis is not
orthogonal, the center at the Cartesian point ``(1/2, 1/2)`` is given most easily as an
`EuclideanVector`:

```@example symmetries
cs = symmetries(FiniteLattice(trellis, [4 0; -2 4], true); origin=EuclideanVector([0.5, 0.5]))
cs.origin   # in the lattice basis
```

### Three-dimensional lattices

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

### Format of the TOML file

After the `Coordinates` and `Interactions` (see [Reading and writing files](io.md)), a file
written with `symmetries=true` contains:

1. **Comments** describing the symmetries: the symmetry center (lattice basis and Cartesian), the
   symmetry groups of the cluster and of the infinite lattice, and all momenta with their labels,
   little co-groups and Cartesian coordinates in the first Brillouin zone; the representatives of
   the stars are marked with `*`.
2. **`Symmetries`**: one site permutation per distinct symmetry operation, the identity first.
   `Symmetries[k][i]` is the site that site `i` is mapped to by operation `k`, as in
   [`site_permutations`](@ref). Site indices follow `zero_based`.
3. **One table per sector**, named by its label. In TOML, the header `[K.C3.Ea]` defines the
   nested tables `K`, `K.C3` and `K.C3.Ea`. Each contains
   - `allowed_symmetries`: the indices (into `Symmetries`, following `zero_based`) of the
     operations forming the group of the sector, i.e. the little group of the momentum, or the
     subgroup ``\ker(\det E)`` for the partners of a two-dimensional irrep ``E``;
   - `characters`: the character of each allowed symmetry as a `[real, imag]` pair, in the order
     of `allowed_symmetries`;
   - `momentum`: the Cartesian momentum in the first Brillouin zone.

```toml
Symmetries = [
  [0, 1, 2, 3, 4, 5, 6, 7, 8],
  [1, 2, 0, 4, 5, 3, 7, 8, 6],
  ...
]

# Irreducible representations
[Gamma.C6v.A1]
characters = [
  [1.0000000000000000, 0.0000000000000000],
  ...
]
allowed_symmetries = [0, 1, 2, ...]
momentum = [0.0000000000000000, 0.0000000000000000]
...
```

Further comments appear where needed:

- a comment before the partners of a two-dimensional irrep, e.g.
  `# Gamma.C6v.E1 is two-dimensional: its partners E1a and E1b ...`;
- a banner `# !!! WARNING: three-dimensional irreducible representations were SKIPPED ...` at the top
  of the file and at the start of the irreps, listing the skipped irreps;
- on clusters with operations acting trivially on the sites, blocks marked
  `# PLACEHOLDER (TOML format to be decided)` listing the omitted operations and the irreps that
  vanish because of them. The format of these blocks is provisional.

## Limitations

- Only fully periodic finite lattices are supported.
- The couplings of the Hamiltonian are not taken into account, see the warning at the top of this
  page.
- Irreducible representations are implemented for symmorphic space groups only. For
  non-symmorphic ones (e.g. [`shastry_sutherland_non_symmorphic`](@ref),
  [`hyperhoneycomb`](@ref), [`diamond`](@ref), [`pyrochlore`](@ref), [`hcp`](@ref)), only the site
  permutations are written.
- Three-dimensional irreps (`T` irreps of the little co-groups T, Th, O, Td and Oh) are skipped, so
  on cubic clusters the written sectors do not span the full Hilbert space.
- Time reversal is not used: momenta ``\mathbf{k}`` and ``-\mathbf{k}`` belong to the same star only
  if a symmetry of the cluster relates them. Otherwise both are listed, and their sectors are
  degenerate.
- If the unit cell of the lattice is not primitive, only translations by lattice vectors are used
  (with a warning), and the momenta refer to the Brillouin zone of the given unit cell.
- Operations acting trivially on the sites are written only once, and irreps that vanish because
  of them are left out; the format of the comments reporting this is provisional.
- [`sector_dimension`](@ref) counts the states of spin-1/2 only.

## Conventions in detail

### Choice of axes

The axes that the Mulliken names refer to are chosen from the momentum, the lattice and the
cluster, so that every momentum of a star gets the same names (carried over by the symmetries of
the cluster), and the names do not depend on how the cluster is described: on the Cartesian
orientation, the basis of the lattice, the boundary vectors or the order of the sites. If the
cluster has less symmetry than the lattice, it distinguishes between equivalent axes of the
lattice (e.g. the two axes of the square lattice on a ``4 \times 2`` torus). The conventional
basis ``\mathbf{a}, \mathbf{b}, (\mathbf{c})`` is then the one in which the cluster has the
lexicographically smallest description: the Hermite normal form of the lattice of torus vectors,
then the positions of the atoms relative to the symmetry center. Where the little co-group alone
does not fix the axes:

- `C2v` (`B1` even under ``\sigma(xz)``): the mirror containing ``\mathbf{k}``. If the twofold axis is
  along ``\mathbf{k}``, the mirror perpendicular to the principal axis of the lattice (for cubic
  lattices: the mirror whose normal is a cubic axis). For a twofold axis along the principal axis,
  the mirror whose normal is a lattice vector (hexagonal lattices) or the mirror containing ``\mathbf{a}``.
- `D2`, `D2h` (`B1`, `B2`, `B3` even under the twofold rotations about ``z``, ``y``, ``x``): ``z``
  along the principal axis of the lattice, ``x`` along ``\mathbf{k}`` or in the plane of
  ``\mathbf{k}`` and ``z``.

### Representatives and numbering of the momenta

The representative of a star is the momentum whose image in the first Brillouin zone has the
largest Cartesian coordinates ``(k_x, k_y, k_z)``, compared in this order, as in earlier files
for exact diagonalization. Its `momentum` is the one written to the TOML file. The labels do not
depend on this choice. Repeated momentum labels (`Sigma0`, `Sigma1`, …, `GP0`, `GP1`, …) are
numbered in an order fixed by the conventional basis described above, so that the numbering does
not depend on the orientation of the lattice either.

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
