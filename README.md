# Latlib.jl

[![Build Status](https://github.com/awietek/Latlib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/awietek/Latlib.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/awietek/Latlib.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://awietek.github.io/Latlib.jl)

**Latlib.jl** builds finite lattice clusters and lattice models for quantum
many-body simulations. Starting from a Bravais lattice with a basis, it
generates finite clusters with periodic or open boundaries, finds neighbors
under the periodic metric, assembles interaction terms into an operator sum,
plots the result, and writes everything to a TOML file that can be consumed by
exact diagonalization or other many-body codes such as
[XDiag](https://github.com/awietek/xdiag).

📖 **[Documentation](https://awietek.github.io/Latlib.jl)**

## Features

- **Lattices in 2D and 3D**: define a `Lattice` from its Bravais vectors and
  atom positions, or use one of the predefined lattices (`square`,
  `triangular`, `honeycomb`, `kagome`, `shastry_sutherland`,
  `shastry_sutherland_non_symmorphic`, `lieb`, `trellis`,
  `maple_leaf`, `hyperhoneycomb`, `simple_cubic`, `bcc`, `fcc`, `diamond`,
  `pyrochlore`, `simple_hexagonal`, `hcp`).
- **Finite clusters**: a `FiniteLattice` is a lattice cut by a boundary box
  given as integer combinations of the lattice vectors, with periodic or open
  boundaries in each direction.
- **Periodic distances**: compute distances, distance matrices, and k-th
  nearest neighbors either in plain Euclidean space or modulo the periodic
  boundaries of a cluster.
- **Symmetries**: determine the plane group or space group of a lattice and the
  symmetry group of a periodic cluster, including site permutations, with
  [spglib](https://spglib.readthedocs.io) (`spacegroup`). Symmorphic and
  non-symmorphic groups are supported.
- **Irreducible representations**: momenta, little groups and characters of the
  irreducible representations of the space group of a periodic cluster, labelled
  in the conventions of the Bilbao Crystallographic Server (`symmetries`, `irreps`),
  for symmorphic space groups in 2D and 3D, written to TOML files for exact
  diagonalization codes like [XDiag](https://github.com/awietek/xdiag).
- **Interactions**: build an `OpSum` of two-body operators from nearest
  neighbor rules (`neighbor_interaction`) or from explicit bonds between atoms
  in the unit cell (`lattice_interaction`), e.g. for Kitaev-type models.
- **TOML I/O**: write coordinates and interactions to a TOML file
  (`write_toml`) and read interactions back (`read_toml_interaction`).
- **Plotting**: interactive 2D and 3D plots of clusters, bonds, and spin
  configurations with [GLMakie](https://docs.makie.org).

## Installation

Latlib.jl is not yet registered. Install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/awietek/Latlib.jl")
```

Julia 1.9 or newer is required. Plotting uses GLMakie, which needs a working
OpenGL display.

## Quick start

Build a 4×4 periodic triangular cluster, define a nearest neighbor Heisenberg
model, and write it to a TOML file:

```julia
using Latlib

# 4x4 cluster of the triangular lattice, periodic in both directions
fl = FiniteLattice(triangular, [4 0; 0 4], true)

# Heisenberg interaction between nearest neighbors
H = OpSum()
H += neighbor_interaction("SdotS", "J", fl; num_distance=1)

# write site coordinates and interactions to a TOML file
write_toml(fl, H, "triangular-N-16.toml"; zero_based=true)

# interactive plot of the cluster with its bonds
f, ax = plot_opsum(H, fl)
wait(display(f))   # keeps the window open when run as a script
```

The written TOML file contains the coordinates of all sites and one entry per
interaction:

```toml
Coordinates = [
  [0.0, 0.0],
  [0.5, 0.8660254],
  ...
]

Interactions = [
  ['J', 'SdotS', 0, 1],
  ['J', 'SdotS', 0, 3],
  ...
]
```

### Custom lattices and bond-resolved interactions

Lattices are defined by a matrix whose rows are the Bravais vectors and a
matrix whose rows are the atom positions in the lattice basis. Bonds between
specific atoms of the unit cell can be added with `lattice_interaction`, which
repeats the bond over all Bravais cells of the cluster:

```julia
using Latlib

# honeycomb lattice: two atoms per unit cell
A = [cos(pi/6)  sin(pi/6);
     cos(pi/6) -sin(pi/6)]
pos = [0.0 0.0;
       1/3 1/3]
hc = Lattice(A, pos)

# cluster spanned by the boundary vectors t1 = 2 a1 - 2 a2 and t2 = a1 + a2
fl = FiniteLattice(hc, [2 -2; 1 1], true)

# Kitaev model: one bond type per direction
H = OpSum()
H += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0])   # atom 1 to atom 2 in the same cell
H += lattice_interaction("SySy", "KY", fl, 1, 2, [0, -1])  # atom 1 to atom 2 in cell [0, -1]
H += lattice_interaction("SzSz", "KZ", fl, 2, 1, [1, 0])   # atom 2 to atom 1 in cell [1, 0]
```

### 3D lattices

The same workflow applies in three dimensions. Boundary vectors can also be
given as `LatticeVector`s:

```julia
using Latlib

t = [LatticeVector(hyperhoneycomb, [-1, 1, 1]),
     LatticeVector(hyperhoneycomb, [1, 1, -1]),
     LatticeVector(hyperhoneycomb, [-1, 1, -1])]
fl = FiniteLattice(t, true)          # 16 sites, fully periodic

H = neighbor_interaction("SdotS", "J", fl)
f, ax = plot_3d(fl, H; cpl_dict=Dict("J" => :black), annotate_sites=true)
wait(display(f))
```

More complete scripts, including Kitaev interactions on the hyperhoneycomb
lattice, a Shastry-Sutherland cylinder, and plotting of spin configurations,
are in the [`examples`](examples) directory and are walked through on the
[Examples](https://awietek.github.io/Latlib.jl/examples/) page of the documentation.

## Site ordering

Sites in a `FiniteLattice` are enumerated by `atoms(fl)`. By default, all
copies of the first atom of the unit cell come first, ordered by their Bravais
coordinates, then all copies of the second atom, and so on. Both orderings can
be customized with the `bravais_order` and `atom_order` keyword arguments of
`FiniteLattice`; the helpers `order_xy`, `order_yx`, and `order_xyz` sort by
Cartesian coordinates. The documentation contains a worked example of
[choosing the MPS path on square, triangular, and kagome cylinders](https://awietek.github.io/Latlib.jl/mps_ordering/).

Site indices in `Op` and `OpSum` are 1-based. Use `zero_based=true` in
`write_toml` to produce 0-based indices for C++ or Python codes.

## Running the tests

```julia
using Pkg
Pkg.test("Latlib")
```

## License

Latlib.jl is released under the MIT License. See [LICENSE](LICENSE).
