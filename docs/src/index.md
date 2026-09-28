# Latlib.jl

*Finite lattices and lattice models for quantum many-body simulations.*

Latlib.jl builds finite lattice clusters and lattice models. Starting from a
Bravais lattice with a basis of atoms, it

- generates finite clusters with periodic or open boundaries,
- finds neighbors under the periodic metric of the cluster,
- assembles interaction terms into an operator sum,
- plots clusters, bonds, and spin configurations,
- and writes site coordinates and interactions to a TOML file that can be
  consumed by exact diagonalization or other many-body codes such as
  [XDiag](https://github.com/awietek/xdiag).

## Installation

Latlib.jl is not yet registered. Install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/awietek/Latlib.jl")
```

Julia 1.10 or newer is required. Plotting uses
[GLMakie](https://docs.makie.org), which needs a working OpenGL display.

## Quick start

The typical workflow consists of four steps:

1. Choose a [`Lattice`](@ref), either predefined or custom.
2. Cut a [`FiniteLattice`](@ref) out of it by specifying boundary vectors.
3. Build an [`OpSum`](@ref) of interactions with [`neighbor_interaction`](@ref)
   and [`lattice_interaction`](@ref).
4. Write the result to a file with [`write_toml`](@ref) or plot it.

```@example quickstart
using Latlib

# 4x4 cluster of the triangular lattice, periodic in both directions
fl = FiniteLattice(triangular, [4 0; 0 4], true)
```

The sites of the cluster are enumerated by [`atoms`](@ref):

```@example quickstart
sites = atoms(fl)
length(sites)
```

A nearest neighbor Heisenberg model is generated with
[`neighbor_interaction`](@ref):

```@example quickstart
H = OpSum()
H += neighbor_interaction("SdotS", "J", fl; num_distance=1)
length(H.ops)
```

Finally, the cluster and its interactions are written to a TOML file
(here we print the file content instead):

```@example quickstart
toml = write_toml(fl, H, "triangular-N-16.toml"; zero_based=true, return_string=true)
println(join(split(toml, "\n")[1:20], "\n"))  # first 20 lines
```

To visualize the cluster interactively, call `plot_opsum(H, fl)`; see
[Plotting](@ref).

## Site ordering

Sites in a `FiniteLattice` are enumerated by [`atoms`](@ref). By default, all
copies of the first atom of the unit cell come first, ordered by their Bravais
coordinates, then all copies of the second atom, and so on. Both orderings can
be customized with the `bravais_order` and `atom_order` keyword arguments of
[`FiniteLattice`](@ref); the helpers [`order_xy`](@ref), [`order_yx`](@ref),
and [`order_xyz`](@ref) sort by Cartesian coordinates.

Site indices in [`Op`](@ref) and [`OpSum`](@ref) are 1-based. Use
`zero_based=true` in [`write_toml`](@ref) to produce 0-based indices for C++
or Python codes.

## Contents

```@contents
Pages = ["lattice.md", "finite_lattice.md", "metric.md", "opsum.md", "io.md", "plots.md"]
Depth = 2
```

## Index

```@index
```
