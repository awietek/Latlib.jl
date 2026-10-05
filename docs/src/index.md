# Latlib.jl

*Finite lattices and lattice models for quantum many-body simulations.*

Latlib.jl builds finite lattice clusters and lattice models. Starting from a
Bravais lattice with a basis of atoms, it

- generates finite clusters with periodic or open boundaries,
- finds neighbors under the periodic metric of the cluster,
- assembles interaction terms into an operator sum,
- determines the symmetry group of a cluster and the irreducible representations of its space
  group (momenta, little groups, characters) for symmetry-resolved exact diagonalization, see
  [Symmetries](symmetries.md),
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

The cluster and its interactions are written to a TOML file with
[`write_toml`](@ref). Site indices in the file are 1-based by default; pass
`zero_based=true` for C++ or Python codes. Here we print the file content
instead of writing it:

```@example quickstart
toml = write_toml(fl, H, "triangular-N-16.toml"; zero_based=true, return_string=true)
println(join(split(toml, "\n")[1:20], "\n"))  # first 20 lines
```

Finally, [`plot_opsum`](@ref) opens an interactive window showing the sites
of the cluster, its boundary box, and the bonds of the `OpSum`, colored by
coupling. In the REPL, the window stays open. In a script, it would close as
soon as the script ends, so wait until the user closes it:

```julia
f, ax = plot_opsum(H, fl)
wait(display(f))   # keeps the window open; not needed in the REPL
```

To draw into an existing Makie axis instead, e.g. to save the figure to a
file, pass it as the `ax` keyword argument:

```@example quickstart
using GLMakie: Figure, Axis, DataAspect, hidedecorations!, save

f = Figure(size=(600, 400))
ax = Axis(f[1, 1], aspect=DataAspect())
plot_opsum(H, fl; ax=ax)
hidedecorations!(ax)
save("triangular-N-16.png", f)
f
```

Three-dimensional lattices are drawn with [`plot_3d`](@ref); see
[Plotting](@ref). More complete models, including Kitaev interactions on the
honeycomb and hyperhoneycomb lattices, are shown on the [Examples](@ref) page.

## Contents

```@contents
Pages = ["lattice.md", "finite_lattice.md", "mps_ordering.md", "metric.md", "opsum.md", "io.md", "plots.md", "examples.md"]
Depth = 2
```

## Index

```@index
```
