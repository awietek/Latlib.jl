# Plotting

Latlib uses [GLMakie](https://docs.makie.org) for interactive plots of finite
lattices. Plotting requires a working OpenGL display.

## 2D lattices

[`plot`](@ref) draws the sites, nearest neighbor bonds, and boundary box of a
two-dimensional finite lattice. [`plot_opsum`](@ref) additionally draws the
bonds of an [`OpSum`](@ref), colored by coupling:

```julia
using Latlib

fl = FiniteLattice(honeycomb, [2 -2; 1 1], true)
H = OpSum()
H += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0])
H += lattice_interaction("SySy", "KY", fl, 1, 2, [0, -1])
H += lattice_interaction("SzSz", "KZ", fl, 2, 1, [1, 0])

plot_opsum(H, fl)
```

```@docs
plot
plot_opsum
```

## 3D lattices

[`plot_3d`](@ref) draws three-dimensional finite lattices, optionally with
the bonds of an `OpSum` and with spin arrows at each site:

```julia
using Latlib

t = [LatticeVector(hyperhoneycomb, [-1, 1, 1]),
     LatticeVector(hyperhoneycomb, [1, 1, -1]),
     LatticeVector(hyperhoneycomb, [-1, 1, -1])]
fl = FiniteLattice(t, true)

H = neighbor_interaction("SdotS", "J", fl)
plot_3d(fl, H; cpl_dict=Dict("J" => :black), annotate_sites=true)
```

See `examples/spins_3d_plot` for plotting a spin configuration.

```@docs
plot_3d
```
