# Plotting

Latlib uses [GLMakie](https://docs.makie.org) for interactive plots of finite
lattices. Plotting requires a working OpenGL display.

All plotting functions display the figure in an interactive window and return
the Makie figure and axis. In the REPL, the window stays open. When a plot is
made from a script, however, the window closes as soon as the script ends.
To keep it open until it is closed by the user, wait for its display:

```julia
f, ax = plot_opsum(H, fl)
wait(display(f))
```

To draw into your own Makie figure, e.g. to combine several plots or to save
a figure to a file, pass an axis with the `ax` keyword argument of
[`plot`](@ref) and [`plot_opsum`](@ref):

```julia
using GLMakie
f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
plot_opsum(H, fl; ax=ax)
save("cluster.png", f)
```

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

f, ax = plot_opsum(H, fl)
wait(display(f))
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
f, ax = plot_3d(fl, H; cpl_dict=Dict("J" => :black), annotate_sites=true)
wait(display(f))
```

See `examples/spins_3d_plot` for plotting a spin configuration.

```@docs
plot_3d
```
