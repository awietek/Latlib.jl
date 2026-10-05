# Trellis lattice: trellis model (coupled two-leg ladders) on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/trellis/trellis.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The trellis lattice consists of two-leg ladders (legs along a1, atoms 1 and 2), coupled by zigzag
# bonds between neighboring ladders; all bonds have length 1. The model has the couplings Jleg,
# Jrung and Jinter. The space group c2mm is symmorphic, so every TOML file contains the symmetry
# operations and the irreducible representations of the cluster. Of the two inequivalent symmetry
# centers (center of a rung and center of a plaquette of a ladder), the center of a rung is used.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = trellis
d = 1.0   # length of all bonds

# periodic clusters of L rungs × W ladders: (version, boundary vectors as rows in the lattice basis);
# the boundary vectors (L, 0) and (-W/2, W) span a rectangular torus
clusters = [
    (1, [4 0; -1 2]),   # N = 16: L = 4, W = 2, point group 2mm
    (1, [6 0; -1 2]),   # N = 24: L = 6, W = 2, 2mm
    (1, [4 0; -2 4]),   # N = 32: L = 4, W = 4, 2mm
    (1, [6 0; -2 4]),   # N = 48: L = 6, W = 4, 2mm
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # trellis model: legs, rungs and zigzag bonds between neighboring ladders
    H = OpSum()
    H += lattice_interaction("SdotS", "Jleg", fl, 1, 1, [1, 0])     # lower leg
    H += lattice_interaction("SdotS", "Jleg", fl, 2, 2, [1, 0])     # upper leg
    H += lattice_interaction("SdotS", "Jrung", fl, 1, 2, [0, 0])    # rung
    H += lattice_interaction("SdotS", "Jinter", fl, 1, 2, [0, -1])  # to the upper leg of the ladder below
    H += lattice_interaction("SdotS", "Jinter", fl, 1, 2, [1, -1])
    check_couplings(fl, H, Dict("Jleg" => (2, d), "Jrung" => (1, d), "Jinter" => (2, d)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "trellis", "Trellis", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, origin=EuclideanVector([0.0, 0.5]))
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(H, fl)
        show_figure(f)
    end
end
