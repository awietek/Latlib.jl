# Triangular lattice: J1-J2 Heisenberg model on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/triangular/triangular.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The space group p6mm is symmorphic, so every TOML file contains the symmetry operations and the
# irreducible representations of the cluster (symmetry center: a site).

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = triangular
d1, d2 = 1.0, sqrt(3)   # distances of nearest and next-nearest neighbors

# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [4 0; 0 4]),   # N = 16: 4 × 4, point group 6mm
    (1, [-1 -4; 4 -5]),   # N = 21: tilted, 6 (chiral)
    (1, [3 3; -6 3]),   # N = 27: tilted, 6mm
    (1, [6 0; 0 6]),   # N = 36: 6 × 6, 6mm
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # J1-J2 Heisenberg model
    H = OpSum()
    H += neighbor_interaction("SdotS", "J1", fl; num_distance=1)
    H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)
    check_couplings(fl, H, Dict("J1" => (3, d1), "J2" => (3, d2)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "triangular", "J1J2", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(H, fl)
        show_figure(f)
    end
end
