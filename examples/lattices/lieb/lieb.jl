# Lieb lattice: J1-J2 Heisenberg model on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/lieb/lieb.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# J1 couples the corner sites (four neighbors) to the edge sites (two neighbors), J2 couples
# neighboring edge sites. The space group p4mm is symmorphic, so every TOML file contains the
# symmetry operations and the irreducible representations of the cluster. Of the two inequivalent
# symmetry centers (corner site and center of a square), the corner site is used.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = lieb
d1, d2 = 0.5, sqrt(2) / 2   # distances of nearest and next-nearest neighbors

# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [2 0; 0 2]),   # N = 12: 2 × 2, point group 4mm
    (1, [2 2; -2 2]),   # N = 24: tilted, 4mm
    (1, [3 0; 0 3]),   # N = 27: 3 × 3, 4mm
    (1, [4 0; 0 4]),   # N = 48: 4 × 4, 4mm
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # J1-J2 Heisenberg model
    H = OpSum()
    H += neighbor_interaction("SdotS", "J1", fl; num_distance=1)
    H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)
    check_couplings(fl, H, Dict("J1" => (4, d1), "J2" => (4, d2)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "lieb", "J1J2", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, origin=LatticeVector(lattice, [0.0, 0.0]))
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(H, fl)
        show_figure(f)
    end
end
