# Pyrochlore lattice: J1-J2 Heisenberg model on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/pyrochlore/pyrochlore.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The space group Fd-3m is non-symmorphic, so the TOML files contain the symmetry operations of
# the cluster but no irreducible representations (`irreps=false`, as they are not implemented
# for non-symmorphic space groups yet).

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = pyrochlore
d1, d2 = sqrt(2) / 4, sqrt(6) / 4   # distances of nearest and next-nearest neighbors

# The conventional cubic cell (16 sites) is too small for J2.
# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [2 0 0; 0 2 0; 0 0 2]),   # N = 32: 2 × 2 × 2 primitive cells, point group m-3m
    (1, [3 -1 -1; -1 3 -1; -1 -1 3]),   # N = 64: torus vectors (-1, 1, 1), (1, -1, 1), (1, 1, -1) (Cartesian), m-3m
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # J1-J2 Heisenberg model
    H = OpSum()
    H += neighbor_interaction("SdotS", "J1", fl; num_distance=1)
    H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)
    check_couplings(fl, H, Dict("J1" => (12, d1), "J2" => (24, d2)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "pyrochlore", "J1J2", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, irreps=false)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_3d(fl, H; cpl_dict=Dict("J1" => :black, "J2" => :red), annotate_sites=true, annotate_sites_zero_based=true)
        show_figure(f)
    end
end
