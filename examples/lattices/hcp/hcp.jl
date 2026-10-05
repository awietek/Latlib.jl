# Hexagonal close-packed lattice: J1-J2 Heisenberg model on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/hcp/hcp.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The space group P6_3/mmc is non-symmorphic, so the TOML files contain the symmetry operations of
# the cluster but no irreducible representations (`irreps=false`, as they are not implemented
# for non-symmorphic space groups yet).
# With the ideal ratio c/a = sqrt(8/3), every site has twelve neighbors at distance 1.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = hcp
d1, d2 = 1.0, sqrt(2)   # distances of nearest and next-nearest neighbors

# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [3 0 0; 0 3 0; 0 0 2]),   # N = 36: 3 × 3 × 2, point group 6/mmm
    (1, [2 4 0; -4 -2 0; 0 0 2]),   # N = 48: tilted in the plane × 2, 6/mmm
    (1, [3 0 0; 0 3 0; 0 0 3]),   # N = 54: 3 × 3 × 3, 6/mmm
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # J1-J2 Heisenberg model
    H = OpSum()
    H += neighbor_interaction("SdotS", "J1", fl; num_distance=1)
    H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)
    check_couplings(fl, H, Dict("J1" => (12, d1), "J2" => (6, d2)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "hcp", "J1J2", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, irreps=false)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_3d(fl, H; cpl_dict=Dict("J1" => :black, "J2" => :red), annotate_sites=true, annotate_sites_zero_based=true)
        show_figure(f)
    end
end
