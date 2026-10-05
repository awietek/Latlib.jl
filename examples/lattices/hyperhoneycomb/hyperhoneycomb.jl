# Hyperhoneycomb lattice: Heisenberg and Kitaev-Heisenberg models on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/hyperhoneycomb/hyperhoneycomb.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# Every site has three nearest neighbors, as on the honeycomb lattice. The clusters are those of
# examples/Hyperhoneycomb/hyperhoneycomb.jl, with the same versions. The space group Fddd is
# non-symmorphic, so the TOML files contain the symmetry operations of the cluster but no irreducible
# representations (`irreps=false`, as they are not implemented for non-symmorphic space groups yet).
# For the Kitaev-Heisenberg model, these are the symmetries of the lattice: the translations are
# symmetries of the Hamiltonian, but the point-group operations exchange the bonds KX, KY and KZ, so
# they are symmetries only together with a corresponding transformation of the spin components.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = hyperhoneycomb
d = sqrt(2)   # distance of nearest neighbors

# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [-1 1 1; 0 1 -1; -1 1 -1]),    # N = 8
    (1, [-1 1 1; 1 1 -1; -1 1 -1]),    # N = 16
    (2, [-1 1 1; 0 1 -1; -2 2 -2]),    # N = 16
    (2, [-1 1 1; 1 1 -1; -2 2 -2]),    # N = 32: two copies of N = 16, version 1, in front of each other
    (3, [-2 2 2; 1 1 -1; -1 1 -1]),    # N = 32: two copies of N = 16, version 1, on top of each other
    (4, [-1 1 1; 2 2 -2; -1 1 -1]),    # N = 32: two copies of N = 16, version 1, side by side
    (5, [-2 2 2; 0 2 -2; -1 1 -1]),    # N = 32: four copies of N = 8, up and to the right
    (6, [-1 1 1; 0 2 -2; -2 2 -2]),    # N = 32: four copies of N = 8, to the front and to the right
    (7, [-2 2 2; 0 1 -1; -2 2 -2]),    # N = 32: four copies of N = 8, up and to the front
    (1, [-2 2 2; 1 1 -1; -2 2 -2]),    # N = 64
    (2, [-2 2 2; 0 2 -2; -2 2 -2]),    # N = 64
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # Heisenberg model on the nearest-neighbor bonds
    H = neighbor_interaction("SdotS", "J", fl; num_distance=1)
    check_couplings(fl, H, Dict("J" => (6, d)))   # bonds per unit cell, length
    file = toml_name(@__DIR__, "hyperhoneycomb", "Heisenberg", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, irreps=false)
    println("wrote ", relpath(file))

    # Kitaev-Heisenberg model: Heisenberg coupling J and Kitaev couplings KX, KY, KZ (X being the
    # "symmetry axis", Y and Z interchangeable)
    HK = neighbor_interaction("SdotS", "J", fl; num_distance=1)
    HK += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0, 0])
    HK += lattice_interaction("SxSx", "KX", fl, 3, 4, [0, 0, 0])
    HK += lattice_interaction("SySy", "KY", fl, 2, 3, [0, 0, 0])
    HK += lattice_interaction("SySy", "KY", fl, 4, 1, [0, 1, 0])   # [0, 1, 0] is Alex's convention, others use [1, 0, 0] here
    HK += lattice_interaction("SzSz", "KZ", fl, 3, 2, [0, 0, 1])
    HK += lattice_interaction("SzSz", "KZ", fl, 4, 1, [1, 0, 0])   # [1, 0, 0] is Alex's convention, others use [0, 1, 0] here
    check_couplings(fl, HK, Dict("J" => (6, d), "KX" => (2, d), "KY" => (2, d), "KZ" => (2, d)))
    file = toml_name(@__DIR__, "hyperhoneycomb", "KitaevHeisenberg", fl, ver)
    write_toml(fl, HK, file; zero_based=true, symmetries=true, irreps=false)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_3d(fl, HK; cpl_dict=Dict("KX" => :blue, "KY" => :red, "KZ" => :green, "J" => :black), annotate_sites=true, annotate_sites_zero_based=true)
        show_figure(f)
    end
end
