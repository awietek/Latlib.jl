# Honeycomb lattice: J1-J2 Heisenberg and Kitaev-Heisenberg models on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/honeycomb/honeycomb.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The space group p6mm is symmorphic, so every TOML file contains the symmetry operations and the
# irreducible representations of the cluster (symmetry center: the center of a hexagon). For the
# Kitaev-Heisenberg model, these are the symmetries of the lattice: the translations are symmetries
# of the Hamiltonian, but the point-group operations exchange the bonds KX, KY and KZ, so they are
# symmetries only together with a corresponding transformation of the spin components.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = honeycomb
d1, d2 = 1 / sqrt(3), 1.0   # distances of nearest and next-nearest neighbors

# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [3 0; 0 3]),    # N = 18: 3 × 3, point group 6mm
    (1, [2 2; -2 4]),   # N = 24: tilted, 6mm
    (1, [4 0; 0 4]),    # N = 32: 4 × 4, 6mm
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # J1-J2 Heisenberg model
    H = OpSum()
    H += neighbor_interaction("SdotS", "J1", fl; num_distance=1)
    H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)
    check_couplings(fl, H, Dict("J1" => (3, d1), "J2" => (6, d2)))   # bonds per unit cell, length
    file = toml_name(@__DIR__, "honeycomb", "J1J2", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true)
    println("wrote ", relpath(file))

    # Kitaev-Heisenberg model: Heisenberg coupling J and Kitaev couplings KX, KY, KZ on the three
    # nearest-neighbor bonds of every site
    HK = OpSum()
    HK += neighbor_interaction("SdotS", "J", fl; num_distance=1)
    HK += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0])    # atom 1 to atom 2 in the same cell
    HK += lattice_interaction("SySy", "KY", fl, 1, 2, [0, -1])   # atom 1 to atom 2 in cell [0, -1]
    HK += lattice_interaction("SzSz", "KZ", fl, 2, 1, [1, 0])    # atom 2 to atom 1 in cell [1, 0]
    check_couplings(fl, HK, Dict("J" => (3, d1), "KX" => (1, d1), "KY" => (1, d1), "KZ" => (1, d1)))
    file = toml_name(@__DIR__, "honeycomb", "KitaevHeisenberg", fl, ver)
    write_toml(fl, HK, file; zero_based=true, symmetries=true)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(HK, fl)
        show_figure(f)
    end
end
