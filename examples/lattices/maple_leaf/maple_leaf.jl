# Maple-leaf lattice: Heisenberg model with hexagon, triangle and dimer couplings
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/maple_leaf/maple_leaf.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# Every site has five nearest neighbors at distance 1, on bonds of three kinds: bonds of the
# hexagons (Jhexagon), of the triangles (Jtriangle) and dimers (Jdimer), as in the TOML files of
# test/data/legacy. Each kind is invariant under the space group p6, which is symmorphic, so every
# TOML file contains the symmetry operations and the irreducible representations of the cluster
# (symmetry center: the center of a hexagon).

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = maple_leaf
d = 1.0   # length of all nearest-neighbor bonds

# periodic clusters: (version, boundary vectors as rows in the lattice basis); the versions of the
# clusters with 12 to 54 sites are those of the files in test/data/legacy
clusters = [
    (1, [1 0; 0 2]),    # N = 12: point group 2
    (1, [1 1; 1 -2]),   # N = 18: 6
    (2, [1 1; 2 -2]),   # N = 24: 2
    (3, [1 1; 3 -3]),   # N = 36: 2
    (1, [1 2; -3 1]),   # N = 42: 6
    (1, [3 0; 0 3]),    # N = 54: 6
]

hexagon = [(1, 4, [0, -1]), (1, 6, [0, 0]), (2, 3, [0, 0]), (2, 5, [-1, 0]), (3, 6, [-1, 1]), (4, 5, [0, 0])]
triangle = [(1, 3, [0, -1]), (1, 5, [-1, 0]), (2, 4, [0, 0]), (2, 6, [0, 0]), (3, 5, [-1, 1]), (4, 6, [0, 0])]
dimer = [(1, 2, [0, 0]), (3, 4, [0, 0]), (5, 6, [0, 0])]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # Heisenberg model with couplings Jhexagon, Jtriangle and Jdimer: (atom 1, atom 2, cell of atom 2)
    H = OpSum()
    for (coupling, bonds) in [("Jhexagon", hexagon), ("Jtriangle", triangle), ("Jdimer", dimer)]
        for (atom1, atom2, cell) in bonds
            H += lattice_interaction("SdotS", coupling, fl, atom1, atom2, cell)
        end
    end
    check_couplings(fl, H, Dict("Jhexagon" => (6, d), "Jtriangle" => (6, d), "Jdimer" => (3, d)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "maple_leaf", "JhexagonJtriangleJdimer", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(H, fl)
        show_figure(f)
    end
end
