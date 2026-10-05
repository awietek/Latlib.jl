# Shastry-Sutherland lattice with the symmetry p4gm of the model: Shastry-Sutherland model
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/shastry_sutherland_non_symmorphic/shastry_sutherland_non_symmorphic.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# In this geometry, the bonds of the squares and the dimers all have length 1, and the lattice has
# the symmetry p4gm of the Shastry-Sutherland model. Since p4gm is non-symmorphic, the TOML files
# contain the symmetry operations of the cluster but no irreducible representations (`irreps=false`,
# as they are not implemented for non-symmorphic space groups yet).

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = shastry_sutherland_non_symmorphic
d = 1.0   # length of all bonds

# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [2 0; 0 2]),    # N = 16: 2 × 2 unit cells, point group 4mm
    (1, [1 2; -2 1]),   # N = 20: tilted, 4
    (1, [2 2; -2 2]),   # N = 32: tilted, 4mm
    (1, [3 0; 0 3]),    # N = 36: 3 × 3 unit cells, 4mm
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # Shastry-Sutherland model: J on the bonds of the squares, Jd on the dimers (atoms 1-2 and 3-4)
    H = OpSum()
    for (atom1, atom2, cell) in [(1, 3, [-1, 0]), (1, 3, [0, 0]), (1, 4, [0, -1]), (1, 4, [0, 0]),
                                 (2, 3, [-1, -1]), (2, 3, [-1, 0]), (2, 4, [-1, -1]), (2, 4, [0, -1])]
        H += lattice_interaction("SdotS", "J", fl, atom1, atom2, cell)
    end
    H += lattice_interaction("SdotS", "Jd", fl, 1, 2, [0, 0])
    H += lattice_interaction("SdotS", "Jd", fl, 3, 4, [0, 0])
    check_couplings(fl, H, Dict("J" => (8, d), "Jd" => (2, d)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "shastry_sutherland_non_symmorphic", "ShastrySutherland", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, irreps=false)
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(H, fl)
        show_figure(f)
    end
end
