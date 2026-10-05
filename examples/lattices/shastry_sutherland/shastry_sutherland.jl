# Shastry-Sutherland lattice: Shastry-Sutherland model on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/shastry_sutherland/shastry_sutherland.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The sites of `shastry_sutherland` form a square lattice with a unit cell of four sites; the
# dimers of the Shastry-Sutherland model exist only in the couplings. The TOML files contain the
# symmetry operations and irreducible representations of the sites (space group p4mm, with the
# translations restricted to the lattice vectors of the four-site unit cell, and a site as symmetry
# center). Only part of these operations are symmetries of the Shastry-Sutherland model, whose space
# group is p4gm: the fourfold rotations about a site, for example, are not. For a geometry with the
# symmetry of the model, see ../shastry_sutherland_non_symmorphic.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = shastry_sutherland
d, dd = 0.5, sqrt(2) / 2   # lengths of the square bonds and of the dimer bonds

# periodic clusters: (version, boundary vectors as rows in the lattice basis, in units of the
# unit cell of four sites)
clusters = [
    (1, [2 0; 0 2]),    # N = 16: 2 × 2 unit cells
    (1, [1 2; -2 1]),   # N = 20: tilted
    (1, [2 2; -2 2]),   # N = 32: tilted
    (1, [3 0; 0 3]),    # N = 36: 3 × 3 unit cells
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # Shastry-Sutherland model: J on the bonds of the square lattice, Jd on the dimers
    H = OpSum()
    H += neighbor_interaction("SdotS", "J", fl; num_distance=1)
    H += lattice_interaction("SdotS", "Jd", fl, 1, 4, [0, 0])    # dimer of atoms 1 and 4 in the same cell
    H += lattice_interaction("SdotS", "Jd", fl, 3, 2, [1, -1])   # dimer of atom 3 and atom 2 in cell [1, -1]
    check_couplings(fl, H, Dict("J" => (8, d), "Jd" => (2, dd)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "shastry_sutherland", "ShastrySutherland", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, origin=LatticeVector(lattice, [0.0, 0.0]))
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_opsum(H, fl)
        show_figure(f)
    end
end
