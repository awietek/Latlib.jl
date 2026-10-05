# Face-centered cubic lattice: J1-J2 Heisenberg model on periodic clusters
#
# Run from the root of the repository with
#     julia --project=. examples/lattices/fcc/fcc.jl
# It writes a TOML file for every cluster and model below and plots every cluster if show_plot is true.
#
# The space group Fm-3m is symmorphic, so every TOML file contains the symmetry operations and the
# irreducible representations of the cluster (symmetry center: a site; the other inequivalent
# center is the octahedral hole).
# The three-dimensional irreducible representations of the cubic little co-groups are skipped,
# with a warning and a banner in the file.

include(joinpath(@__DIR__, "..", "common.jl"))
show_plot = true   # plot every cluster (close the window to continue)

lattice = fcc
d1, d2 = sqrt(2) / 2, 1.0   # distances of nearest and next-nearest neighbors

# The conventional 2 × 2 × 2 cube (32 sites) is too small for J2: its torus vector (2, 0, 0) is twice
# a J2 bond, so the J2 bonds along it would wrap around the torus.
# periodic clusters: (version, boundary vectors as rows in the lattice basis)
clusters = [
    (1, [3 0 0; 0 3 0; 0 0 3]),   # N = 27: 3 × 3 × 3 primitive cells, point group m-3m
]

for (ver, boundary) in clusters
    fl = FiniteLattice(lattice, boundary, true)

    # J1-J2 Heisenberg model
    H = OpSum()
    H += neighbor_interaction("SdotS", "J1", fl; num_distance=1)
    H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)
    check_couplings(fl, H, Dict("J1" => (6, d1), "J2" => (3, d2)))   # bonds per unit cell, length

    file = toml_name(@__DIR__, "fcc", "J1J2", fl, ver)
    write_toml(fl, H, file; zero_based=true, symmetries=true, origin=LatticeVector(lattice, [0.0, 0.0, 0.0]))
    println("wrote ", relpath(file))

    if show_plot
        f, ax = plot_3d(fl, H; cpl_dict=Dict("J1" => :black, "J2" => :red), annotate_sites=true, annotate_sites_zero_based=true)
        show_figure(f)
    end
end
