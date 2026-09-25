using Revise
using Latlib
using GLMakie

# infinite Bravais lattice for honeycomb lattice is predefined in lattice/predefined_lattices.jl
infinite_lat = honeycomb

# pick number of sites (atoms) in finite cluster (e.g. 16 or 32)
N = 8

# for most N, there are multiple finite clusters, pick a "version" here starting form 1
ver = 1




fl_vecs = nothing
# -------------------------------------------------
#                   N = 6 cluster               
# -------------------------------------------------
if (N, ver) == (6, 1)
    fl_vecs = [
        LatticeVector(honeycomb, [2, -1]),  # t1
        LatticeVector(honeycomb, [1, 1]),  # t2
    ]
end

# -------------------------------------------------
#                   N = 8 cluster               
# -------------------------------------------------
if (N, ver) == (8, 1)
    fl_vecs = [
        LatticeVector(honeycomb, [2, -2]),  # t1
        LatticeVector(honeycomb, [1, 1]),  # t2
    ]
end





fl = FiniteLattice(fl_vecs, true)

# define nearest neighbor OpSum of Heisenberg-Kitaev model
H = OpSum()
H += neighbor_interaction("SdotS", "J", fl; num_distance = 1)
H += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0]) # connect 1st and 2nd site in the same unit cell
H += lattice_interaction("SySy", "KY", fl, 1, 2, [0, -1]) # connect 1st site in [0, 0] cell to 2nd site in [0, -1] cell
H += lattice_interaction("SzSz", "KZ", fl, 2, 1, [1, 0]) # connect 2nd site in [0, 0] cell to 1st site in [1, 0] cell

# write to TOML
write_toml(fl, H, (@__DIR__) * "/honeycomb-N-$N-ver-$ver.toml"; zero_based=true)

# print
GLMakie.activate!()
plot_opsum(H, fl)

