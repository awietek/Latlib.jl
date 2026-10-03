using Revise
using Latlib
using GLMakie

# infinite Bravais lattice for square lattice is predefined in lattice/predefined_lattices.jl
infinite_lat = square

# pick number of sites (atoms) in finite cluster (e.g. 16 or 32)
N = 16

# for most N, there are multiple finite clusters, pick a "version" here starting form 1
ver = 1




fl_vecs = nothing
# -------------------------------------------------
#                   N = 4 cluster               
# -------------------------------------------------
if (N, ver) == (4, 1)
    fl_vecs = [
        LatticeVector(square, [2, 0]),  # t1
        LatticeVector(square, [0, 2]),  # t2
    ]
end

# -------------------------------------------------
#                   N = 8 cluster               
# -------------------------------------------------
if (N, ver) == (8, 1)
    fl_vecs = [
        LatticeVector(square, [2, 0]),  # t1
        LatticeVector(square, [0, 4]),  # t2
    ]
end

# -------------------------------------------------
#                   N = 16 cluster               
# -------------------------------------------------
if (N, ver) == (16, 1)
    fl_vecs = [
        LatticeVector(square, [4, 0]),  # t1
        LatticeVector(square, [0, 4]),  # t2
    ]
end



fl = FiniteLattice(fl_vecs, true)

# define nearest neighbor OpSum of J_vert, J_horiz interactions
H = OpSum()
# H += neighbor_interaction("SdotS", "J", fl; num_distance = 1)
H += lattice_interaction("SdotS", "Jv", fl, 1, 1, [0, 1]) 
H += lattice_interaction("SdotS", "Jh", fl, 1, 1, [1, 0])

# write to TOML
write_toml(fl, H, (@__DIR__) * "/square-N-$N-ver-$ver.toml"; zero_based=true)

# print
GLMakie.activate!()
plot_opsum(H, fl)

