using Revise
using Latlib

# pick number of sites (atoms) in finite cluster (e.g. 16 or 32)
N = 16

# for most N, there are multiple finite clusters, pick a "version" here starting form 1
ver = 1

fl_vecs = nothing

# -------------------------------------------------
#                   N = 16 cluster                
# -------------------------------------------------
if (N, ver) == (16, 1)
    fl_vecs = [
        LatticeVector(hyperhoneycomb, [-1, 1, 1]),  # t1
        LatticeVector(hyperhoneycomb, [1, 1, -1]),  # t2
        LatticeVector(hyperhoneycomb, [-1, 1, -1]), # t3
    ]
end


# -------------------------------------------------
#                   N = 32 clusters               
# -------------------------------------------------

if (N, ver) == (32, 2) # starts with 2 for "historic" reasons
    # two copies of 16v1 stacked in front of each other
    fl_vecs = [
        LatticeVector(hyperhoneycomb, [-1, 1, 1]), # t1
        LatticeVector(hyperhoneycomb, [1, 1, -1]), # t2
        LatticeVector(hyperhoneycomb, [-2, 2, -2]),# t3
    ]
elseif (N, ver) == (32, 3)
    # two copies of 16v1 stacked on top of each other
    fl_vecs = [
        LatticeVector(hyperhoneycomb, [-2, 2, 2]), # t1
        LatticeVector(hyperhoneycomb, [1, 1, -1]), # t2
        LatticeVector(hyperhoneycomb, [-1, 1, -1]),# t3
    ]
elseif (N, ver) == (32, 4)
    # two copies of 16v1 stacked side by side
    fl_vecs = [
        LatticeVector(hyperhoneycomb, [-1, 1, 1]), # t1
        LatticeVector(hyperhoneycomb, [2, 2, -2]), # t2
        LatticeVector(hyperhoneycomb, [-1, 1, -1]),# t3
    ]
end

# -------------------------------------------------
#                   N = 64 clusters               
# -------------------------------------------------


if (N, ver) == (64, 1)
    fl_vecs = [
        LatticeVector(hyperhoneycomb, [-2, 2, 2]), # t1
        LatticeVector(hyperhoneycomb, [1, 1, -1]), # t2
        LatticeVector(hyperhoneycomb, [-2, 2, -2]),# t3
    ]
end



for (i, v) in enumerate(fl_vecs)
    println("t$i = ", to_euclidean_basis(v))
end



# throw error message for undefined N and version
if isnothing(fl_vecs)
    error("No finite lattice defined for N = $N and version = $ver. Try e.g. N=16, or N=32 and ver=1.")
end

# -------------------------------------------------
#              FINITE LATTICE & PLOT                
# -------------------------------------------------
fl = FiniteLattice(fl_vecs, true)

#=
plot_3d(fl; 
    show_neighbors=true,
    #show_unit_cell=true,
    annotate_sites=true,
    annotate_sites_zero_based=true,
    #draw_periodic_flattice=true,
    #draw_periodic_flattice_shifts=[(1,0,0),(0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1, 1, 1)],
    )
=#


# -------------------------------------------------
#                   INTERACTIONS                
# -------------------------------------------------

# define Heisenberg interaction between all nearest neighbors
opsum = neighbor_interaction("SdotS", "J", fl)
# define Kitaev X, Y, Z interactions (X being the "symmetry axis", and Y-Z being interachangable)
opsum += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0, 0])
opsum += lattice_interaction("SxSx", "KX", fl, 3, 4, [0, 0, 0])
opsum += lattice_interaction("SySy", "KY", fl, 2, 3, [0, 0, 0])
opsum += lattice_interaction("SySy", "KY", fl, 4, 1, [0, 1, 0]) # [0, 1, 0] is Alex's convention, others use [1, 0, 0] here
opsum += lattice_interaction("SzSz", "KZ", fl, 3, 2, [0, 0, 1])
opsum += lattice_interaction("SzSz", "KZ", fl, 4, 1, [1, 0, 0]) # [1, 0, 0] is Alex's convention, others use [0, 1, 0] here

# write to TOML
write_toml(fl, opsum, (@__DIR__) * "/hyperhoneycomb-N-$N-ver-$ver.toml"; zero_based=true)

# draw opsum into lattice


plot_3d(fl, opsum; 
    # ----- keywords for plot_3d(fl, opsum)
    cpl_dict = Dict("KX" => :blue, "KY" => :red, "KZ" => :green),
    # ----- keywords for plot_3d(fl)
    #show_unit_cell=true,
    annotate_sites=true,
    annotate_sites_zero_based=true,
    draw_periodic_flattice=true,
    draw_periodic_flattice_shifts=[(1,0,0),(0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1, 1, 1)],
    scale_factor=2.0,
    )
