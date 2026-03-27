using Revise
using Latlib
using HDF5 # for retrieving spin data from h5 file (not required in general)


# this example demonstrates how to plot 3D spins into a 3D lattice structure


function main()

    # ----- define finite lattice
    N = 16 
    fl_vecs = [
        LatticeVector(hyperhoneycomb, [-1, 1, 1]),  # t1
        LatticeVector(hyperhoneycomb, [1, 1, -1]),  # t2
        LatticeVector(hyperhoneycomb, [-1, 1, -1]), # t3
    ]
    fl = FiniteLattice(fl_vecs, true)

    # ---- retrieve spin data (must give Vector{N, Vector{3, Float64}} for 3D plotting)
    spin_data = get_spin_data(N)

    # ---- define relevant Hamiltonian (this makes the final plot more accessible)
    opsum = neighbor_interaction("SdotS", "J", fl)
    # define Kitaev X, Y, Z interactions (X being the "symmetry axis", and Y-Z being interachangable)
    opsum += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0, 0])
    opsum += lattice_interaction("SxSx", "KX", fl, 3, 4, [0, 0, 0])
    opsum += lattice_interaction("SySy", "KY", fl, 2, 3, [0, 0, 0])
    opsum += lattice_interaction("SySy", "KY", fl, 4, 1, [0, 1, 0]) # [0, 1, 0] is Alex's convention, others use [1, 0, 0] here
    opsum += lattice_interaction("SzSz", "KZ", fl, 3, 2, [0, 0, 1])
    opsum += lattice_interaction("SzSz", "KZ", fl, 4, 1, [1, 0, 0]) # [1, 0, 0] is Alex's convention, others use [0, 1, 0] here

    # ---- plot spins in 3D lattice with Hamiltonian interactions
    plot_3d(fl, opsum, spin_data; 
        # ----- keywords for plot_3d(fl, opsum)
        cpl_dict = Dict("KX" => :blue, "KY" => :red, "KZ" => :green),
        # ----- keywords for plot_3d(fl)
        #show_unit_cell=true,
        annotate_sites=true,
        annotate_sites_zero_based=true,
        spin_length_multiplier=1.5, # control size of spin arrows to match with lattice
        spin_color=:gray,
        site_marksize = 0.0, # hide lattice sites for better visibility of spins
        #draw_periodic_flattice=true,
        #draw_periodic_flattice_shifts=[(1,0,0),(0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1, 1, 1)],
        scale_factor=2.0,
        )

end

# this is an example function to retrieve spin data and needs to be adjusted by the user!
function get_spin_data(N::Int64) :: Vector{Vector{Float64}}
    h5_path = "local_spins_N-16-n-144-D-6.h5" # this file is not included in the repo!
    i = 125# pick a point in phaes diagram of hyperhoneycomb HB-Kitaev model
    f = h5open(h5_path, "r")
    spin_data = Vector{Vector{Float64}}()
    for j in 0:(N-1)
        push!(spin_data, real.(read(f["$i/$j"])))
    end
    return spin_data
end







main()