using Revise
using Latlib

# Hyperhoneycomb: 3D lattice with 4 atoms per unit cell
flat_2x2 = FiniteLattice(hyperhoneycomb, [2 0 0; 0 2 0; 0 0 2])
#plot_3d(flat_2x2; show_neighbors=true, show_unit_cell=true)

# N=16 cluster
flat_16_vecs = [
    LatticeVector(hyperhoneycomb, [-1, 1, 1]),  # t1
    LatticeVector(hyperhoneycomb, [1, 1, -1]),  # t2
    LatticeVector(hyperhoneycomb, [-1, 1, -1]), # t3
]

flat_16 = FiniteLattice(flat_16_vecs, true)
#plot_3d(flat_16; show_neighbors=true, show_unit_cell=true)


# N=32 cluster
flat_32_vecs = [
    LatticeVector(hyperhoneycomb, [-1, 1, 1]), # t1
    LatticeVector(hyperhoneycomb, [1, 1, -1]), # t2
    LatticeVector(hyperhoneycomb, [-2, 1, -2]),# t3
]
flat_32 = FiniteLattice(flat_32_vecs, true)
plot_3d(flat_32; show_neighbors=true, show_unit_cell=true)



# print to stdout what would be written to TOML file
opsum = neighbor_interaction("SdotS", "J", fl)
toml_string = write_toml(fl, opsum, "placeholder"; zero_based=true, return_string=true)
#println(toml_string)
