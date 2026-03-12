using Revise
using Latlib

# Hyperhoneycomb: 3D lattice with 4 atoms per unit cell
fl = FiniteLattice(hyperhoneycomb, [2 0 0; 0 2 0; 0 0 2])

#plot_3d(fl; show_neighbors=true, show_unit_cell=true)

# define Heisenberg interactions between nearest neighbors
opsum = neighbor_interaction("SdotS", "J", fl)

toml_string = write_toml(fl, opsum, "placeholder"; index_zero=true, return_string=true)
println(toml_string)
