using Revise
using Latlib

# infinite Bravais lattice for triangular lattice is predefined in lattice/predefined_lattices.jl
infinite_lat = triangular

# define finite cluster
L = 12
W = 4
boundary = [L 0; 0 W]
finite_lat = FiniteLattice(infinite_lat, boundary, true)

# define nearest neighbor OpSum
H = OpSum()
H += neighbor_interaction("HB", "J1", finite_lat; num_distance = 1)

# print
GLMakie.activate!()
plot_opsum(H, finite_lat)


