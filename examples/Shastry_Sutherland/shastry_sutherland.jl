using Revise
using Latlib
using GLMakie

### Lattice Dimensions
L = 6
W = 4

# Shastry-Sutherland: square Bravais lattice with 4 atoms per unit cell
A_ss = [1.0 0.0; 0.0 1.0]
positions_ss = [0.0 0.0;
                0.0 0.5;
                0.5 0.0;
                0.5 0.5]
ss_lattice = Lattice(A_ss, positions_ss)

# Finite cluster with periodic boundaries
boundary = [L÷2 0; 0 W÷2]
finite_lat = FiniteLattice(ss_lattice, boundary, true)

# Define Interactions
H = OpSum()
H += neighbor_interaction("HB", "J", finite_lat; num_distance=1)
# Jd dimer bond 1: atom 1 and atom 4 within the same unit cell
H += lattice_interaction("HB", "Jd", finite_lat, 1, 4, [0, 0])
# Jd dimer bond 2: atom 3 (origin cell) to atom 2 (cell offset [1,-1])
H += lattice_interaction("HB", "Jd", finite_lat, 3, 2, [1, -1])

# Plot
GLMakie.activate!(inline=true)
plot_opsum(H, finite_lat)

# write_opsum_to_toml!(H, "shastry_sutherland.toml", index_zero=false)
