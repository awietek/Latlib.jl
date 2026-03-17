using Revise
using Latlib
using GLMakie

# Shastry-Sutherland lattice is among predefined constants
ss_lattice = shastry_sutherland


# Finite cluster with open boundaries (periodicity false) on one side
L = 6
W = 4
boundary = [L÷2 0; 0 W÷2]
fl = FiniteLattice(ss_lattice, boundary, [false, true], order_xfirst)


# Define Interactions
H = OpSum()
H += neighbor_interaction("SdotS", "J", fl; num_distance=1)
# Jd dimer bond 1: atom 1 and atom 4 within the same unit cell
H += lattice_interaction("SdotS", "Jd", fl, 1, 4, [0, 0])
# Jd dimer bond 2: atom 3 (origin cell) to atom 2 (cell offset [1,-1])
H += lattice_interaction("SdotS", "Jd", fl, 3, 2, [1, -1])


# write cylinder to TOML file
toml_path = joinpath(@__DIR__, "shastry_sutherland-cylinder-L-$L-W-$W.toml")
write_toml(fl, H, toml_path; zero_based=false)


# Plot
GLMakie.activate!()
plot_opsum(H, fl)