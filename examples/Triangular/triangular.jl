using Revise
using Latlib
using GLMakie

function order_a2first_triangular(A::Vector{Float64},B::Vector{Float64})

    BravaisVec = [1 0.5; 0.0 sqrt(3)/2]
    P1 = round.(BravaisVec \ A)
    P2 = round.(BravaisVec \ B)
    return (P1[1] < P2[1]) || ((P1[1] == P2[1]) && (P1[2] < P2[2]))
end

# infinite Bravais lattice for triangular lattice is predefined in lattice/predefined_lattices.jl
infinite_lat = Latlib.triangular

# define finite cluster
L = 12
W = 4
boundary = [L 0; 0 W]
finite_lat = Latlib.FiniteLattice(infinite_lat, boundary, true; order=order_a2first_triangular)

# define nearest neighbor OpSum
H = Latlib.OpSum()
H += Latlib.neighbor_interaction("HB", "J1", finite_lat; num_distance = 1)

# print
GLMakie.activate!(inline=true)
Latlib.plot_opsum(H, finite_lat)

#write_opsum_to_toml!(H,"triangular_J1.toml",index_zero=false)

