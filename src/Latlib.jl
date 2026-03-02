module Latlib
using LinearAlgebra

include("vectors.jl")
export EuclideanVector

include("lattice/lattice.jl")
export Lattice,
        LatticeVector,
        dim,
        natoms,
        to_lattice_basis,
        to_euclidean_basis,
        in_lattice


#=
include("lattice/orders.jl")
export order_xfirst

include("lattice/predefined_lattices.jl")
export square, kagome

include("lattice/finite_lattice.jl")
export FiniteLattice, bravais_coordinates, coordinates, boundary_vectors, periodicity_vectors

include("metric.jl")
export distance, distance_matrix, distances, neighbors, periodicity

include("opsum.jl")
export Op, OpSum, unique_ops!, neighbor_bonds, lattice_bonds, write_opsum_to_toml!

include("plots.jl")
export plot, plot_opsum
=#

end
