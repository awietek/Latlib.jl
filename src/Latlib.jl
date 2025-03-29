module Latlib
using LinearAlgebra

include("lattice/lattice.jl")
export Lattice, dimension, natoms, inlattice, vector_position

include("lattice/predefined_lattices.jl")
export square, kagome

include("lattice/finite_lattice.jl")
export FiniteLattice, bravais_coordinates, coordinates, boundary_vectors, periodicity_vectors

include("metric.jl")
export distance, distance_matrix, distances, neighbors, periodicity

include("opsum.jl")
export Op, OpSum, unique_ops!, neighbor_bonds, lattice_bonds, write_OpSum_to_toml!

include("plots.jl")
export plot

end
