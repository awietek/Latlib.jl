module Latlib
using LinearAlgebra

include("lattice/lattice.jl")
export Lattice, dimension, natoms

include("lattice/predefined_lattices.jl")
export square, kagome

include("lattice/finite_lattice.jl")
export FiniteLattice, bravais_coordinates, coordinates, boundary_vectors, periodicity_vectors

include("metric.jl")
export distance, distance_matrix, distances, neighbors, periodicity

include("bonds.jl")
export Bond, nearest_neighbor_bonds, lattice_bonds

include("plots.jl")
export plot

end
