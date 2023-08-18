module Latlib

using LinearAlgebra

include("lattice.jl")
export Lattice, dimension, coordinates

include("metric.jl")
export distance, neighbors, periodicity

include("bonds.jl")
export Bond, nearest_neighbor_bonds, lattice_bonds

include("plots.jl")
export plot

end
