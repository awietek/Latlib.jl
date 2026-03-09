module Latlib
using LinearAlgebra

include("utils.jl")


include("vectors.jl")
export EuclideanVector,
        equal,
        plus,
        minus,
        mul,
        norm


include("lattice/lattice.jl")
export Lattice,
        LatticeVector,
        dim,
        natoms,
        to_lattice_basis,
        to_euclidean_basis,
        in_lattice



include("lattice/orders.jl")
export order_xfirst


include("lattice/finite_lattice.jl")
export FiniteLattice,
        dim,
        natoms,
        FiniteLatticeVector, 
        to_lattice_basis,
        to_euclidean_basis,
        to_finite_lattice_vector,
        bravais_cells,
        atom_coords,
        boundary_vecs


include("metric.jl")
export EuclideanMetric,
        PeriodicEuclideanMetric,
        distance,
        distance_vector,
        distance_matrix,
        distances,
        neighbors


#=
include("lattice/predefined_lattices.jl")
export square, kagome




include("opsum.jl")
export Op, OpSum, unique_ops!, neighbor_bonds, lattice_bonds, write_opsum_to_toml!

include("plots.jl")
export plot, plot_opsum
=#

end
