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
        periodicity,
        boundary,
        FiniteLatticeVector, 
        to_lattice_basis,
        to_euclidean_basis,
        to_finite_lattice_vector,
        bravais_cells,
        atoms


include("metric.jl")
export EuclideanMetric,
        PeriodicEuclideanMetric,
        distance,
        distance_vector,
        distance_matrix,
        distances,
        neighbors

        
include("opsum.jl")
export Op,
        OpSum,
        isequal,
        isless,
        unique_ops!,
        neighbor_interaction,
        lattice_interaction


include("lattice/predefined_lattices.jl")
export square,
        triangular,
        kagome


include("plots.jl")
export plot,
        plot_opsum



end
