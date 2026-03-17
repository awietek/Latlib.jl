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
        positions,
        get_position,
        lattice_vecs,
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
        periodic_boundary,
        lattice_vecs,
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
        kagome,
        shastry_sutherland,
        hyperhoneycomb


include("plots.jl")
export plot,
        plot_3d,
        plot_opsum


include("write.jl")
export write_toml,
        toml_coordinates,
        toml_interactions,
        toml_lattice,
        toml_metadata

include("read.jl")
export read_toml_interaction


end
