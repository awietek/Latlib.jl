using Printf
using IterTools

@doc raw"""
    FiniteLattice

A lattice with a finite number of sites confined by a boundary box.
The boundary box vectors``\mathcal{B}_i`` are specified as integer
multiples of the lattice vectors. The cartesian boundary box vectors
``\mathbf{b}_i`` are then given by,

``\mathbf{b}_i = \mathbf{A}\mathcal{B}_i``

# Arguments
- `lattice::Lattice`: instance of the underlying lattice geometry.
- `boundary::AbstractMatrix`: ``D \times D`` integer matrix whose columns define the boundary box
- `periodicity::Vector{Bool}`: vector defining which boundary direction is periodic
- `order`: functions setting 
"""
struct FiniteLattice
    lattice::Lattice
    boundary::Matrix{Int64}
    periodicity::Vector{Bool}
    order

    function FiniteLattice(lattice::Lattice, boundary::Matrix{Int64}, periodicity::Vector{Bool}, order=nothing)
        dim_vectors = size(lattice.vectors)
        dim_boundary = size(boundary)
        if dim_vectors != dim_boundary
            error(@sprintf "Incompatible dimension lattice vectors %s and boundary matrix %s" dim_vectors dim_boundary)
        end

        dim = dimension(lattice)
        
        if length(periodicity) != dim
            error("Periodicity vector is not of same length as dimension of lattice")
        end
        if order === nothing
            new(lattice, boundary, periodicity, isless)
        else
            new(lattice, boundary, periodicity, order)
        end
    end
end

@doc raw"""
    FiniteLattice(lattice::Lattice, boundary::AbstractMatrix; periodic=true)

Create a finite lattice given an underlying lattice, its boundary, and optionally 
whether or not fully periodic boundary conditions need to be applied

# Arguments
- `lattice::Lattice`: instance of the underlying lattice geometry.
- `boundary::AbstractMatrix`: ``D \times D`` integer matrix whose columns define the boundary box

# Keyword arguments
- `periodic::Bool=true`: flag to set if periodic boundary conditions are used
"""
FiniteLattice(lattice::Lattice, boundary::AbstractMatrix; periodic::Bool=true) =    
    FiniteLattice(lattice, boundary, fill(periodic, dimension(lattice)))


"""
    dimension(flattice::FiniteLattice)

Obtain the dimension of the lattice.
"""
dimension(flattice::FiniteLattice) = size(flattice.lattice.vectors)[1]

"""
    natoms(flattice::FiniteLattice)

Obtain the number of atomic positions
"""
natoms(lattice::FiniteLattice) = size(flattice.lattice.positions)[2]


function Base.show(io::IO, flattice::FiniteLattice)
    dim = dimension(flattice)
    println(io, "FiniteLattice")
    println(io, @sprintf "D = %d" dim)
    for i in 1:dim
        println(io, @sprintf "a%d = %s" i flattice.lattice.vectors[:,i])
    end
    for i in 1:natoms(flattice.lattice)
        println(io, @sprintf "x%d = %s (type: %d)" i flattice.lattice.positions[:,i] flattice.lattice.types[i])
    end
    for i in 1:dim
        println(io, @sprintf "b%d = %s" i flattice.boundary[:,i])
    end
    println(io, @sprintf "periodicity = %s" flattice.periodicity)
end


"""
    bravais_coordinates(flattice::FiniteLattice)

Computes the Bravais coordinates of the finite lattice within
the boundary box, i.e. only the coordinates of the Bravais lattice
but not the coordinates of the atomic positions
"""
function bravais_coordinates(flattice::FiniteLattice)
    lattice = flattice.lattice
    dim = dimension(lattice)

    boundary_cartesian = lattice.vectors * flattice.boundary

    # compute Bravais coordinates within the boundary
    bravais_coords = Vector{Float64}[]
    max_mult = sum(abs.(flattice.boundary), dims=1)
    ranges = [-mm:mm for mm in max_mult]
    for c in product(ranges...)
        bravais_coord = lattice.vectors * convert.(Float64, collect(c))

        # Check if proposed coordinate is in the bounding box
        fractional_coord = boundary_cartesian \ bravais_coord
        incell = all(fractional_coord .> -1e-6) && all(fractional_coord .< 1 - 1e-6)

        if incell
            push!(bravais_coords, bravais_coord)
        end
    end
    sort!(bravais_coords, lt=flattice.order)
    return hcat(bravais_coords...)
end

"""
    coordinates(flattice::FiniteLattice)

Computes all coordinates of the finite lattice within
the boundary box.
"""
function coordinates(flattice::FiniteLattice)
    lattice = flattice.lattice
    bravais_coords = bravais_coordinates(flattice)

    # generate (unsorted) list of all coordinates
    all_coords = Vector{Float64}[]
    for bravais_coord in eachcol(bravais_coords)
        for fractional_basis_coord in eachcol(lattice.positions)
            basis_coord = lattice.vectors * fractional_basis_coord
            push!(all_coords, bravais_coord + basis_coord)
        end
    end

    # sort the coordinates w.r.t. given ordering
    sort!(all_coords, lt=flattice.order)
    return hcat(all_coords...)
end

"""
    boundary_vectors(flattice::FiniteLattice)

Computes the vectors of the boundary box in Cartesian coordinates
"""
function boundary_vectors(flattice::FiniteLattice)
    return flattice.lattice.vectors * flattice.boundary
end

"""
    periodicity_vectors(flattice::FiniteLattice)

Computes the vectors of the periodic directions in Cartesian coordinates
"""
function periodicity_vectors(flattice::FiniteLattice)
    return boundary_vectors(flattice)[:, flattice.periodicity]
end


"""
    neighbors(flattice::FiniteLattice; num_distance::Integer=1)

Computes which pairs of coordinates are neighbors

# Arguments
- `flattice::FiniteLattice`: finite lattice 

# Keyword arguments
- `num_distance::Integer=1`: at which distance neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.
"""
function neighbors(flattice::FiniteLattice; num_distance::Integer=1)
    return neighbors(coordinates(flattice);
                     num_distance=num_distance,
                     periodicity_vectors=periodicity_vectors(flattice))
end


"""
    distance(x1::AbstractVector, x2::AbstractVector, flattice::FiniteLattice)

Computes the distance between two points given the periodicity of the
finite lattice.
"""
function distance(x1::AbstractVector, x2::AbstractVector, flattice::FiniteLattice)
    return distance(x1, x2; periodicity_vectors=periodicity_vectors(flattice))
end

function distance_vector(x1::AbstractVector, x2::AbstractVector, flattice::FiniteLattice)
    return distance_vector(x1, x2; periodicity_vectors=periodicity_vectors(flattice))
end

"""
    distances(flattice::FiniteLattice)

Computes which unique values of distances are present on the lattice
"""
function distances(flattice::FiniteLattice)
    return distances(coordinates(flattice);
                     periodicity_vectors=periodicity_vectors(flattice))
end

