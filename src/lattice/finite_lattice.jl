import Base: round
using Printf
using LinearAlgebra
#using IterTools

@doc raw"""
    FiniteLattice

A lattice with a finite number of sites confined by a boundary box.
The boundary box vectors``\mathcal{B}_i`` are specified as integer
multiples of the lattice vectors. The cartesian boundary box vectors
``\mathbf{b}_i`` are then given by,

``\mathbf{b}_i = \mathbf{A}\mathcal{B}_i``

# Arguments
- `lattice::Lattice`: The underlying (infinite) lattice geometry.
- `boundary::Matrix{Int64}`: ``D \times D`` integer matrix whose rows define the boundary box vectors in terms of the basis of the lattice.
- `periodicity::Vector{Bool}`: Vector defining which boundary direction is periodic.
- `order`: Functions determining how sites are ordered. 
"""
struct FiniteLattice
    lattice::Lattice
    boundary::Matrix{Int64}
    periodicity::Vector{Bool}
    tol::Float64
    order

    function FiniteLattice(lattice::Lattice, boundary::Matrix{Int64}, periodicity::Vector{Bool}, order=nothing)
        lat_dim = lattice.dim

        # check if boundary is D x D matrix
        if size(boundary) != (lat_dim, lat_dim)
            error(@sprintf "Boundary matrix for `FiniteLattice` instances must be square with dimension matching the lattice (%d). Got boundary of size %s" lat_dim size(boundary))
        end
        
        # check if periodicity is length D vector
        if length(periodicity) != lat_dim
            error(@sprintf "Length of periodicity vector (%d) for `FiniteLattice` instances must match the lattice dimension (%d)" length(periodicity) lat_dim)
        end

        # determine tolerance based on tolarance of lattice and size of largest boundary vector
        max_boundary_vec_length = maximum([LinearAlgebra.norm(boundary[i, :]) for i in 1:lat_dim])
        tol = lattice.tol * max_boundary_vec_length

        # check if order is defined, else default to "isless"
        if isnothing(order)
            return new(lattice, boundary, periodicity, tol, isless)
        else
            return new(lattice, boundary, periodicity, tol, order)
        end
    end

    # alternative constructor to allow setting all boundary conditions at once (true by default)
    function FiniteLattice(lattice::Lattice, boundary::Matrix{Int64}, periodic::Bool=true; order=nothing)
        return FiniteLattice(lattice, boundary, fill(periodic, lattice.dim), order)
    end

    # alternative constructor for use of LatticeVector types for boundary
    function FiniteLattice(boundary::Vector{LatticeVector}, periodicity::Vector{Bool}, order=nothing)
        lattice = boundary[1].lattice
        lat_dim = lattice.dim
        for v in boundary
            if v.lattice != lattice
                error("All boundary vectors must belong to the same lattice")
            end
            if !in_lattice(v)
                error("One of the `LatticeVector` objects provided as boundary is not contained in the lattice (has no integer components).")
            end
        end
        # after this check dimension of LatticeVectors is automatically D

        # check if lattice dimension matches periodicity
        if length(periodicity) != lat_dim
            error(@sprintf "Length of periodicity vector (%d) for `FiniteLattice` instances must match the lattice dimension (%d)" length(periodicity) lat_dim)
        end

        boundary_mat = convert(Matrix{Int64}, hcat([v.coords for v in boundary]...)')
        return FiniteLattice(lattice, boundary_mat, periodicity, order)
    end

    # alternative constructor with LatticeVectors and periodicity set everywhere
    function FiniteLattice(boundary::Vector{LatticeVector}, periodicity::Bool=true; order=nothing)
        return FiniteLattice(boundary, fill(periodicity, length(boundary)), order)
    end

end

@doc """
    dim(flattice::FiniteLattice)

Obtain the dimension of the lattice.
"""
dim(flattice::FiniteLattice) = flattice.lattice.dim

@doc """
    natoms(flattice::FiniteLattice)

Obtain the number of atomic positions
"""
natoms(flattice::FiniteLattice) = flattice.lattice.natoms

@doc """
    positions(flattice::FiniteLattice)

Obtain Vector{LatticeVector} objects containing positions of all atoms IN THE UNIT CELL of the underlying lattice.
"""
positions(flattice::FiniteLattice) = positions(flattice.lattice)

@doc """
    get_position(flattice::FiniteLattice, atom_index::Int64)

Obtain the position of a specific atom IN THE UNIT CELL of the underlying lattice as a `LatticeVector` object.
"""
get_position(flattice::FiniteLattice, idx::Int64) = get_position(flattice.lattice, idx)

@doc """
    periodicity(flattice::FiniteLattice)

Obtain the periodicity (::Bool) for each boundary direction.
"""
periodicity(flattice::FiniteLattice) = flattice.periodicity

@doc """
    boundary(flattice::FiniteLattice)

Obtain the boundary vectors of the finite lattice as Vector{LatticeVector}.
"""
boundary(flattice::FiniteLattice) = [LatticeVector(flattice.lattice, flattice.boundary[i, :]) for i in 1:dim(flattice)]

# nice printing function
function Base.show(io::IO, flattice::FiniteLattice)
    d = flattice.lattice.dim
    println(io, "FiniteLattice")
    println(io, @sprintf "D = %d" d)
    for i in 1:d
        println(io, @sprintf "a%d = %s" i flattice.lattice.A[i, :])
    end
    for i in 1:natoms(flattice.lattice)
        println(io, @sprintf "x%d = %s (type: %d)" i flattice.lattice.positions[i, :] flattice.lattice.types[i])
    end
    for i in 1:d
        println(io, @sprintf "b%d = %s" i flattice.boundary[i, :])
    end
    println(io, @sprintf "periodicity = %s" flattice.periodicity)
end








@doc raw"""
    FiniteLatticeVector(flattice::FiniteLattice, coords::Vector{Float64})

    A vector given in terms of the boundary (not lattice) vectors of a `FiniteLattice`.

    Inputs:
    - `flattice::FiniteLattice`: the finite lattice to which the vector belongs, defining the boundary vectors in terms of which the coordinates are given.
    - `coords::Vector{Float64}`: the coordinates of the vector in terms of the boundary vectors of the finite lattice.

"""
struct FiniteLatticeVector
    flattice::FiniteLattice
    coords::Vector{Float64} 

    function FiniteLatticeVector(flattice::FiniteLattice, coords::Vector{Float64})
        if length(coords) != dim(flattice)
            error(@sprintf "Length of coordinates vector (%d) must match the dimension of the finite lattice (%d)" length(coords) dim(flattice))
        end
        return new(flattice, coords)
    end

    # alternative constructor to allow for integer coordinates
    function FiniteLatticeVector(flattice::FiniteLattice, coords::Vector{Int64})
        return FiniteLatticeVector(flattice, convert(Vector{Float64}, coords))
    end

    # alternative constructor from LatticeVector
    function FiniteLatticeVector(flattice::FiniteLattice, v::LatticeVector)
        if !(v.lattice == flattice.lattice)
            error("The `LatticeVector` provided to construct a `FiniteLatticeVector` must belong to the same lattice as the `FiniteLattice`.")
        end
        return FiniteLatticeVector(flattice, inv(float.(flattice.boundary)') * v.coords)
    end

    # alternative constructor from LatticeVector with boundary only
    function FiniteLatticeVector(v::LatticeVector, boundary::Matrix{Int64})
        flattice = FiniteLattice(v.lattice, boundary)
        return FiniteLatticeVector(flattice, v)
    end

    # alternative constructor from EuclideanVector
    function FiniteLatticeVector(flattice::FiniteLattice, v::EuclideanVector)
        return FiniteLatticeVector(flattice, to_lattice_basis(flattice.lattice, v))
    end
end


@doc """
    to_lattice_vector(v::FiniteLatticeVector)

    Convert a `FiniteLatticeVector` to a `LatticeVector`.
"""
function to_lattice_basis(v::FiniteLatticeVector) ::LatticeVector
    return LatticeVector(v.flattice.lattice, v.flattice.boundary' * v.coords)
end

@doc """
    to_euclidean_basis(v::FiniteLatticeVector)

    Convert a `FiniteLatticeVector` to a `EuclideanVector`.
"""
function to_euclidean_basis(v::FiniteLatticeVector) ::EuclideanVector
    return to_euclidean_basis(to_lattice_basis(v))
end

@doc """
    to_finite_lattice_basis(v::LatticeVector, flattice::FiniteLattice)

    Convert a `LatticeVector` to a `FiniteLatticeVector`.
"""
function to_finite_lattice_basis(flattice::FiniteLattice, v::LatticeVector) ::FiniteLatticeVector
    return FiniteLatticeVector(flattice, v)
end

@doc """
    round(v::FiniteLatticeVector)

    Returns the nearest `FiniteLatticeVector` with integer coordinates to the input vector `v`,
    i.e., the best approximation of v as as multiple of boundary vectors.
"""
function round(v::FiniteLatticeVector) ::FiniteLatticeVector
    return FiniteLatticeVector(v.flattice, round.(v.coords))
end

@doc """
    on_boundary(v::FiniteLatticeVector)

    Checks if a `FiniteLatticeVector` is a multiple of the boundary vectors, i.e., if it corresponds to a torus vector.
"""
function on_boundary(v::FiniteLatticeVector) :: Bool
    return all(is_whole.(v.coords; atol=v.flattice.tol))
end

@doc """
    in_lattice(v::FiniteLatticeVector)

    Checks if a `FiniteLatticeVector` is inside the underlying lattice.
"""
function in_lattice(v::FiniteLatticeVector) :: Bool
    return in_lattice(to_lattice_basis(v))
end


# printing function for `FiniteLatticeVector`
function Base.show(io::IO, v::FiniteLatticeVector)
    print(io, "FiniteLatticeVector(", v.coords, ")")
end

















@doc """
    bravais_cells(flattice::FiniteLattice)

    Computes the Bravais coordinates (integer multiples of the lattice vectors)
    inside the boundary box (not the real space coordinates of the atoms).
"""
function bravais_cells(flattice::FiniteLattice) :: Vector{LatticeVector}
    lattice = flattice.lattice
    ndim = dim(lattice)

    unit_cell_area = 1.0 # work in lattice units here
    flattice_area = abs(det(flattice.boundary))
    num_cells = flattice_area / unit_cell_area
    if !is_whole(num_cells; atol=lattice.tol)
        error(@sprintf "The area of the finite lattice must be an integer multiple of the unit cell area! Finite lattice area = %f, unit cell area = %f, ratio = %f." flattice_area unit_cell_area num_cells)
    end
    num_cells = round(Int, num_cells) # this is how many cells we should find inside the boundaries

    # get rectangular bounding box of finite lattice (all 2^D corners of the parallelepiped)
    flattice_corners = LatticeVector[]
    for bits in Iterators.product(ntuple(_ -> 0:1, ndim)...)
        corner_coords = sum(bits[i] * flattice.boundary[i, :] for i in 1:ndim)
        push!(flattice_corners, LatticeVector(lattice, float.(corner_coords)))
    end
    dim_min = Int64[]
    dim_max = Int64[]
    for i in 1:ndim
        push!(dim_min, round(Int64, floor(minimum([v.coords[i] for v in flattice_corners]))))
        push!(dim_max, round(Int64, ceil(maximum([v.coords[i] for v in flattice_corners]))))
    end
    ranges = [dim_min[i]:dim_max[i] for i in 1:ndim]

    # generate trial points inside the bounding box
    bravais_coords = Vector{LatticeVector}()
    for brav_tuple in Iterators.product(ranges...)
        trial_point = LatticeVector(lattice, collect(brav_tuple))
        if !in_lattice(trial_point) # consistency check, should always be true for integer coordinates
            throw("bravais_cells: generated trial point that is not in the lattice, this is a bug, please report!")
        end
        # check if trial point is inside the finite lattice
        trial_point_fl = to_finite_lattice_basis(flattice, trial_point)
        tol = flattice.lattice.tol
        is_in_fl = all(trial_point_fl.coords .> -tol) && all(trial_point_fl.coords .< 1 - tol) # include (0, 0) but not torus vectors themselves
        if is_in_fl
            push!(bravais_coords, trial_point)
        end
    end

    # take care of sorting 
    sort!(bravais_coords; lt = (v_lat1, v_lat2) -> flattice.order(v_lat1.coords, v_lat2.coords))

    return bravais_coords
end

"""
    atoms(flattice::FiniteLattice)

    Computes all atom coordinates inside the finite lattice in Euclidean coordinates.
    Coordinates are ordered such that identical atoms in different unit cells appear
    next to each other, while the order of unit cells follows flattice.order.
"""
function atoms(flattice::FiniteLattice) :: Vector{EuclideanVector}
    bravais_coords = bravais_cells(flattice) # respects flattice.order already!
    pos = positions(flattice.lattice) # Vector{LatticeVector} of atom positions in (0, 0) cell
 
    # generate list of all coordiantes (same atoms appearing next to each other)
    lattice_coords = Vector{LatticeVector}()
    for position_lat_vec in pos
        for bravais_coord in bravais_coords
            push!(lattice_coords, position_lat_vec + bravais_coord) # sum of `LatticeVector` objects
        end
    end

    # convert to Euclidean coordinates
    return [to_euclidean_basis(v) for v in lattice_coords]
end

"""
    boundary_vecs(flattice::FiniteLattice)

    Returns the vectors defining the boundary as Vector{LatticeVector}
"""
function boundary_vecs(flattice::FiniteLattice) :: Vector{LatticeVector}
    return [LatticeVector(flattice.lattice, flattice.boundary[i]) for i in 1:dim(flattice)]
end


#=


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

=#