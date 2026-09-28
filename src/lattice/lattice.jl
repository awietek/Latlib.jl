using Printf
using LinearAlgebra
using Base


@doc raw"""
    Lattice(A::Matrix{Float64}, positions::Matrix{Float64}; types=ones(Int64, P), tol=1e-8)
    Lattice(A::Matrix{Float64}, positions::Vector{EuclideanVector}; types, tol)
    Lattice(vs::Vector{EuclideanVector}, positions; types, tol)
    Lattice(A::Matrix{Float64})

An infinite lattice, i.e., a Bravais lattice with a basis of atoms.

A lattice is defined by its:
1.  Bravais lattice vectors, 
2.  atom positions,
3.  and optionally the types of atoms at the atom positions.

The Bravais lattice is defined at a set of points

``\mathbf{R} = \sum_{i=1}^D n_i \mathbf{a}_i``

where ``D`` denotes the dimension of the lattice, 
``n_i`` are integers and ``\mathbf{a}_i`` are the
lattice vectors. For example, in two dimensions,

``\mathbf{R} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2``

with 

``\mathbf{a}_1 = \begin{pmatrix} a_1^x \\ a_1^y \end{pmatrix}``
and 
``\mathbf{a}_2 = \begin{pmatrix} a_2^x \\ a_2^y \end{pmatrix}``

When defining the lattice, the vectors ``\mathbf{a}_i`` are
represented via the `EuclideanVector` type, i.e., real-space
Euclidean vectors relative to the standard basis.
Lattice vectors are conveniently stored in a matrix, where
each row corresponds to a lattice vector:

``\mathbf{A} = \begin{pmatrix} \mathbf{a}_1 \\ \cdots \\ \mathbf{a}_D \end{pmatrix}``.

The atom positions can be given in Cartesian standard coordiantes,
i.e., as a vector of `EuclideanVector`s, or as a simple vector of
Float vectors in which case the positions are interpreted as being
expressed in the basis of the lattice vectors.
The real space coordinates are then given by

``\mathbf{x} = \mathbf{A}^\top \mathcal{X}``,

where ``\mathcal{X}`` is the vector in the basis of lattice vectors.

Each atom position can be given a type to further model symmetries of the lattice.
By default, the `types` parameter defaults to the vector of ones, i.e., all positions
are occupied by identical atoms.

# Arguments
- `A::Matrix{Float64}`: ``D \times D`` matrix whose rows are the Bravais lattice vectors.
- `positions::Matrix{Float64}`: ``P \times D`` matrix whose rows are the atom positions in fractional (lattice basis) coordinates. If given as a `Vector{EuclideanVector}` instead, the positions are interpreted as Cartesian coordinates and converted to the lattice basis.
- `types::Vector{Int64}`: ``P`` dimensional integer vector defining the types of the atoms. Defaults to vector of ones.
- `tol::Float64`: Tolerance for checking whether a point is part of the lattice. Defaults to `1e-8`.

If only `A` is given, a single atom at the origin is assumed.

# Examples
```julia
# square lattice with one atom per unit cell
sq = Lattice([1.0 0.0; 0.0 1.0])

# honeycomb lattice with two atoms per unit cell (positions in lattice basis)
A = [cos(pi/6) sin(pi/6); cos(pi/6) -sin(pi/6)]
hc = Lattice(A, [0.0 0.0; 1/3 1/3])
```

See also [`FiniteLattice`](@ref), [`LatticeVector`](@ref), [`EuclideanVector`](@ref).
"""
struct Lattice
    A::Matrix{Float64}
    positions::Matrix{Float64}
    types::Vector{Int64}
    dim::Int64
    natoms::Int64
    tol::Float64

    # base constructor
    function Lattice(A::Matrix{Float64}, positions::Matrix{Float64}; types::Vector{Int64}=ones(Int64, size(positions)[1]), tol::Float64=1e-8)
        # check for D x D matrix
        dim = size(A)[1]
        if dim != size(A)[2]
            error("Invalid Bravais lattice vectors! The number of rows and columns of the matrix must be the same, i.e., D lattice vectors in D dimensions.")
        end

        # check for linear independence of lattice vectors
        if isapprox(det(A), 0.0, atol=1e-10) 
            error("Invalid Bravais lattice vectors! The lattice vectors must be linearly independent.")
        end

        # check if positions are D-dimensional
        if dim != size(positions)[2]
            error(@sprintf "Dimension of the atom positions (%s) does not match the dimension of the lattice vectors (%s)" size(A) size(positions))
        end

        # check number of atoms = P
        natoms = size(positions)[1]
        if natoms == 0
            error("At least one atom position must be specified.")
        end

        # check if there are repeating positions
        for i in 1:natoms
            for j in i+1:natoms
                if isapprox(positions[i, :], positions[j, :], atol=tol)
                    error(@sprintf "Atom positions %d and %d are approximately the same! Please remove duplicates." i j)
                end
            end
        end

        # check if types are specified consistently (length P vector)
        if length(types) != natoms
            error(@sprintf "Length of 'types' vector (%d) does not match the number of atom positions (%d)." length(types) natoms)
        end

        new(A, positions, types, dim, natoms, tol)
    end

    # ----- Definition of alternative constructors for convenience -----

    # RULE: If positions are given as EuclideanVectors, they are interpreted
    #       as Cartesian coordinates, i.e., not in terms of the lattice basis!

    # alternative constructor using `EuclideanVector` type for lattice vectors but not for positions
    function Lattice(vs::Vector{EuclideanVector}, positions::Matrix{Float64}; types::Vector{Int64}=ones(Int64, size(positions)[1]), tol::Float64=1e-8)
        A = Matrix(hcat([v.coords for v in vs]...)')
        Lattice(A, positions; types = types, tol = tol)
    end

    # alternative constructor using `EuclideanVector` type for lattice vectors and for positions
    function Lattice(vs::Vector{EuclideanVector}, positions::Vector{EuclideanVector}; types::Vector{Int64}=ones(Int64, length(positions)), tol::Float64=1e-8)
        A = Matrix(hcat([v.coords for v in vs]...)')
        pos_vecs = [inv(A') * p.coords for p in positions] # converts Cartesian basis to lattice basis
        pos = Matrix(hcat(pos_vecs...)')
        Lattice(A, pos; types = types, tol = tol)
    end

    # alternative constructor using `EuclideanVector` type for positions but not for lattice vectors
    function Lattice(A::Matrix{Float64}, positions::Vector{EuclideanVector}; types::Vector{Int64}=ones(Int64, length(positions)), tol::Float64=1e-8)
        pos_vecs = [inv(A') * p.coords for p in positions] # converts Cartesian basis to lattice basis
        pos = Matrix(hcat(pos_vecs...)')
        Lattice(A, pos; types = types, tol = tol)
    end

    # alternative constructor without atom positions
    function Lattice(A::Matrix{Float64})
        positions = zeros(Float64, 1, size(A)[2])
        Lattice(A, positions)
    end

    # alternative constructor without atom positions and using `EuclideanVector` type for lattice vectors
    function Lattice(vs::Vector{EuclideanVector})
        A = Matrix(hcat([v.coords for v in vs]...)')
        positions = zeros(Float64, 1, size(A)[2])
        Lattice(A, positions)
    end

end

"""
    dim(lattice::Lattice)
    dim(flattice::FiniteLattice)

Obtain the dimension of the lattice.
"""
dim(lattice::Lattice) = lattice.dim

"""
    natoms(lattice::Lattice)
    natoms(flattice::FiniteLattice)

Obtain the number of atoms in the unit cell of the lattice.
"""
natoms(lattice::Lattice) = lattice.natoms

"""
    positions(lattice::Lattice)
    positions(flattice::FiniteLattice)

Obtain the positions of the atoms in the unit cell at the origin as a `Vector{LatticeVector}`
(i.e., in the basis of the lattice vectors, not in Cartesian coordinates).
"""
function positions(lattice::Lattice) :: Vector{LatticeVector}
    return [LatticeVector(lattice, lattice.positions[i, :]) for i in 1:lattice.natoms]
    # no transformation to lattice basis needed since positions::Matrix{Float64} already assumes lattice basis implicitly
end

@doc """
    get_position(lattice::Lattice, atom_index::Int64)

Obtain the position of a specific atom as a `LatticeVector` object.
"""
function get_position(lattice::Lattice, atom_index::Int64) :: LatticeVector
    if atom_index < 1 || atom_index > lattice.natoms
        error(@sprintf "Specified atom index %s is out of bounds. Must be between 1 and %s according to specified `Lattice`." atom_index lattice.natoms)
    end
    return LatticeVector(lattice, lattice.positions[atom_index, :])
end

@doc """
    lattice_vecs(lattice::Lattice) :: Vector{EuclideanVector}

Returns list of lattice vectors as `EuclideanVector`s.
"""
function lattice_vecs(lattice::Lattice) :: Vector{EuclideanVector}
    return [EuclideanVector(lattice.A[i, :]) for i in 1:lattice.dim]
end


# print lattice information in a nice format
function Base.show(io::IO, lattice::Lattice)
    println(io, "---- Lattice ----")
    println(io, @sprintf "dim = %d" lattice.dim)
    println(io, "Direct lattice vectors:")
    for i in 1:lattice.dim
        println(io, @sprintf "a%d = %s" i lattice.A[i, :])
    end
    println(io, "Atoms:")
    for i in 1:lattice.natoms
        println(io, @sprintf "x%d = %s (type: %d)" i lattice.positions[i, :] lattice.types[i])
    end
end

# function to check if two lattices are equal
# TO-DO: we may want to take into account lattice.tol or atoms matching up to permutations in the future?
function Base.:(==)(l1::Lattice, l2::Lattice)
    # check if dimensions match
    if l1.dim != l2.dim
        return false
    end

    # check if lattice vectors match
    A1 = l1.A
    A2 = l2.A
    if !isapprox(A1, A2)
        return false
    end

    # check if atom positions and types match (up to permutation)
    if l1.natoms != l2.natoms
        return false
    end
    return isapprox(l1.positions, l2.positions) && isapprox(l1.types, l2.types)
end











@doc raw"""
    LatticeVector(lattice::Lattice, coords::Vector{<:Real})

A vector expressed in terms of the basis of Bravais lattice vectors of a [`Lattice`](@ref),
i.e., ``\mathbf{x} = \sum_i c_i \mathbf{a}_i`` with coordinates ``c_i``.

`LatticeVector`s of the same lattice support addition and subtraction.
Use [`to_euclidean_basis`](@ref) to convert to Cartesian coordinates and
[`in_lattice`](@ref) to check whether the vector is a Bravais lattice vector
(i.e., has integer coordinates).

# Fields
- `lattice::Lattice`: the lattice in whose basis the vector is expressed.
- `coords::Vector{Float64}`: the coordinates of the vector in the lattice basis.
- `dim::Int`: dimension of the vector.

# Examples
```julia
t1 = LatticeVector(honeycomb, [2, -1])
to_euclidean_basis(t1)   # EuclideanVector([0.866..., 1.5])
in_lattice(t1)           # true
```
"""
struct LatticeVector
    lattice::Lattice
    coords::Vector{Float64}
    dim::Int

    # standard constructor with float inputs
    function LatticeVector(lattice::Lattice, coords::Vector{Float64})
        dim = length(coords)
        if !(dim == lattice.dim)
            throw(ArgumentError("Dimension coordinate vector in lattice basis must match the dimension of the lattice."))
        end
        new(lattice, coords, dim)
    end

    # alternative constructor with integer inputs for convenience
    function LatticeVector(lattice::Lattice, coords::Vector{Int})
        LatticeVector(lattice, float.(coords))
    end

end

# equality check for LatticeVectors
function Base.:(==)(v1::LatticeVector, v2::LatticeVector)
    v1.lattice == v2.lattice && v1.dim == v2.dim && isapprox(v1.coords, v2.coords)
end

# addition for LatticeVectors
function Base.:+(v1::LatticeVector, v2::LatticeVector)
    if v1.lattice != v2.lattice
        throw(ArgumentError("Cannot add LatticeVectors from different lattices."))
    end
    if v1.dim != v2.dim
        throw(ArgumentError("Cannot add LatticeVectors of different dimensions."))
    end
    return LatticeVector(v1.lattice, v1.coords + v2.coords)
end

# subtraction for LatticeVectors
function Base.:-(v1::LatticeVector, v2::LatticeVector)
    return v1 + LatticeVector(v1.lattice, -v2.coords)
end

@doc raw"""
    in_lattice(v::LatticeVector)
    in_lattice(lattice::Lattice, v::EuclideanVector)
    in_lattice(v::FiniteLatticeVector)

Returns whether a vector is a Bravais lattice vector of the underlying lattice,
i.e., whether its coordinates in the lattice basis are integers (up to the tolerance of the lattice).
"""
function in_lattice(v_lat::LatticeVector)
    v_lat_diff = v_lat.coords - Base.round.(v_lat.coords)
    return all(is_whole.(v_lat_diff; atol=v_lat.lattice.tol))
end

@doc raw"""
    to_lattice_basis(lattice::Lattice, v::EuclideanVector) -> LatticeVector
    to_lattice_basis(v::FiniteLatticeVector) -> LatticeVector

Convert a vector in Cartesian coordinates (or in the basis of boundary vectors of a
finite lattice) to a [`LatticeVector`](@ref) expressed in the basis of the lattice vectors.
"""
function to_lattice_basis(lattice::Lattice, v::EuclideanVector) :: LatticeVector
    v_lat = inv(lattice.A') * v.coords
    return LatticeVector(lattice, v_lat)
end


@doc raw"""
    to_euclidean_basis(v::LatticeVector) -> EuclideanVector
    to_euclidean_basis(v::FiniteLatticeVector) -> EuclideanVector

Convert a vector expressed in the lattice basis (or in the basis of boundary vectors of a
finite lattice) to an [`EuclideanVector`](@ref) in Cartesian coordinates.
"""
function to_euclidean_basis(v_lat::LatticeVector) :: EuclideanVector
    return EuclideanVector(v_lat.lattice.A' * v_lat.coords)
end


function in_lattice(lattice::Lattice, v::EuclideanVector)
    v_lat = to_lattice_basis(lattice, v)
    return in_lattice(v_lat)
end







#=

    

"""
    vector_position(lattice::Lattice, point::AbstractVector)

Returns which multiple of the spanning vectors and which position of
the lattice a point is represented by, if it is part of the lattice
"""
function vector_position(lattice::Lattice, point::AbstractVector)
    r = lattice.vectors \ point

    vector = floor.(Int, round.(r, digits=8))
    
    r = r - vector
    point0 = lattice.vectors * r
    positions_cartesian = lattice.vectors *lattice.positions

    for (idx, a) in enumerate(eachcol(positions_cartesian))
        if isapprox(point0, a, rtol=1e-8, atol=1e-8)
            return vector, idx
        end
    end

    for (idx, a) in enumerate(eachcol(positions_cartesian))
        @show point0, a
    end
    exit()
    return [], 0
end

=#


