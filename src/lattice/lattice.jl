using Printf

@doc raw"""
    Lattice

A lattice is defined by its Bravais lattice vectors, 
atomic positions and optionally the type of 
the atomic lattice sites.

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

Here, the lattice vectors are given in Cartesian
coordinates.

The atomic positions are stored in fractional values 
``\mathcal{X}`` relative to lattice vectors ``\mathbf{a}_i``. 
This means, if we denote the matrix of lattice vectors by

``\mathbf{A} = \begin{pmatrix} \mathbf{a}_1 | \cdots | \mathbf{a}_D \end{pmatrix}``

the cartesian coordinates ``\mathbf{x}`` of an atomic position 
are given by

``\mathbf{x} = \mathbf{A} \mathcal{X}``

Each atomic position is also given a type, which is simply an
integer number. This can be useful if there are different 
atomic species in a unit cell.

# Arguments
- `vectors::AbstractMatrix`: ``D \times D`` matrix whose columns are the Bravais lattice vectors.
- `positions::AbstractMatrix`: ``D \times P`` matrix whose columns are the atomic positions.
- `types::AbstractVector`: ``P`` dimensional integer vector defining the types of the atoms.
"""
struct Lattice
    vectors::Matrix{Float64}
    positions::Matrix{Float64}
    types::Vector{Int64}

    function Lattice(vectors::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector)
        dim = size(vectors)[1]
        if dim != size(vectors)[2]
            error("Invalid Bravais vectors specified. \"vectors\" must be a square matrix")
        end

        if dim != size(positions)[1]
            error(@sprintf "Incompatible dimension of Bravais \"vectors\" %s and atomic \"positions\" %s" size(vectors) size(positions))
        end

        nbasis = size(positions)[2]
        if length(types) == 0
            types = ones(Integer, nbasis)
        end

        if length(types) != nbasis
            error("Number of \"types\" must be same as number of \"positions\"")
        end
        new(vectors, positions, types)
    end
end

@doc raw""" 
    Lattice(vectors::AbstractMatrix, positions::AbstractMatrix)

Create a lattice from Bravais lattice vectors and atomic positions.

The types of the atomic positions are automatically set to \"1\".

# Arguments
- `vectors::AbstractMatrix`: ``D \times D`` matrix whose columns are the Bravais lattice vectors.
- `positions::AbstractMatrix`: ``P \times D`` matrix whose rows are the atomic positions
"""
Lattice(vectors::AbstractMatrix, positions::AbstractMatrix) = Lattice(vectors, positions, Integer[])

@doc raw""" 
    Lattice(vectors::AbstractMatrix)

Create a Bravais lattice from the Bravais lattice vectors. 

The (single) atomic position is set to be zero and of type \"1\".

# Argument
- `vectors::AbstractMatrix`: ``D \times D`` matrix whose columns are the Bravais lattice vectors.
"""
Lattice(vectors::AbstractMatrix) = Lattice(vectors, zeros(Real, 1, size(vectors)[1]))

"""
    dimension(lattice::Lattice)

Obtain the dimension of the lattice.
"""
dimension(lattice::Lattice) = size(lattice.vectors)[1]

"""
    natoms(lattice::Lattice)

Obtain the number of atomic positions
"""
natoms(lattice::Lattice) = size(lattice.positions)[2]

"""
    inlattice(lattice::Lattice, point::AbstractVector)

Returns whether a given point is part of the lattice
"""
function inlattice(lattice::Lattice, point::AbstractVector)
    r = lattice.vectors \ point
    r = r - floor.(Int, round.(r, digits=8))
    point0 = lattice.vectors * r
    positions_cartesian = lattice.vectors * lattice.positions
    for a in eachcol(positions_cartesian)
        if isapprox(point0, a, rtol=1e-8, atol=1e-8)
            return true
        end
    end
    return false
end

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

function Base.show(io::IO, lattice::Lattice)
    dim = dimension(lattice)
    println(io, "Lattice")
    println(io, @sprintf "D = %d" dim)
    for i in 1:dim
        println(io, @sprintf "a%d = %s" i lattice.vectors[:,i])
    end
    for i in 1:natoms(lattice)
        println(io, @sprintf "x%d = %s (type: %d)" i lattice.positions[:,i] lattice.types[i])
    end
end
