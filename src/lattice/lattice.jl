using IterTools

@doc raw"""
    Lattice

A lattice is defined by its Bravais lattice, 
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
- `lattice::Matrix{Real}`: ``D \times D`` matrix whose columns are the Bravais lattice vectors.
- `position::Matrix{Real}`: ``P \times D`` matrix whose rows are the atomic positions.
- `types::Vector{Integer}`: ``P`` dimensional integer vector defining the types of the atoms.
"""
struct Lattice
    lattice::Matrix{Real}
    position::Matrix{Real}
    types::Vector{Integer}

    function Lattice(lattice, position, types)
        dim = size(lattice)[1]
        if dim != size(lattice)[2]
            error("Invalid Bravais lattice specified. \"lattice\" must be a square matrix")
        end

        if dim != size(position)[2]
            error("Incompatible dimension of lattice and atomic positions")
        end

        nbasis = size(position)[1]
        if size(types) == 0
            types = ones(Integer, nbasis)
        end

        if size(types) != nbasis
            error("Number of \"types\" must be same as number of \"positions\"")
        end
        new(lattice, position, types)
    end
end

@doc raw""" 
    Lattice(lattice::Matrix{Real}, position::Matrix{Real})

Create a lattice from Bravais lattice vectors and atomic positions.

The types of the atomic positions are automatically set to \"1\".

# Arguments
- `lattice::Matrix{Real}`: ``D \times D`` matrix whose columns are the Bravais lattice vectors.
- `position::Matrix{Real}`: ``P \times D`` matrix whose rows are the atomic positions
"""
Lattice(lattice::Matrix{Real}, position::Matrix{Real}) = Lattice(lattice, position, Integer[])

@doc raw""" 
    Lattice(lattice::Matrix{Real})

Create a Bravais lattice from the Bravais lattice vectors. 

The (single) atomic position is set to be zero and of type \"1\".

# Argument
- `lattice::Matrix{Real}`: ``D \times D`` matrix whose columns are the Bravais lattice vectors.
"""
Lattice(lattice::Matrix{Real}) = Lattice(lattice, zeros(Real, 1, size(lattice)[1]))

"""
    dimension(lattice::Lattice)

Obtain the dimension of the lattice.
"""
dimension(lattice::Lattice) = size(lattice.lattice)[1]

"""
    natoms(lattice::Lattice)

Obtain the number of atomic positions
"""
natoms(lattice::Lattice) = size(lattice.position)[1]
