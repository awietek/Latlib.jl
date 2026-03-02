#=
    Latlib defines the following types of vectors:
        1. 'EuclideanVector' for vectors in Euclidean (real) space, relative to standard basis.
        2. 'LatticeVector' for vectors expressed in terms of the basis of a lattice.

=#

@doc raw"""
    `EuclideanVector` is a type representing vectors in Euclidean (real) space, relative to the standard basis. It consists of two fields:
    
    - `coords::Vector{Float64}`: A vector containing the coordinates of the Euclidean vector.
    - `dim::Int`: Dimension of the vector.
"""
struct EuclideanVector
    coords::Vector{Float64}
    dim::Int

    # standard constructor with float inputs
    function EuclideanVector(coords::Vector{Float64})
        dim = length(coords)
        if !(dim in (2, 3))
            throw(ArgumentError("EuclideanVector only supports 2D and 3D vectors."))
        end
        new(coords, dim)
    end

    # alternative constructor with integer inputs for convenience
    function EuclideanVector(coords::Vector{Int})
        EuclideanVector(float.(coords))
    end

end 

# equality check
function Base.:(==)(v1::EuclideanVector, v2::EuclideanVector)
    v1.dim == v2.dim && isapprox(v1.coords, v2.coords)
end

# addition
function Base.:+(v1::EuclideanVector, v2::EuclideanVector)
    if v1.dim != v2.dim
        throw(ArgumentError("Cannot add EuclideanVectors of different dimensions."))
    end
    return EuclideanVector(v1.coords + v2.coords)
end

# subtraction
function Base.:-(v1::EuclideanVector, v2::EuclideanVector)
    return v1 + EuclideanVector(-v2.coords)
end










