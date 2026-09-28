#=
    Latlib defines the following types of vectors:
        1. 'EuclideanVector' for vectors in Euclidean (real) space, relative to standard basis.
        2. 'LatticeVector' for vectors expressed in terms of the basis of a lattice.

=#

@doc raw"""
    EuclideanVector(coords::Vector{<:Real})

A vector in Euclidean (real) space, given relative to the standard Cartesian basis.
Only two- and three-dimensional vectors are supported.

`EuclideanVector`s support addition, subtraction, multiplication by a scalar, the
inner product `v1 * v2`, and `LinearAlgebra.norm`.

# Fields
- `coords::Vector{Float64}`: coordinates of the vector.
- `dim::Int`: dimension of the vector (2 or 3).

# Examples
```julia
v = EuclideanVector([1.0, 2.0])
w = EuclideanVector([3, 4])   # integer input is converted to Float64
v + w                         # EuclideanVector([4.0, 6.0])
2 * v                         # EuclideanVector([2.0, 4.0])
v * w                         # 11.0 (inner product)
```
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

# product with scalar
function Base.:*(scalar::Number, v::EuclideanVector)
    return EuclideanVector(scalar * v.coords)
end

# inner product
function Base.:*(v1::EuclideanVector, v2::EuclideanVector)
    if v1.dim != v2.dim
        throw(ArgumentError("Cannot compute inner product of EuclideanVectors of different dimensions."))
    end
    return dot(v1.coords, v2.coords)
end

# norm (extends LinearAlgebra.norm)
function LinearAlgebra.norm(v::EuclideanVector)
    return sqrt(v * v)
end

# print
function Base.show(io::IO, v::EuclideanVector)
    print(io, "EuclideanVector(", v.coords, ")")
end







