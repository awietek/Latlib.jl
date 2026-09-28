using LinearAlgebra

# Collections of functions used for sorting.
# order(a, b) = true iff `a` should be ordered before `b`.


# --------------------------------------------------------
#             Sorting routines for 2D vectors
# --------------------------------------------------------

"""
    order_xy(a::Vector{Float64}, b::Vector{Float64}) -> Bool

Ordering of 2D vectors by x-coordinate first, then y-coordinate (smaller values first).
Returns `true` if `a` should be ordered before `b`. Can be passed as `bravais_order` or
`atom_order` to [`FiniteLattice`](@ref).
"""
function order_xy(a::Vector{Float64}, b::Vector{Float64})
    if (length(a) != length(b)) || (length(a) != 2)
        error("Vectors a and b must be vectors of length 2 for order_xy!")
    end
    return (a[1] < b[1]) || ((a[1] == b[1]) && (a[2] < b[2]))
end

"""
    order_yx(a::Vector{Float64}, b::Vector{Float64}) -> Bool

Ordering of 2D vectors by y-coordinate first, then x-coordinate (smaller values first).
Returns `true` if `a` should be ordered before `b`. Can be passed as `bravais_order` or
`atom_order` to [`FiniteLattice`](@ref).
"""
function order_yx(a::Vector{Float64}, b::Vector{Float64})
    if (length(a) != length(b)) || (length(a) != 2)
        error("Vectors a and b must be vectors of length 2 for order_yx!")
    end
    return (a[2] < b[2]) || ((a[2] == b[2]) && (a[1] < b[1]))
end

# --------------------------------------------------------
#             Sorting routines for 3D vectors
# --------------------------------------------------------

"""
    order_xyz(a::Vector{Float64}, b::Vector{Float64}) -> Bool

Ordering of 3D vectors by x-coordinate first, then y, then z (smaller values first).
Returns `true` if `a` should be ordered before `b`. Can be passed as `bravais_order` or
`atom_order` to [`FiniteLattice`](@ref).
"""
function order_xyz(a::Vector{Float64}, b::Vector{Float64})
    if (length(a) != length(b)) || (length(a) != 3)
        error("Vectors a and b must be vectors of length 3 for order_xyz!")
    end
    return order_xy(a[1:2], b[1:2]) || ((a[1] == b[1]) && (a[2] == b[2]) && (a[3] < b[3]))
end

