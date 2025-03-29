struct Op
    type::AbstractString
    coupling::Union{AbstractString,Number}
    sites::Vector{Int64}
end

# Methods of for the struct Op:
function Base.:(==)(b1::Op, b2::Op)
    return b1.type == b2.type && b1.coupling == b2.coupling && b1.sites == b2.sites
end

function Base.isequal(b1::Op, b2::Op)
    return b1 == b2
end

function Base.hash(b::Op, h::UInt64)
    return hash(b.type, h) + hash(b.coupling, h) + hash(b.sites, h)
end

function Base.isless(b1::Op, b2::Op)
    if b1.type < b2.type
        return true
    elseif b1.coupling < b2.coupling
        return true
    else
        s1 = b1.sites
        s2 = b2.sites
        sort!(s1)
        sort!(s2)
        if s1 < s2
            return true
        else
            return false
        end
    end
end

# Define a struct for collecting op
mutable struct OpSum
    ops::Vector{Op}
end

#Initialize OpSum
function OpSum()
    OpSum(Op[])
end

#Overload the += operator for OpSum
function Base.:(+)(sumOp::OpSum, op::Op)
    return OpSum(push!(sumOp.ops, op))
end

function Base.:(+)(sumOp1::OpSum, sumOp2::OpSum)
    return OpSum(push!(sumOp1.ops, sumOp2.ops...))
end

function unique_ops!(opsum::OpSum)
    opsum.ops = sort(unique(opsum.ops))
end


function neighbor_bonds(type::AbstractString, coupling::AbstractString,
    lattice::FiniteLattice; num_distance::Integer=1)
    nbors = neighbors(lattice; num_distance)

    # Create an OpSum:
    OpSumngb = OpSum()
    for nbor in eachcol(nbors)
        OpSumngb += Op(type, coupling, [nbor[1],nbor[2]])
    end
    return unique_ops!(OpSumngb)
end




# function lattice_bonds(type::AbstractString, coupling::AbstractString,
#                        lattice::Lattice, b1::Vector{<:Integer}, b2::Vector{<:Integer};
#                        periodic_dims::Vector{<:Integer}=Integer[],
#                        zero_indexed::Bool=false)
#     bravais_coords = bravais_coordinates(lattice)
#     coords = coordinates(lattice)
#     period = periodicity(lattice, periodic_dims)
#     if period == nothing
#         metric = EuclideanMetric()
#     else
#         metric = PeriodicMetric(period)
#     end
#     dist(x, y) = metric(x, y)

#     dim = dimension(lattice)
#     if length(b1) != dim + 1
#         error("Argument b1 has wrong length. Must be dim( of lattice) + 1!")
#     elseif length(b2) != dim + 1
#         error("Argument b2 has wrong length. Must be dim( of lattice) + 1!")
#     end

#     bonds = Bond[]
#     for bcoord in eachcol(bravais_coords)
#         p1 = bcoord + lattice.basis[:, b1[1]] + lattice.lattice * b1[2:end]
#         p2 = bcoord + lattice.basis[:, b2[1]] + lattice.lattice * b2[2:end]

#         # Check whether the points are in lattice
#         s1 = 0
#         for (idx, c) in enumerate(eachcol(coords))
#             if dist(p1, c) < 1e-6
#                 s1 = idx
#                 break
#             end
#         end

#         s2 = 0
#         for (idx, c) in enumerate(eachcol(coords))
#             if dist(p2, c) < 1e-6
#                 s2 = idx
#                 break
#             end
#         end

#         # append if both are found
#         if s1 != 0 && s2 != 0
#             ss1 = min(s1, s2)
#             ss2 = max(s1, s2)

#             if zero_indexed
#                 push!(bonds, Bond(type, coupling, [ss1-1, ss2-1]))
#             else
#                 push!(bonds, Bond(type, coupling, [ss1, ss2]))
#             end
#         end
#     end
#     return sort(bonds)
# end
