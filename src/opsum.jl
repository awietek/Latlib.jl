struct Op
    type::AbstractString
    coupling::Union{AbstractString, Number}
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

# Initializer
# function OpSum(Opsv::Vector{Op})
#     OpSum(Opsv)
# end

function OpSum()
    OpSum(Op[])
end

#Overload the += operator for OpSum
    function Base.:(+)(sumOp::OpSum, op::Op)
    return OpSum(push!(sumOp.ops, op))
end

# function nearest_neighbor_bonds(type::AbstractString, coupling::AbstractString,
#                                 lattice::Lattice, num_distance::Integer;
#                                 periodic_dims::Vector{<:Integer}=Integer[],
#                                 zero_indexed::Bool=false)
#     coords = coordinates(lattice)
#     period = periodicity(lattice, periodic_dims)
#     nbors = neighbors(coords, num_distance; periodicity=period)

#     # Create vector of bonds
#     bonds = Bond[]
#     for (site, nbor) in enumerate(nbors)
#         for nsite in nbor
#             s1 = min(site, nsite)
#             s2 = max(site, nsite)
#             if zero_indexed
#                 push!(bonds, Bond(type, coupling, [s1-1, s2-1]))
#             else
#                 push!(bonds, Bond(type, coupling, [s1, s2]))
#             end                
#         end
#     end
#     return sort(unique(bonds))
# end


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
