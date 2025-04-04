using TOML

@doc raw"""
    Op
    
Defines the information regarding an operator ``\mathbf{O}``.

# Arguments
- `type::AbstractString`: String specifing the type of operator: "HB", "Cup", etc ...
- `coupling::Union{AbstractString,Number}`: String or number specifying the coupling constant;
- `sites::Vector{Int64}`: Vector specifing the sites in each the operator acts on the Hilbert Space;
"""
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


@doc raw"""
    OpSum

Constructor to collect all operators that appear in the Hamiltonian.

One can use the `+=` operator to add operators to the OpSum:

```julia
    opsum = OpSum()
    opsum += Op("HB", "Jd", [1,2])
    opsum += Op("HB", "Jd", [2,3])
    opsum += Op("HB", "Jd", [3,4])
```
"""
mutable struct OpSum
    ops::Vector{Op}
end


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


@doc raw"""
Returns the OpSum() for corresponding to a two-body operator acting between two neighboring sites.

# Arguments
- `type::AbstractString`: String specifing the type of operator: "HB", "Cup", etc ...
- `coupling::AbstractString`: String or number specifying the coupling constant;
- `lattice::FiniteLattice`: The lattice on which the operator acts;
- `num_distance::Integer=1`: at which distance neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.
"""
function neighbor_bonds(type::AbstractString, coupling::AbstractString,
    lattice::FiniteLattice; num_distance::Integer=1)
    nbors = neighbors(lattice; num_distance)
    # Create an OpSum:
    OpSumngb = OpSum()
    for nbor in eachcol(nbors)
        OpSumngb += Op(type, coupling, [nbor[1], nbor[2]])
    end
    unique_ops!(OpSumngb)

    return OpSumngb
end

@doc raw"""
Returns the OpSum() for the corresponding two-body operator acting on the sites specified by the vectors b1 and b2.

    The vectors b1 and b2 are given in the following way:
    - b_[1] = index of the atom in the unit cell (using the same order as given to the lattice instance)
    - b_[2] = Coordinate of the unit cell along the first Bravais vector
    - b_[3] = Coordinate of the unit cell along the second Bravais vector
    - b_[4] = Coordinate of the unit cell along the third Bravais vector
    - ...

# Arguments
- `type::AbstractString`: String specifing the type of operator: "HB", "Cup", etc ...
- `coupling::AbstractString`: String or number specifying the coupling constant;
- `lattice::FiniteLattice`: The lattice on which the operator acts;
- `b1::Vector{<:Integer},`: Vector of integers specifying the first site;
- `b2::Vector{<:Integer},`: Vector of integers specifying the second site;
"""
function lattice_bonds(type::AbstractString, coupling::AbstractString, Flattice::FiniteLattice, b1::Vector{<:Integer}, b2::Vector{<:Integer})

    bravais_coords = bravais_coordinates(Flattice)
    coords = coordinates(Flattice)
    period = Flattice.periodicity
    periodicvec = periodicity_vectors(Flattice)

    if period == [0 0]
        metric = EuclideanMetric()
    else
        metric = PeriodicMetric(periodicvec)
    end

    dist(x, y) = metric(x, y)

    dim = dimension(Flattice)

    if length(b1) != dim + 1
        error("Argument b1 has wrong length. Must be dim( of lattice) + 1!")
    elseif length(b2) != dim + 1
        error("Argument b2 has wrong length. Must be dim( of lattice) + 1!")
    end

    lattice = Flattice.lattice

    OpSumbonds = OpSum()
    for bcoord in eachcol(bravais_coords)

        p1 = bcoord + lattice.positions[:, b1[1]] + lattice.vectors * b1[2:end]
        p2 = bcoord + lattice.positions[:, b2[1]] + lattice.vectors * b2[2:end]


        # Check whether the points are in lattice
        s1 = 0
        for (idx, c) in enumerate(eachcol(coords))
            if dist(p1, c) < 1e-6
                s1 = idx
                break
            end
        end

        s2 = 0
        for (idx, c) in enumerate(eachcol(coords))
            if dist(p2, c) < 1e-6
                s2 = idx
                break
            end
        end

        # append if both are found
        if s1 != 0 && s2 != 0
            ss1 = min(s1, s2)
            ss2 = max(s1, s2)
            OpSumbonds += Op(type, coupling, [ss1, ss2])
        end
    end

    unique_ops!(OpSumbonds)
    return OpSumbonds
end


@doc raw"""
Writes into a TOML file all the operators in the OpSum.
    The file is written in the following format:

```TOML
    Interactions = [
        ["HB", "J1", [1, 2]],
        ["HB", "J2", [2, 3]],
        ...
    ]
```
            
# Arguments
- `opsum::OpSum`: OpSum object containing the operators;
- `filename::String`: Name of the file to write to;
# Keyword arguments
- `index_zero::Bool=false`: If true, the indices of the sites are written starting from 0 instead of 1;
"""
function write_opsum_to_toml!(opsum::OpSum, filename::String; index_zero::Bool=false)

    open(filename, "w") do io
        println(io, "Interactions = [")
        for (i, op) in enumerate(opsum.ops)
            if index_zero == true
                @show op.sites .- 1,
                println(io, "[", " \"$(op.type)\"  ,  ", " \"$(op.coupling)\"  ,  ", op.sites .- 1, " ],")
            else
                println(io, "[", " \"$(op.type)\"  ,  ", " \"$(op.coupling)\"  ,  ", op.sites, " ],")
            end
        end
        println(io, "]")
    end
end