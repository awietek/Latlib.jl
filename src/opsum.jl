using TOML

@doc raw"""
    Op(type::String, cpl::Union{String, Number}, sites::Vector{Int64})

Represents a single term of a many-body Hamiltonian, e.g. a two-body interaction between
two sites.

# Arguments
- `type::String`: string specifying the type of operator, e.g. "SdotS", "HB", "Cup", ...
- `cpl::Union{String, Number}`: string or number specifying the coupling constant, e.g. "J" or 1.0
- `sites::Vector{Int64}`: 1-based indices of the sites on which the operator acts, as enumerated by [`atoms`](@ref).

```julia
Op("SdotS", "J", [1, 2])
```
"""
struct Op
    type::String
    cpl::Union{String, Number}
    sites::Vector{Int64}

    # standard constructor
    function Op(type::String, cpl::Union{String, Number}, sites::Vector{Int64})
        # check that sites are not empty
        if length(sites) == 0
            error("Vector of sites specified in Op() is empty.")
        end

        # check that sites does not contain duplicates
        if length(unique(sites)) != length(sites)
            error("The vector of sites in the constructor of Op must not contain duplicates.")
        end

        # check that sites are positive integers
        for s in sites
            if s < 1
                error(@sprintf "The vector of sites in the constructor of Op must contain positive integers (got %s)" s)
            end
        end

        new(type, cpl, sites)
    end
end

# Methods of for the struct Op:
function Base.:(==)(b1::Op, b2::Op)
    return b1.type == b2.type && b1.cpl == b2.cpl && b1.sites == b2.sites
end

function Base.isequal(b1::Op, b2::Op)
    return b1 == b2
end

function Base.hash(b::Op, h::UInt64)
    return hash(b.type, h) + hash(b.cpl, h) + hash(b.sites, h)
end

function Base.isless(b1::Op, b2::Op)
    if b1.type < b2.type
        return true
    elseif b1.cpl < b2.cpl
        return true
    else
        s1 = b1.sites
        s2 = b2.sites
        s1_sort = sort(s1)
        s2_sort = sort(s2)
        if s1_sort < s2_sort
            return true
        else
            return false
        end
    end
end


@doc raw"""
    OpSum()
    OpSum(ops::Vector{Op})

A sum of operators [`Op`](@ref), representing a Hamiltonian.

Operators and other `OpSum`s can be added with `+=`:

```julia
opsum = OpSum()
opsum += Op("HB", "Jd", [1, 2])
opsum += Op("HB", "Jd", [2, 3])
opsum += neighbor_interaction("HB", "J", flattice)
```

Typically, an `OpSum` is generated on a [`FiniteLattice`](@ref) with
[`neighbor_interaction`](@ref) or [`lattice_interaction`](@ref) and written to a file
with [`write_toml`](@ref).
"""
mutable struct OpSum
    ops::Vector{Op}

    # standard constructor
    function OpSum(ops::Vector{Op})
        new(ops)
    end

    # constructor with empty vector of operators
    function OpSum()
        new(Op[])
    end

end

# define OpSum += Op
function Base.:(+)(sum::OpSum, op::Op)
    if length(op.sites) > 0
        return OpSum(push!(sum.ops, op))
    else
        error("Cannot add an operator `Op` with empty vector of sites to `OpSum`.")
    end
end

# define OpSum += OpSum
function Base.:(+)(sum1::OpSum, sum2::OpSum)
    return OpSum(push!(sum1.ops, sum2.ops...))
end

"""
    unique_ops!(opsum::OpSum)

Remove duplicate operators from `opsum` and sort the remaining operators in place.
"""
function unique_ops!(opsum::OpSum)
    opsum.ops = sort(unique(opsum.ops))
end

# equality check for two OpSums
function Base.:(==)(opsum1::OpSum, opsum2::OpSum)
    unique_ops!(opsum1)
    unique_ops!(opsum2)
    return opsum1.ops == opsum2.ops
end


@doc raw"""
    neighbor_interaction(type::String, cpl::Union{String, Number}, flattice::FiniteLattice; num_distance::Int64=1)

Constructs the `OpSum` for a two-body interaction of k-th
nearest neighbors on a finite lattice where k = num_distance.

# Arguments
- `type::String`: String specifing the type of operator: "HB", "Cup", etc ...;
- `cpl::Union{String, Number}`: String or number specifying the coupling constant;
- `flattice::FiniteLattice`: The lattice on which the operator acts;

# Keyword arguments
- `num_distance::Int64=1`: Distance at which neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.

# Examples
```julia
fl = FiniteLattice(square, [4 0; 0 4], true)
H = neighbor_interaction("SdotS", "J1", fl)                  # nearest neighbors
H += neighbor_interaction("SdotS", "J2", fl; num_distance=2)  # next-nearest neighbors
```
"""
function neighbor_interaction(type::String,
                        cpl::Union{String, Number},
                        flattice::FiniteLattice;
                        num_distance::Int64=1)
    if num_distance < 1
        error("Two-body interaction defined by neighbor_interaction must be called with num_distance >= 1, i.e., we do not consider self-interactions.")
    end

    nbs = neighbors(flattice; num_distance) # returns Vector{Tuple{Int64, Int64}} of index pairs (i, j) with i<j
    # Create an OpSum:
    ops = OpSum()
    for (i, j) in nbs
        ops += Op(type, cpl, [i, j])
    end
    unique_ops!(ops)

    return ops
end

@doc raw"""
    lattice_interaction(type::String, cpl::Union{String, Number}, flattice::FiniteLattice, atom1::Int64, atom2::Int64, cell2::Union{Vector{Int64}, LatticeVector})

Returns the OpSum() for the corresponding two-body operator acting
between two sites that may (or may not) be in different Bravais cells.
The interaction is repeated for all Bravais cells of the finite lattice!

# Arguments
- `type::String`: String specifing operator type: "HB", "Cup", etc ...;
- `cpl::Union{String, Number}`: String or number specifying the coupling constant.;
- `flattice::FiniteLattice`: The lattice on which the operator acts.;
- `atom1::Int64`: Index (as defined by `positions(flattice)`) of the first atom taking part in the interaction;
- `atom2::Int64`: Index (as defined by `positions(flattice)`) of the second atom taking part in the interaction;
- `cell2::Union{Vector{Int64}, LatticeVector}`: Bravais cell of atom2 (atom1 is always assumed in the origin Bravais cell).

# Returns
- `OpSum`: contains the corresponding two-body operator acting between the specified sites in all Bravais cells of the finite lattice. The site indices of the operators refer to the enumeration of sites given by [`atoms`](@ref).

# Examples
```julia
# Kitaev model on the honeycomb lattice (two atoms per unit cell)
fl = FiniteLattice(honeycomb, [2 -2; 1 1], true)
H = OpSum()
H += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0])   # atom 1 to atom 2 in the same cell
H += lattice_interaction("SySy", "KY", fl, 1, 2, [0, -1])  # atom 1 to atom 2 in cell [0, -1]
H += lattice_interaction("SzSz", "KZ", fl, 2, 1, [1, 0])   # atom 2 to atom 1 in cell [1, 0]
```
"""
function lattice_interaction(type::String,
                       cpl::Union{String, Number},
                       flattice::FiniteLattice,
                       atom1::Int64,
                       atom2::Int64,
                       cell2::Union{Vector{Int64}, LatticeVector})

    # assume that cell2 is of `LatticeVector` type`
    if isa(cell2, Vector{Int64})
        cell2 = LatticeVector(flattice.lattice, cell2)
    end

    # check if cell2 is a valid Bravais cell
    if !in_lattice(cell2)
        error("Specified value of cell2 is not a valid Bravais vector of the lattice. Please specify a vector of integers of a `LatticeVector` with integer coefficients.")
    end
    lat = flattice.lattice
    if !( cell2.lattice == lat)
        error("Specified `LatticeVector` for cell2 in lattice_bonds() must have the same underlying lattice instance as the finite lattice.")
    end
    lat_dim = dim(flattice)

    # check if atom1 and atom2 are valid
    Natoms = natoms(flattice)
    atom_vec = [atom1, atom2]
    for a in atom_vec
        if a < 1 || a > Natoms
            error(@sprintf "Specified atom index %s is out of bounds. Must be between 1 and %s according to specified `FiniteLattice`." a Natoms)
        end
    end

    # get Bravais cells (as Vector{LatticeVector})
    brav_cells = bravais_cells(flattice)

    # get Euclidean coordinates (as Vector{EuclideanVector}) of all atoms in each Bravais cell 
    atom_coords = atoms(flattice)
    
    # take care of periodicity of finite lattice
    period = periodicity(flattice)

    # define metric and distance function
    metric = nothing
    if all(period .== false)
        metric = EuclideanMetric()
    else
        metric = PeriodicEuclideanMetric(flattice)
    end
    dist(x, y) = metric(x, y)

    # construct OpSum()
    a1_lvec = get_position(flattice, atom1)
    a2_lvec = get_position(flattice, atom2) + cell2
    a1_evec = to_euclidean_basis(a1_lvec)
    a2_evec = to_euclidean_basis(a2_lvec)
    ops = OpSum()
    for cell in brav_cells

        cell_evec = to_euclidean_basis(cell)

        # get coordinates of atom1 and atom2 in terms of lattice basis
        a1_evec_cell = a1_evec + cell_evec
        a2_evec_cell = a2_evec + cell_evec

        # identify these coordinates with atom indices
        site1 = nothing
        for (idx, atom_evec) in enumerate(atoms(flattice)) # gives Vector{EuclideanVector}
            if dist(a1_evec_cell, atom_evec) < flattice.lattice.tol
                site1 = idx
                break
            end
        end
        if isnothing(site1)
            continue
            #error(@sprintf "Could not match atom1 %s in cell %s to any of the vectors in flattice.lattice.positions." a1_lvec cell)
        end

        site2 = nothing
        for (idx, atom_evec) in enumerate(atoms(flattice)) # gives Vector{EuclideanVector}
            if dist(a2_evec_cell, atom_evec) < flattice.lattice.tol
                site2 = idx
                break
            end
        end
        if isnothing(site2)
            continue
            #error(@sprintf "Could not match atom2 %s in cell %s to any of the vectors in flattice.lattice.positions." a2_lvec cell)
        end

        # append tuples (i, j) where i < j
        v1 = min(site1, site2)
        v2 = max(site1, site2)
        ops += Op(type, cpl, [v1, v2])
    end

    return ops
end

