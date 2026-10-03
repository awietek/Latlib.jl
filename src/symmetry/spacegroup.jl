import Spglib
using LinearAlgebra
using Printf

# ----------------------------------------------------------------------
#                        Space-group type tables
# ----------------------------------------------------------------------

# spglib only works in three dimensions. A two-dimensional lattice is therefore
# handed to spglib as a single planar layer (all atoms at z = 0) stacked with a
# large vacuum along z. The space group of such a layer is the plane group times
# {1, σ_h}, and each of the 17 plane groups corresponds to exactly one
# three-dimensional space-group type:
#   3D ITA number => (plane group number, symbol, point group (HM), point group (Schoenflies))
const LAYER_TO_PLANE_GROUP = Dict{Int, Tuple{Int, String, String, String}}(
    6   => (1,  "p1",   "1",   "C1"),
    10  => (2,  "p2",   "2",   "C2"),
    25  => (3,  "pm",   "m",   "Cs"),
    26  => (4,  "pg",   "m",   "Cs"),
    38  => (5,  "cm",   "m",   "Cs"),
    47  => (6,  "p2mm", "2mm", "C2v"),
    51  => (7,  "p2mg", "2mm", "C2v"),
    55  => (8,  "p2gg", "2mm", "C2v"),
    65  => (9,  "c2mm", "2mm", "C2v"),
    83  => (10, "p4",   "4",   "C4"),
    123 => (11, "p4mm", "4mm", "C4v"),
    127 => (12, "p4gm", "4mm", "C4v"),
    174 => (13, "p3",   "3",   "C3"),
    187 => (14, "p3m1", "3m",  "C3v"),
    189 => (15, "p31m", "3m",  "C3v"),
    175 => (16, "p6",   "6",   "C6"),
    191 => (17, "p6mm", "6mm", "C6v"),
)

const SYMMORPHIC_PLANE_GROUPS = Set([1, 2, 3, 5, 6, 9, 10, 11, 13, 14, 15, 16, 17])

const SYMMORPHIC_SPACE_GROUPS = Set([
    1, 2, 3, 5, 6, 8, 10, 12, 16, 21, 22, 23, 25, 35, 38, 42, 44, 47, 65, 69, 71,
    75, 79, 81, 82, 83, 87, 89, 97, 99, 107, 111, 115, 119, 121, 123, 139,
    143, 146, 147, 148, 149, 150, 155, 156, 157, 160, 162, 164, 166,
    168, 174, 175, 177, 183, 187, 189, 191,
    195, 196, 197, 200, 202, 204, 207, 209, 211, 215, 216, 217, 221, 225, 229,
])


# ----------------------------------------------------------------------
#                         Symmetry operations
# ----------------------------------------------------------------------

@doc raw"""
    SymmetryOperation(W::Matrix{Int}, w::Vector{Float64})

A space-group operation acting on coordinates ``\mathcal{X}`` expressed in the basis
of the lattice vectors (see [`LatticeVector`](@ref)),

``\mathcal{X} \mapsto W \mathcal{X} + \mathbf{w}``,

where `W` is an integer matrix (the rotational part) and `w` the translational part.
In Cartesian coordinates ``\mathbf{x} = \mathbf{A}^\top \mathcal{X}`` the rotational
part reads ``\mathbf{A}^\top W \mathbf{A}^{-\top}``, see [`cartesian_rotation`](@ref).

Operations can be applied to coordinate vectors, `op(x)`, and composed, `op1 * op2`
(first `op2`, then `op1`).
"""
struct SymmetryOperation
    W::Matrix{Int}
    w::Vector{Float64}
end

(op::SymmetryOperation)(x::AbstractVector{<:Real}) = op.W * x + op.w

Base.:*(op1::SymmetryOperation, op2::SymmetryOperation) =
    SymmetryOperation(op1.W * op2.W, op1.W * op2.w + op1.w)

function Base.show(io::IO, op::SymmetryOperation)
    print(io, "SymmetryOperation(W = ", op.W, ", w = ", op.w, ")")
end

"""
    cartesian_rotation(op::SymmetryOperation, lattice::Lattice) -> Matrix{Float64}

Rotational part of `op` in Cartesian coordinates, an orthogonal matrix.
"""
function cartesian_rotation(op::SymmetryOperation, lattice::Lattice) :: Matrix{Float64}
    return lattice.A' * op.W / lattice.A'
end

# reduce the entries of a vector to [0, 1), mapping values close to 1 onto 0
function _reduce_mod1(v::AbstractVector{<:Real}; atol::Float64=1e-10) :: Vector{Float64}
    r = v .- floor.(v .+ atol)
    r[abs.(r) .< atol] .= 0.0
    return r
end

# whether two operations agree modulo Bravais translations
_equal_mod_lattice(op1::SymmetryOperation, op2::SymmetryOperation; atol::Float64=1e-8) =
    op1.W == op2.W && all(is_whole.(op1.w - op2.w; atol=atol))

# For an operation mapping the lattice onto itself, find for every atom `a` of the
# unit cell the atom `b = perm[a]` of the same type and the integer vector
# `offsets[a]` with W x_a + w = x_b + offsets[a], up to a Cartesian distance `tol`.
# Returns `nothing` if `op` is not a symmetry of the lattice.
function _atom_map(lattice::Lattice, op::SymmetryOperation, tol::Float64)
    P = natoms(lattice)
    perm = zeros(Int, P)
    offsets = Vector{Vector{Int}}(undef, P)
    for a in 1:P
        y = op(lattice.positions[a, :])
        for b in 1:P
            lattice.types[b] == lattice.types[a] || continue
            d = y - lattice.positions[b, :]
            n = round.(Int, d)
            if norm(lattice.A' * (d - n)) < tol
                perm[a] = b
                offsets[a] = n
                break
            end
        end
        perm[a] == 0 && return nothing
    end
    return isperm(perm) ? (perm, offsets) : nothing
end

# Replace the translation part found by spglib (accurate to about `symprec`) by
# the least-squares translation mapping the atoms onto each other, reduced to [0, 1).
function _refine(lattice::Lattice, op::SymmetryOperation, tol::Float64) :: SymmetryOperation
    amap = _atom_map(lattice, op, tol)
    isnothing(amap) && error("spglib returned an operation that does not map the lattice onto itself: $op. Try a different `symprec`.")
    perm, offsets = amap
    P = natoms(lattice)
    w = sum(lattice.positions[perm[a], :] + offsets[a] - op.W * lattice.positions[a, :] for a in 1:P) / P
    return SymmetryOperation(op.W, _reduce_mod1(w))
end

# identity first, then pure (fractional) translations, then the remaining operations
_sortkey(op::SymmetryOperation) = (op.W == I ? 0 : 1, vec(op.W)..., round.(op.w; digits=8)...)


# ----------------------------------------------------------------------
#                    Space group of an infinite lattice
# ----------------------------------------------------------------------

@doc raw"""
    SpaceGroup

Symmetry group of an infinite [`Lattice`](@ref) as determined by
[spglib](https://spglib.readthedocs.io), see [`spacegroup`](@ref).

All atom positions and their `types` are taken into account. For a two-dimensional
lattice the plane group (wallpaper group) is determined.

# Fields
- `lattice::Lattice`: the lattice.
- `operations::Vector{SymmetryOperation}`: one operation per coset of the Bravais
  translations of `lattice`, i.e. translation parts reduced to ``[0, 1)``. The identity
  comes first. If the unit cell of `lattice` is not primitive, pure translations by
  fractions of the lattice vectors (rotational part = identity) are included.
- `number::Int`: ITA number of the plane group (1–17) in 2D or of the space group
  (1–230) in 3D.
- `symbol::String`: Hermann–Mauguin symbol, e.g. `"p6mm"` or `"Fd-3m"`.
- `pointgroup::String`: Hermann–Mauguin symbol of the point group, e.g. `"6mm"` or `"m-3m"`.
- `schoenflies::String`: Schoenflies symbol of the point group, e.g. `"C6v"` or `"Oh"`.
- `symmorphic::Bool`: whether the group is symmorphic.
- `origin::Union{Nothing, Vector{Float64}}`: for symmorphic groups, a point (in the lattice
  basis) whose site-symmetry group is the full point group, namely the standard origin
  chosen by spglib. `nothing` for non-symmorphic groups.
- `ntranslations::Int`: number of pure translations modulo the Bravais lattice, i.e. 1 if the
  unit cell is primitive.
- `symprec::Float64`: tolerance (Cartesian distance) used to find the symmetries.
"""
struct SpaceGroup
    lattice::Lattice
    operations::Vector{SymmetryOperation}
    number::Int
    symbol::String
    pointgroup::String
    schoenflies::String
    symmorphic::Bool
    origin::Union{Nothing, Vector{Float64}}
    ntranslations::Int
    symprec::Float64
end

# Convert a Latlib lattice into an spglib cell (two-dimensional lattices become a
# planar layer stacked with a vacuum along z). Returns the cell and the 3x3 matrix
# whose columns are the lattice vectors of the cell.
function _spglib_cell(lattice::Lattice)
    A = lattice.A
    P = natoms(lattice)
    if dim(lattice) == 3
        latmat = Matrix(A')
        positions = [lattice.positions[i, :] for i in 1:P]
    elseif dim(lattice) == 2
        # the vacuum must be much longer than (and incommensurate with) the in-plane
        # vectors, otherwise spglib could find spurious symmetries mixing in- and out-of-plane directions
        vacuum = 10 * sqrt(2) * maximum(norm(A[i, :]) for i in 1:2)
        latmat = [A' zeros(2); 0.0 0.0 vacuum]
        positions = [vcat(lattice.positions[i, :], 0.0) for i in 1:P]
    else
        error(@sprintf "Symmetries are only implemented for 2D and 3D lattices (got dimension %d)." dim(lattice))
    end
    return Spglib.SpglibCell(latmat, positions, lattice.types), latmat
end

# whether the pure rotations about `origin` are part of the group of `ops`
function _is_symmorphic_origin(ops::Vector{SymmetryOperation}, origin::Vector{Float64}) :: Bool
    for W in unique(op.W for op in ops)
        target = SymmetryOperation(W, (I - W) * origin)
        any(op -> _equal_mod_lattice(op, target; atol=1e-6), ops) || return false
    end
    return true
end

@doc raw"""
    spacegroup(lattice::Lattice; symprec=1e-5) -> SpaceGroup
    spacegroup(flattice::FiniteLattice; symprec=1e-5) -> FiniteSpaceGroup

Determine the symmetry group of an infinite lattice (a [`SpaceGroup`](@ref)) or of a
periodic finite lattice (a [`FiniteSpaceGroup`](@ref)) using
[spglib](https://spglib.readthedocs.io).

Atom positions and the atom `types` of the [`Lattice`](@ref) are taken into account;
couplings between the atoms are not. Two-dimensional lattices are treated as a planar
layer, so that the plane group is found. Both symmorphic and non-symmorphic groups are
supported.

For a [`FiniteLattice`](@ref), the symmetry group consists of all operations of the
infinite lattice whose rotational part maps the torus spanned by the boundary vectors
onto itself, combined with all Bravais translations of the cluster. Only fully periodic
finite lattices are supported.

# Keyword arguments
- `symprec::Float64=1e-5`: tolerance (Cartesian distance) up to which atoms are
  considered to be mapped onto each other.

# Examples
```julia
sg = spacegroup(honeycomb)
sg.symbol        # "p6mm"
sg.origin        # [2/3, 2/3]: center of a hexagon

fsg = spacegroup(FiniteLattice(honeycomb, [2 -2; 1 1], true))
fsg.symbol       # "c2mm": the rectangular torus breaks the sixfold rotation
length(fsg)      # 16 operations = 4 point-group operations x 4 translations
```
"""
function spacegroup(lattice::Lattice; symprec::Float64=1e-5) :: SpaceGroup
    D = dim(lattice)
    cell, _ = _spglib_cell(lattice)
    dataset = Spglib.get_dataset(cell, symprec)
    tol = 10 * symprec

    # collect the operations (in 2D: restricted to the plane) without duplicates
    ops = SymmetryOperation[]
    for (W3, w3) in zip(dataset.rotations, dataset.translations)
        if D == 2 && !(all(W3[1:2, 3] .== 0) && all(W3[3, 1:2] .== 0) && is_whole(w3[3]; atol=tol))
            error("spglib returned an operation that does not map the plane onto itself. This is a bug, please report!")
        end
        op = _refine(lattice, SymmetryOperation(Matrix{Int}(W3[1:D, 1:D]), Vector{Float64}(w3[1:D])), tol)
        any(o -> _equal_mod_lattice(o, op), ops) || push!(ops, op)
    end
    sort!(ops; by=_sortkey)

    # type of the group
    if D == 2
        haskey(LAYER_TO_PLANE_GROUP, dataset.spacegroup_number) || error(@sprintf "Unexpected space group %d of a planar layer. This is a bug, please report!" dataset.spacegroup_number)
        number, symbol, pointgroup, schoenflies = LAYER_TO_PLANE_GROUP[dataset.spacegroup_number]
        symmorphic = number in SYMMORPHIC_PLANE_GROUPS
    else
        sgtype = Spglib.get_spacegroup_type(dataset.hall_number)
        number = dataset.spacegroup_number
        symbol = dataset.international_symbol
        pointgroup = sgtype.pointgroup_international
        schoenflies = sgtype.pointgroup_schoenflies
        symmorphic = number in SYMMORPHIC_SPACE_GROUPS
    end

    # standard origin of spglib: x_std = P x + p, so that x = -P^{-1} p is mapped to the origin
    origin = nothing
    if symmorphic
        origin = _reduce_mod1(-(dataset.transformation_matrix \ dataset.origin_shift)[1:D]; atol=1e-8)
        _is_symmorphic_origin(ops, origin) || error("Could not verify the symmorphic origin found by spglib. This is a bug, please report!")
    end

    ntranslations = count(op -> op.W == I, ops)
    return SpaceGroup(lattice, ops, number, symbol, pointgroup, schoenflies, symmorphic, origin, ntranslations, symprec)
end


# ----------------------------------------------------------------------
#                   Space group of a periodic finite lattice
# ----------------------------------------------------------------------

@doc raw"""
    FiniteSpaceGroup

Symmetry group of a periodic [`FiniteLattice`](@ref) (a cluster on a torus), see
[`spacegroup`](@ref).

It consists of all operations of the infinite lattice's space group whose rotational
part maps the torus (the lattice spanned by the boundary vectors) onto itself, combined
with all Bravais translations of the cluster.

# Fields
- `flattice::FiniteLattice`: the finite lattice.
- `spacegroup::SpaceGroup`: symmetry group of the infinite lattice `flattice.lattice`.
- `operations::Vector{SymmetryOperation}`: all operations of the cluster. The rotational
  part runs in the outer loop and the Bravais translations in the inner loop: first the
  zero translation, then the remaining cells in the order of [`bravais_cells`](@ref).
  Hence the identity comes first.
- `permutations::Vector{Vector{Int}}`: site permutations. `permutations[k][i]` is the index
  (in [`atoms`](@ref)`(flattice)`) of the image of site `i` under `operations[k]`.
- `number`, `symbol`, `pointgroup`, `schoenflies`, `symmorphic`: type of the symmetry
  group of the cluster, i.e. of the subgroup of the infinite lattice's space group whose
  rotational parts preserve the torus. Same conventions as in [`SpaceGroup`](@ref).
"""
struct FiniteSpaceGroup
    flattice::FiniteLattice
    spacegroup::SpaceGroup
    operations::Vector{SymmetryOperation}
    permutations::Vector{Vector{Int}}
    number::Int
    symbol::String
    pointgroup::String
    schoenflies::String
    symmorphic::Bool
end

# Identify the space-group type of a set of operations (coset representatives with
# respect to the Bravais lattice of `lattice`) with spglib.
function _identify(lattice::Lattice, ops::Vector{SymmetryOperation}, symprec::Float64)
    _, latmat = _spglib_cell(lattice)
    if dim(lattice) == 3
        rotations = [op.W for op in ops]
        translations = [op.w for op in ops]
    else
        # planar layer: each in-plane operation appears with and without σ_h
        rotations = [[op.W zeros(Int, 2); 0 0 s] for op in ops for s in (1, -1)]
        translations = [vcat(op.w, 0.0) for op in ops for s in (1, -1)]
    end
    sgtype = Spglib.get_spacegroup_type_from_symmetry(rotations, translations, Spglib.Lattice(latmat), symprec)
    if dim(lattice) == 2
        haskey(LAYER_TO_PLANE_GROUP, sgtype.number) || error(@sprintf "Unexpected space group %d of a planar layer. This is a bug, please report!" sgtype.number)
        number, symbol, pointgroup, schoenflies = LAYER_TO_PLANE_GROUP[sgtype.number]
        return number, symbol, pointgroup, schoenflies, number in SYMMORPHIC_PLANE_GROUPS
    else
        return sgtype.number, sgtype.international_short, sgtype.pointgroup_international,
               sgtype.pointgroup_schoenflies, sgtype.number in SYMMORPHIC_SPACE_GROUPS
    end
end

function spacegroup(flattice::FiniteLattice; symprec::Float64=1e-5) :: FiniteSpaceGroup
    if !all(periodicity(flattice))
        throw(ArgumentError("Symmetries can only be determined for fully periodic finite lattices (got periodicity = $(periodicity(flattice)))."))
    end
    lattice = flattice.lattice
    sg = spacegroup(lattice; symprec=symprec)
    P = natoms(lattice)

    # torus: columns of T are the boundary vectors in the lattice basis
    T = Matrix{Int}(flattice.boundary')
    Tinv = inv(Rational{Int}.(T))
    reduce_cell(n) = n - T * floor.(Int, Tinv * n)   # canonical representative modulo the torus (exact)

    # rotational parts mapping the torus onto itself
    reps = [op for op in sg.operations if all(isinteger, Tinv * op.W * T)]

    # Bravais cells of the cluster and index of every site (atom a in cell j)
    cells = [round.(Int, v.coords) for v in bravais_cells(flattice)]
    cell_index = Dict{Vector{Int}, Int}()
    for (j, n) in enumerate(cells)
        reduce_cell(n) == n || error("Unexpected Bravais cell outside of the torus. This is a bug, please report!")
        cell_index[n] = j
    end
    site = zeros(Int, P, length(cells))
    for (s, x) in enumerate(atoms(flattice))
        X = lattice.A' \ x.coords
        for a in 1:P
            d = X - lattice.positions[a, :]
            n = round.(Int, d)
            if norm(lattice.A' * (d - n)) < 1e-6
                site[a, cell_index[reduce_cell(n)]] = s
                break
            end
        end
    end
    any(iszero, site) && error("Could not assign all sites of the finite lattice. This is a bug, please report!")

    # all operations of the cluster and the corresponding site permutations:
    # (W, w + m) maps atom a in cell n onto atom perm[a] in cell offsets[a] + W n + m.
    # The zero translation comes first, so that the identity is the first operation.
    translations = vcat(filter(iszero, cells), filter(!iszero, cells))
    operations = SymmetryOperation[]
    permutations = Vector{Int}[]
    for rep in reps
        perm, offsets = _atom_map(lattice, rep, 10 * symprec)
        for m in translations
            p = zeros(Int, P * length(cells))
            for a in 1:P, (j, n) in enumerate(cells)
                p[site[a, j]] = site[perm[a], cell_index[reduce_cell(offsets[a] + rep.W * n + m)]]
            end
            push!(operations, SymmetryOperation(rep.W, rep.w + m))
            push!(permutations, p)
        end
    end

    number, symbol, pointgroup, schoenflies, symmorphic = _identify(lattice, reps, symprec)
    return FiniteSpaceGroup(flattice, sg, operations, permutations, number, symbol, pointgroup, schoenflies, symmorphic)
end


# ----------------------------------------------------------------------
#                          Accessors and printing
# ----------------------------------------------------------------------

"""
    operations(g::SpaceGroup) -> Vector{SymmetryOperation}
    operations(g::FiniteSpaceGroup) -> Vector{SymmetryOperation}

Symmetry operations of the group, see [`SpaceGroup`](@ref) and [`FiniteSpaceGroup`](@ref).
"""
operations(g::Union{SpaceGroup, FiniteSpaceGroup}) = g.operations

"""
    site_permutations(g::FiniteSpaceGroup) -> Vector{Vector{Int}}

Site permutations of the symmetry operations of a finite lattice: `site_permutations(g)[k][i]`
is the index (in [`atoms`](@ref)) of the image of site `i` under the `k`-th operation.
"""
site_permutations(g::FiniteSpaceGroup) = g.permutations

"""
    issymmorphic(g::SpaceGroup) -> Bool
    issymmorphic(g::FiniteSpaceGroup) -> Bool

Whether the space group (or plane group) is symmorphic.
"""
issymmorphic(g::Union{SpaceGroup, FiniteSpaceGroup}) = g.symmorphic

"""
    pointgroup_operations(g::SpaceGroup) -> Vector{Matrix{Int}}
    pointgroup_operations(g::FiniteSpaceGroup) -> Vector{Matrix{Int}}

Distinct rotational parts (in the lattice basis) of the operations of the group, i.e. the
operations of its point group. The identity comes first.
"""
pointgroup_operations(g::Union{SpaceGroup, FiniteSpaceGroup}) = unique(op.W for op in g.operations)

Base.length(g::Union{SpaceGroup, FiniteSpaceGroup}) = length(g.operations)

function Base.show(io::IO, g::SpaceGroup)
    println(io, "SpaceGroup")
    println(io, @sprintf "group       = %s (#%d)%s" g.symbol g.number (g.symmorphic ? ", symmorphic" : ", non-symmorphic"))
    println(io, @sprintf "point group = %s (%s), %d operations" g.pointgroup g.schoenflies length(pointgroup_operations(g)))
    println(io, @sprintf "operations  = %d (modulo Bravais translations)" length(g))
    g.ntranslations > 1 && println(io, @sprintf "translations: %d pure translations per unit cell (unit cell is not primitive)" g.ntranslations)
    isnothing(g.origin) || println(io, "origin      = ", g.origin)
end

function Base.show(io::IO, g::FiniteSpaceGroup)
    println(io, "FiniteSpaceGroup")
    println(io, @sprintf "group       = %s (#%d)%s" g.symbol g.number (g.symmorphic ? ", symmorphic" : ", non-symmorphic"))
    println(io, @sprintf "point group = %s (%s), %d operations" g.pointgroup g.schoenflies length(pointgroup_operations(g)))
    println(io, @sprintf "operations  = %d on %d sites" length(g) length(first(g.permutations)))
    println(io, @sprintf "lattice     = %s" g.spacegroup.symbol)
end
