using LinearAlgebra
using Printf

# ----------------------------------------------------------------------
#      Symmetries of a cluster prepared for irreducible representations
# ----------------------------------------------------------------------

@doc raw"""
    ClusterSymmetries

Symmetry operations of a periodic [`FiniteLattice`](@ref) prepared for the irreducible
representations of its space group, see [`symmetries`](@ref).

Every operation is written relative to the symmetry center ``\mathbf{c}`` (`origin`) as
``\mathcal{X} \mapsto W(\mathcal{X} - \mathbf{c}) + \mathbf{c} + \mathbf{t}`` with a rotation
``W`` and a Bravais translation ``\mathbf{t}`` (both in the lattice basis).

# Fields
- `spacegroup::FiniteSpaceGroup`: symmetry group of the cluster.
- `origin::Vector{Float64}`: symmetry center (lattice basis), a point with the full point-group
  symmetry of the lattice.
- `operations::Vector{Int}`: indices (into [`operations`](@ref)`(spacegroup)`) of the operations
  used. If the unit cell is not primitive, operations with fractional translations are left out.
- `translations::Vector{Vector{Int}}`: translation ``\mathbf{t}`` of each used operation.
- `permutations::Vector{Vector{Int}}`: distinct site permutations of the used operations,
  identity first. These are written to the `Symmetries` section of TOML files.
- `permutation_index::Vector{Int}`: for each used operation, the index of its site permutation
  in `permutations`.
- `trivial::Vector{Int}`: indices (into `operations`) of the used operations that act trivially
  on the sites, see [`trivial_operations`](@ref).
- `pointgroup::Vector{Matrix{Int}}`: point group of the cluster (rotations in the lattice basis),
  identity first.
- `pointgroup_name::String`: Schoenflies symbol of the point group of the cluster.
"""
struct ClusterSymmetries
    spacegroup::FiniteSpaceGroup
    origin::Vector{Float64}
    operations::Vector{Int}
    translations::Vector{Vector{Int}}
    permutations::Vector{Vector{Int}}
    permutation_index::Vector{Int}
    trivial::Vector{Int}
    pointgroup::Vector{Matrix{Int}}
    pointgroup_name::String
end

# Points with the full point-group symmetry of a symmorphic space group, one per class.
# Two points are in the same class if they differ by a lattice vector plus a vector that is
# invariant under the whole point group; only then do they give the same irrep labels.
function _origin_candidates(sg::SpaceGroup) :: Vector{Vector{Float64}}
    D = dim(sg.lattice)
    Ws = pointgroup_operations(sg)
    full_symmetry(Δ) = all(W -> all(is_whole.((I - W) * Δ; atol=1e-9)), Ws)
    shifts = [collect(n) for n in Iterators.product(ntuple(_ -> -1:1, D)...)]
    equivalent(Δ1, Δ2) = any(ℓ -> all(W -> norm((I - W) * (Δ1 - Δ2 - ℓ)) < 1e-9, Ws), shifts)
    classes = Vector{Float64}[]
    # points with full symmetry relative to the origin of spglib have denominators dividing 12
    for n in Iterators.product(ntuple(_ -> 0:11, D)...)
        Δ = collect(n) ./ 12
        if full_symmetry(Δ) && !any(Δc -> equivalent(Δ, Δc), classes)
            push!(classes, Δ)
        end
    end
    return [_reduce_mod1(sg.origin + Δ) for Δ in classes]
end

function _describe_point(lattice::Lattice, x::Vector{Float64}) :: String
    frac = "(" * join(_fraction_string.(x), ", ") * ")"
    cart = "(" * join([@sprintf("%.6f", v) for v in lattice.A' * x], ", ") * ")"
    atom = findfirst(a -> all(is_whole.(x - lattice.positions[a, :]; atol=1e-6)), 1:natoms(lattice))
    site = isnothing(atom) ? "no site" : "site of atom $atom"
    return "$frac in the lattice basis, $cart Cartesian, $site"
end

# symmetry center: given by the user, or unique up to equivalence
function _symmetry_center(sg::SpaceGroup, origin) :: Vector{Float64}
    lattice = sg.lattice
    if !isnothing(origin)
        if origin isa LatticeVector
            x = origin.coords
        elseif origin isa EuclideanVector
            x = lattice.A' \ origin.coords
        else
            throw(ArgumentError("`origin` must be a `LatticeVector` or an `EuclideanVector`."))
        end
        Ws = pointgroup_operations(sg)
        if !all(W -> all(is_whole.((I - W) * (x - sg.origin); atol=1e-6)), Ws)
            throw(ArgumentError("The origin " * _describe_point(lattice, x) * " does not have the full point-group symmetry $(sg.pointgroup) of the lattice."))
        end
        return x
    end
    candidates = _origin_candidates(sg)
    length(candidates) == 1 && return sg.origin
    msg = "The symmetry center is ambiguous: the lattice ($(sg.symbol)) has $(length(candidates)) inequivalent points with the full point-group symmetry $(sg.pointgroup), which lead to different labels of the irreducible representations away from Gamma:\n"
    for x in candidates
        msg *= "  - " * _describe_point(lattice, x) * "\n"
    end
    msg *= "Choose one with the keyword `origin`, e.g. `origin=LatticeVector(lattice, [" * join(_fraction_string.(candidates[1]), ", ") * "])`."
    throw(ArgumentError(msg))
end

@doc raw"""
    symmetries(flattice::FiniteLattice; origin=nothing, symprec=1e-5) -> ClusterSymmetries
    symmetries(g::FiniteSpaceGroup; origin=nothing) -> ClusterSymmetries

Determine the symmetry operations of a periodic finite lattice (or take them from its
[`FiniteSpaceGroup`](@ref)) and prepare them for the irreducible representations of its space
group, see [`irreps`](@ref) and [`momenta`](@ref).

Only two-dimensional lattices with a symmorphic plane group are supported so far.

Every operation is written relative to a symmetry center, a point with the full point-group
symmetry of the lattice. Its choice affects the labels of the irreducible representations
away from Gamma. If the lattice has a unique such point (up to lattice vectors), it is used
automatically, e.g. the center of a hexagon of the honeycomb lattice. Otherwise, e.g. for the
square lattice (site or plaquette center), it has to be chosen with the keyword `origin`.

If the unit cell of the lattice is not primitive, only the translations by lattice vectors are
used, and a warning is issued. Operations that act on the sites like another operation (see
[`trivial_operations`](@ref)) are represented by a single site permutation.

# Keyword arguments
- `origin=nothing`: symmetry center as a [`LatticeVector`](@ref) or [`EuclideanVector`](@ref).
- `symprec::Float64=1e-5`: tolerance passed to [`spacegroup`](@ref).

# Examples
```julia
cs = symmetries(FiniteLattice(maple_leaf, [1 1; 1 -2], true))
cs = symmetries(FiniteLattice(square, [4 0; 0 4], true); origin=LatticeVector(square, [0.0, 0.0]))
```
"""
function symmetries(flattice::FiniteLattice; origin=nothing, symprec::Float64=1e-5) :: ClusterSymmetries
    if dim(flattice.lattice) != 2
        throw(ArgumentError("Irreducible representations are only implemented for two-dimensional lattices so far."))
    end
    return symmetries(spacegroup(flattice; symprec=symprec); origin=origin)
end

function symmetries(g::FiniteSpaceGroup; origin=nothing) :: ClusterSymmetries
    lattice = g.flattice.lattice
    if dim(lattice) != 2
        throw(ArgumentError("Irreducible representations are only implemented for two-dimensional lattices so far."))
    end
    sg = g.spacegroup
    if !sg.symmorphic
        throw(ArgumentError("The plane group $(sg.symbol) of the lattice is non-symmorphic. Irreducible representations are only implemented for symmorphic groups so far."))
    end
    c = _symmetry_center(sg, origin)
    if sg.ntranslations > 1
        @warn "The unit cell of the lattice is not primitive: spglib finds $(sg.ntranslations) pure translations per unit cell (atom types are taken into account, couplings are not). Only translations by lattice vectors are used, so momenta refer to the Brillouin zone of the given unit cell."
    end

    used = Int[]
    translations = Vector{Int}[]
    for (k, op) in enumerate(operations(g))
        t = op.w - (I - op.W) * c
        if all(is_whole.(t; atol=1e-6))
            push!(used, k)
            push!(translations, round.(Int, t))
        end
    end

    permutations = Vector{Int}[]
    index = Dict{Vector{Int}, Int}()
    permutation_index = Int[]
    for k in used
        p = g.permutations[k]
        if !haskey(index, p)
            push!(permutations, p)
            index[p] = length(permutations)
        end
        push!(permutation_index, index[p])
    end
    trivial = [j for j in eachindex(used) if permutation_index[j] == 1]

    pointgroup = unique(operations(g)[k].W for k in used)
    pointgroup_name = _pointgroup_name_2d([lattice.A' * W / lattice.A' for W in pointgroup])
    return ClusterSymmetries(g, c, used, translations, permutations, permutation_index, trivial, pointgroup, pointgroup_name)
end

"""
    momenta(cs::ClusterSymmetries) -> Vector{ClusterMomentum}

All momenta resolved by the finite lattice of `cs`, with their stars under the point group of
the cluster, the representative of each star, labels and little co-groups, see
[`ClusterMomentum`](@ref).
"""
momenta(cs::ClusterSymmetries) = _momenta(cs.spacegroup.flattice, cs.pointgroup)


# ----------------------------------------------------------------------
#                      Irreducible representations
# ----------------------------------------------------------------------

@doc raw"""
    Irrep

A sector of the symmetry group of a finite lattice: a one-dimensional representation of the
little group of a momentum ``\mathbf{k}``, see [`irreps`](@ref).

The character of an operation ``(W, \mathbf{t})`` (relative to the symmetry center) is
``\chi = \rho(W)\, e^{i \mathbf{k} \cdot \mathbf{t}}``, where ``\rho`` is a one-dimensional
representation of the little co-group of ``\mathbf{k}``.

A two-dimensional irrep `parent` of the little co-group (e.g. `E1` of `C6v`) is represented
by its two partners (`E1a`, `E1b`), one-dimensional representations of the subgroup
``\ker(\det E)`` (e.g. `C6`); the two partners are exactly degenerate.

# Fields
- `label::String`: full label `"<momentum>.<little co-group>.<irrep>"`, e.g. `"K.C3v.A1"`.
- `kpoint::String`: label of the momentum, see [`ClusterMomentum`](@ref).
- `littlegroup::String`: Schoenflies symbol of the little co-group.
- `name::String`: name of the irrep (Mulliken symbol; ASCII, with `p`/`pp` for primes).
- `parent::String`: name of the irrep of the little co-group the sector belongs to.
- `dimension::Int`: dimension of the irrep `parent` (1 or 2).
- `momentum::Vector{Float64}`: Cartesian momentum in the first Brillouin zone.
- `allowed_symmetries::Vector{Int}`: indices into the site permutations of
  [`ClusterSymmetries`](@ref) (the `Symmetries` of the TOML file) forming the group of the sector.
- `characters::Vector{ComplexF64}`: characters of the allowed symmetries.
"""
struct Irrep
    label::String
    kpoint::String
    littlegroup::String
    name::String
    parent::String
    dimension::Int
    momentum::Vector{Float64}
    allowed_symmetries::Vector{Int}
    characters::Vector{ComplexF64}
end

function Base.show(io::IO, irrep::Irrep)
    print(io, "Irrep(", irrep.label, ", momentum = ", round.(irrep.momentum; digits=6), ", ",
          length(irrep.allowed_symmetries), " symmetries)")
end

_clean(x::Real) = abs(x) < 1e-14 ? 0.0 : Float64(x)
_clean(z::Complex) = complex(_clean(real(z)), _clean(imag(z)))

# Gamma first, then labelled momenta, then generic ones; numbers sorted numerically
function _kpoint_sortkey(label::String)
    m = match(r"^(.*?)(\d*)$", label)
    base, number = m[1], isempty(m[2]) ? -1 : parse(Int, m[2])
    return (label == "Gamma" ? 0 : base == "GP" ? 2 : 1, base, number)
end

# all sectors of the representatives of the stars, with a flag whether the sector vanishes
# on the cluster (because it is not trivial on operations acting trivially on the sites)
function _sectors(cs::ClusterSymmetries)
    lattice = cs.spacegroup.flattice.lattice
    holo = _holohedry(lattice)
    ops = operations(cs.spacegroup)
    result = Tuple{Irrep, Bool}[]
    for k in momenta(cs)
        k.representative || continue
        ctx = _NamingContext(holo.type, holo.conventional, holo.shortest, k.label == "Gamma" ? nothing : k.momentum)
        gname, Ws, sectors = _pointgroup_sectors(k.littlegroup, lattice, ctx)
        for s in sectors
            exponent = Dict(Ws[s.elements[j]] => s.exponents[j] for j in eachindex(s.elements))
            characters = Dict{Int, ComplexF64}()
            vanishes = false
            for (j, opindex) in enumerate(cs.operations)
                W = ops[opindex].W
                haskey(exponent, W) || continue
                # χ = ρ(W) exp(i k·t) = exp(2πi (a/12 + κ·t)), computed exactly modulo 1
                phase = mod(exponent[W] // 12 + sum(k.coords .* cs.translations[j]), 1)
                χ = cispi(2 * Float64(phase))
                p = cs.permutation_index[j]
                if haskey(characters, p)
                    isapprox(characters[p], χ; atol=1e-10) || (vanishes = true)
                else
                    characters[p] = χ
                end
            end
            allowed = sort(collect(keys(characters)))
            irrep = Irrep(k.label * "." * gname * "." * s.name, k.label, gname, s.name, s.parent, s.dimension,
                          k.momentum, allowed, [_clean(characters[p]) for p in allowed])
            push!(result, (irrep, vanishes))
        end
    end
    return sort(result; by=r -> (_kpoint_sortkey(r[1].kpoint), r[1].name))
end

@doc raw"""
    irreps(cs::ClusterSymmetries) -> Vector{Irrep}

Irreducible representations (sectors) of the symmetry group of a finite lattice, one set for
the representative of every star of momenta (see [`momenta`](@ref)), with their characters.

For a momentum ``\mathbf{k}`` with little co-group ``P_\mathbf{k}``, the sectors are labelled
`"<momentum>.<P_k>.<irrep>"`, e.g. `"Gamma.C6v.A1"` or `"K.C3v.Ea"`. Two-dimensional irreps of
``P_\mathbf{k}`` are represented by two one-dimensional partners (`a`, `b`), see [`Irrep`](@ref).
Irreps that vanish on the cluster, because they are not trivial on operations acting trivially
on the sites (see [`trivial_operations`](@ref)), are left out.

The character of an operation ``\mathcal{X} \mapsto W(\mathcal{X} - \mathbf{c}) + \mathbf{c} + \mathbf{t}``
is ``\rho(W)\, e^{+i \mathbf{k} \cdot \mathbf{t}}``.

# Examples
```julia
cs = symmetries(FiniteLattice(maple_leaf, [1 1; 1 -2], true))
[irrep.label for irrep in irreps(cs)]   # "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", …, "K.C3.Eb"
```
"""
irreps(cs::ClusterSymmetries) = [irrep for (irrep, vanishes) in _sectors(cs) if !vanishes]

function Base.show(io::IO, cs::ClusterSymmetries)
    g = cs.spacegroup
    println(io, "ClusterSymmetries")
    println(io, @sprintf "group       = %s (#%d), point group %s" g.symbol g.number cs.pointgroup_name)
    println(io, "origin      = ", cs.origin)
    println(io, @sprintf "operations  = %d, %d distinct site permutations on %d sites" length(cs.operations) length(cs.permutations) length(cs.permutations[1]))
end
