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
        # a LatticeVector refers to the basis of its own lattice, an EuclideanVector to Cartesian coordinates
        if !(origin isa LatticeVector || origin isa EuclideanVector)
            throw(ArgumentError("`origin` must be a `LatticeVector` or an `EuclideanVector`."))
        end
        cartesian = origin isa LatticeVector ? to_euclidean_basis(origin).coords : origin.coords
        if length(cartesian) != dim(lattice)
            throw(ArgumentError("`origin` has dimension $(length(cartesian)), the lattice has dimension $(dim(lattice))."))
        end
        x = lattice.A' \ cartesian
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

Two- and three-dimensional lattices with a symmorphic space group are supported so far.

Every operation is written relative to a symmetry center, a point with the full point-group
symmetry of the lattice. Its choice affects the labels of the irreducible representations
away from Gamma. If the lattice has a unique such point (up to lattice vectors), it is used
automatically, e.g. the center of a hexagon of the honeycomb lattice. Otherwise, e.g. for the
square lattice (site or plaquette center), it has to be chosen with the keyword `origin`.

If the unit cell of the lattice is not primitive, only the translations by lattice vectors are
used, and a warning is issued. Operations that act on the sites like another operation (see
[`trivial_operations`](@ref)) are represented by a single site permutation.

# Keyword arguments
- `origin=nothing`: symmetry center, as a [`LatticeVector`](@ref) (coordinates in the basis of its
  lattice) or an [`EuclideanVector`](@ref) (Cartesian coordinates).
- `symprec::Float64=1e-5`: tolerance passed to [`spacegroup`](@ref).

# Examples
```julia
cs = symmetries(FiniteLattice(maple_leaf, [1 1; 1 -2], true))
cs = symmetries(FiniteLattice(square, [4 0; 0 4], true); origin=LatticeVector(square, [0.0, 0.0]))
```
"""
function symmetries(flattice::FiniteLattice; origin=nothing, symprec::Float64=1e-5) :: ClusterSymmetries
    return symmetries(spacegroup(flattice; symprec=symprec); origin=origin)
end

function symmetries(g::FiniteSpaceGroup; origin=nothing) :: ClusterSymmetries
    lattice = g.flattice.lattice
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
    pointgroup_name = _pointgroup_name([lattice.A' * W / lattice.A' for W in pointgroup])
    return ClusterSymmetries(g, c, used, translations, permutations, permutation_index, trivial, pointgroup, pointgroup_name)
end

"""
    momenta(cs::ClusterSymmetries) -> Vector{ClusterMomentum}

All momenta resolved by the finite lattice of `cs`, with their stars under the point group of
the cluster, the representative of each star, labels and little co-groups, see
[`ClusterMomentum`](@ref).
"""
function momenta(cs::ClusterSymmetries)
    lattice = cs.spacegroup.flattice.lattice
    holo = dim(lattice) == 2 ? _holohedry(lattice) : _holohedry3d(lattice)
    return _momenta(cs.spacegroup.flattice, cs.pointgroup, _canonical_frame(cs, holo))
end


# ----------------------------------------------------------------------
#        Conventional frame of the cluster for naming the irreps
# ----------------------------------------------------------------------

# Hermite normal form of the lattice generated by the rows of M (upper triangular, positive
# pivots, entries above a pivot reduced modulo it); unique for a given lattice
function _hermite_normal_form(M::Matrix{Int}) :: Matrix{Int}
    H = copy(M)
    m, n = size(H)
    row = 1
    for col in 1:n
        row > m && break
        while true
            nonzero = [i for i in row:m if H[i, col] != 0]
            isempty(nonzero) && break
            p = nonzero[argmin(abs.(H[nonzero, col]))]
            H[[row, p], :] = H[[p, row], :]
            for i in row+1:m
                H[i, :] -= div(H[i, col], H[row, col]) * H[row, :]
            end
            all(i -> H[i, col] == 0, row+1:m) && break
        end
        H[row, col] == 0 && continue
        H[row, col] < 0 && (H[row, :] = -H[row, :])
        for i in 1:row-1
            H[i, :] -= fld(H[i, col], H[row, col]) * H[row, :]
        end
        row += 1
    end
    return H
end

# The conventional basis (columns, Cartesian) to which the names of the irreps and the numbering of
# repeated momentum labels refer. Equivalent conventional bases of the lattice (images under its
# proper rotations) are distinguished by the cluster: the one is chosen in which the cluster has the
# lexicographically smallest description, given by the Hermite normal form of the torus lattice and
# then by the atom positions relative to the symmetry center (modulo the lattice), both in the
# coordinates of the basis. This does not depend on the Cartesian orientation, the lattice basis or
# the boundary vectors. Two bases with the same description are related by a rotation of the cluster
# and give the same names (the names are covariant under the rotations of the cluster).
function _canonical_frame(cs::ClusterSymmetries, holo) :: Matrix{Float64}
    fl = cs.spacegroup.flattice
    lattice = fl.lattice
    D = dim(lattice)
    A = lattice.A'                         # columns: lattice vectors
    T = A * fl.boundary'                   # columns: torus vectors
    C = Matrix{Float64}(holo.conventional)
    Rs = holo.rotations
    # right-handed: the inversion (3D) or a mirror (2D, if any) of the lattice reverses the handedness
    if det(C) < 0
        flip = findfirst(R -> det(R) < 0, Rs)
        isnothing(flip) || (C = Rs[flip] * C)
    end
    candidates = unique(F -> round.(F; digits=8) .+ 0.0, [R * C for R in Rs if det(R) > 0])
    # 3D: the principal direction of the basis should be invariant (up to sign) under the rotations
    # of the cluster, as it orients the names (automatic except for cubic lattices; impossible only
    # for clusters with a threefold axis, for which every basis is equivalent)
    if D == 3
        P = [A * W / A for W in cs.pointgroup]
        invariant(F) = (p = _directions3d(holo.number, F)[1]; all(Q -> _parallel3(Q * p, p), P))
        any(invariant, candidates) && filter!(invariant, candidates)
    end
    # common denominator of the coordinates of lattice vectors in the conventional basis (centering)
    d = 1
    while !all(is_whole.(d * (C \ A); atol=1e-6))
        d += 1
    end
    # lattice vectors modulo the conventional lattice (centering translations), Cartesian
    centerings = unique(z -> mod.(round.(C \ z; digits=8), 1) .+ 0.0, [A * collect(n) for n in Iterators.product(ntuple(_ -> 0:d, D)...)])
    center = A * cs.origin
    sites = [(lattice.types[j], A * lattice.positions[j, :] - center) for j in 1:size(lattice.positions, 1)]
    # coordinates in the basis F reduced modulo the conventional lattice
    reduced(F, x) = [mod(v, 1) + 0.0 for v in round.(F \ x; digits=8)]
    function key(F)
        hnf = _hermite_normal_form(Matrix(round.(Int, d * (F \ T))'))
        motif = sort([(t, reduced(F, x + z)) for (t, x) in sites for z in centerings])
        return (vec(hnf), motif)
    end
    return candidates[argmin([key(F) for F in candidates])]
end


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

# context to name the irreps of the little co-group of k (holo: holohedry of the lattice, frame:
# conventional basis of the cluster, see `_canonical_frame`)
function _naming_context(lattice::Lattice, holo, frame::AbstractMatrix, k::ClusterMomentum)
    kvec = k.label == "Gamma" ? nothing : _canonical_image(lattice, frame, k.momentum)
    if dim(lattice) == 2
        return _NamingContext(holo.type, frame, holo.shortest, kvec)
    end
    principal, secondary = _directions3d(holo.number, frame)
    return _NamingContext3D(frame, principal, secondary, holo.hexagonal, holo.number in (221, 225, 229), kvec)
end

# All sectors of the representatives of the stars, with a status:
#  :ok         the sector is written,
#  :vanishing  the sector vanishes on the cluster, because it is not trivial on operations
#              acting trivially on the sites,
#  :skipped    a three-dimensional irrep of the little co-group (no characters).
# The sectors are checked against the dimensions of the spin-1/2 Hilbert space, see `_check_sectors`.
function _sectors(cs::ClusterSymmetries)
    lattice = cs.spacegroup.flattice.lattice
    holo = dim(lattice) == 2 ? _holohedry(lattice) : _holohedry3d(lattice)
    ops = operations(cs.spacegroup)
    result = Tuple{Irrep, Symbol}[]
    frame = _canonical_frame(cs, holo)
    ks = _momenta(cs.spacegroup.flattice, cs.pointgroup, frame)
    for k in ks
        k.representative || continue
        gname, Ws, sectors, skipped = _pointgroup_sectors(k.littlegroup, lattice, _naming_context(lattice, holo, frame, k))
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
            push!(result, (irrep, vanishes ? :vanishing : :ok))
        end
        for name in skipped
            irrep = Irrep(k.label * "." * gname * "." * name, k.label, gname, name, name, 3, k.momentum, Int[], ComplexF64[])
            push!(result, (irrep, :skipped))
        end
    end
    sort!(result; by=r -> (_kpoint_sortkey(r[1].kpoint), r[1].name))
    _check_sectors(cs, ks, result)
    return result
end

_skipped_labels(sectors) = [irrep.label for (irrep, status) in sectors if status == :skipped]

function _warn_skipped(skipped::Vector{String})
    isempty(skipped) && return
    @warn "Three-dimensional irreducible representations are skipped, since they cannot be written as one-dimensional characters. The sectors do not span the full Hilbert space. Skipped: " * join(skipped, ", ")
end

@doc raw"""
    irreps(cs::ClusterSymmetries) -> Vector{Irrep}

Irreducible representations (sectors) of the symmetry group of a finite lattice, one set for
the representative of every star of momenta (see [`momenta`](@ref)), with their characters.

For a momentum ``\mathbf{k}`` with little co-group ``P_\mathbf{k}``, the sectors are labelled
`"<momentum>.<P_k>.<irrep>"`, e.g. `"Gamma.C6v.A1"` or `"K.C3v.Ea"`. Two-dimensional irreps of
``P_\mathbf{k}`` are represented by two one-dimensional partners (`a`, `b`), see [`Irrep`](@ref).
Irreps that vanish on the cluster, because they are not trivial on operations acting trivially
on the sites (see [`trivial_operations`](@ref)), are left out. Three-dimensional irreps (`T`
irreps of the cubic point groups) cannot be written as one-dimensional characters; they are
skipped with a warning.

The character of an operation ``\mathcal{X} \mapsto W(\mathcal{X} - \mathbf{c}) + \mathbf{c} + \mathbf{t}``
is ``\rho(W)\, e^{+i \mathbf{k} \cdot \mathbf{t}}``.

# Examples
```julia
cs = symmetries(FiniteLattice(maple_leaf, [1 1; 1 -2], true))
[irrep.label for irrep in irreps(cs)]   # "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", …, "K.C3.Eb"
```
"""
function irreps(cs::ClusterSymmetries)
    sectors = _sectors(cs)
    _warn_skipped(_skipped_labels(sectors))
    return [irrep for (irrep, status) in sectors if status == :ok]
end

function Base.show(io::IO, cs::ClusterSymmetries)
    g = cs.spacegroup
    println(io, "ClusterSymmetries")
    println(io, @sprintf "group       = %s (#%d), point group %s" g.symbol g.number cs.pointgroup_name)
    println(io, "origin      = ", cs.origin)
    println(io, @sprintf "operations  = %d, %d distinct site permutations on %d sites" length(cs.operations) length(cs.permutations) length(cs.permutations[1]))
end
