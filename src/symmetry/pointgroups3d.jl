using LinearAlgebra

# ----------------------------------------------------------------------
#      Three-dimensional crystallographic point groups: identification
#      and Mulliken names of their irreps (Bilbao conventions)
# ----------------------------------------------------------------------

# rotation types (det, trace) -> symbol, as in the International Tables
function _rotation_type_3d(R::AbstractMatrix) :: String
    t = round(Int, tr(R))
    if det(R) > 0
        return Dict(3 => "1", -1 => "2", 0 => "3", 1 => "4", 2 => "6")[t]
    else
        return Dict(-3 => "-1", 1 => "m", 0 => "-3", -1 => "-4", -2 => "-6")[t]
    end
end

const _ROTATION_TYPES_3D = ["1", "2", "3", "4", "6", "-1", "m", "-3", "-4", "-6"]

# numbers of elements of each rotation type -> Schoenflies symbol (all 32 point groups)
const _POINTGROUPS_3D = Dict{Vector{Int}, String}(
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0] => "C1",
    [1, 0, 0, 0, 0, 1, 0, 0, 0, 0] => "Ci",
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0] => "C2",
    [1, 0, 0, 0, 0, 0, 1, 0, 0, 0] => "Cs",
    [1, 1, 0, 0, 0, 1, 1, 0, 0, 0] => "C2h",
    [1, 3, 0, 0, 0, 0, 0, 0, 0, 0] => "D2",
    [1, 1, 0, 0, 0, 0, 2, 0, 0, 0] => "C2v",
    [1, 3, 0, 0, 0, 1, 3, 0, 0, 0] => "D2h",
    [1, 1, 0, 2, 0, 0, 0, 0, 0, 0] => "C4",
    [1, 1, 0, 0, 0, 0, 0, 0, 2, 0] => "S4",
    [1, 1, 0, 2, 0, 1, 1, 0, 2, 0] => "C4h",
    [1, 5, 0, 2, 0, 0, 0, 0, 0, 0] => "D4",
    [1, 1, 0, 2, 0, 0, 4, 0, 0, 0] => "C4v",
    [1, 3, 0, 0, 0, 0, 2, 0, 2, 0] => "D2d",
    [1, 5, 0, 2, 0, 1, 5, 0, 2, 0] => "D4h",
    [1, 0, 2, 0, 0, 0, 0, 0, 0, 0] => "C3",
    [1, 0, 2, 0, 0, 1, 0, 2, 0, 0] => "C3i",
    [1, 3, 2, 0, 0, 0, 0, 0, 0, 0] => "D3",
    [1, 0, 2, 0, 0, 0, 3, 0, 0, 0] => "C3v",
    [1, 3, 2, 0, 0, 1, 3, 2, 0, 0] => "D3d",
    [1, 1, 2, 0, 2, 0, 0, 0, 0, 0] => "C6",
    [1, 0, 2, 0, 0, 0, 1, 0, 0, 2] => "C3h",
    [1, 1, 2, 0, 2, 1, 1, 2, 0, 2] => "C6h",
    [1, 7, 2, 0, 2, 0, 0, 0, 0, 0] => "D6",
    [1, 1, 2, 0, 2, 0, 6, 0, 0, 0] => "C6v",
    [1, 3, 2, 0, 0, 0, 4, 0, 0, 2] => "D3h",
    [1, 7, 2, 0, 2, 1, 7, 2, 0, 2] => "D6h",
    [1, 3, 8, 0, 0, 0, 0, 0, 0, 0] => "T",
    [1, 3, 8, 0, 0, 1, 3, 8, 0, 0] => "Th",
    [1, 9, 8, 6, 0, 0, 0, 0, 0, 0] => "O",
    [1, 3, 8, 0, 0, 0, 6, 0, 6, 0] => "Td",
    [1, 9, 8, 6, 0, 1, 9, 8, 6, 0] => "Oh",
)

# names of the three-dimensional irreps (skipped, since they cannot be written as characters)
const _IRREPS_3D = Dict("T" => ["T"], "Th" => ["Tg", "Tu"], "O" => ["T1", "T2"], "Td" => ["T1", "T2"],
                        "Oh" => ["T1g", "T1u", "T2g", "T2u"])

function _pointgroup_name_3d(Rs::Vector{<:AbstractMatrix}) :: String
    types = [_rotation_type_3d(R) for R in Rs]
    counts = [count(==(t), types) for t in _ROTATION_TYPES_3D]
    haskey(_POINTGROUPS_3D, counts) || error("Unknown point group with rotation types $counts. This is a bug, please report!")
    return _POINTGROUPS_3D[counts]
end

_pointgroup_name(Rs::Vector{<:AbstractMatrix}) = size(first(Rs), 1) == 2 ? _pointgroup_name_2d(Rs) : _pointgroup_name_3d(Rs)

# proper part of an operation and the (unit) axis of a non-trivial proper part
_proper(R::AbstractMatrix) = det(R) > 0 ? R : -R
function _rotation_axis(R::AbstractMatrix) :: Vector{Float64}
    P = _proper(R)
    u = nullspace(P - I; atol=1e-8)
    size(u, 2) == 1 || error("No unique rotation axis. This is a bug, please report!")
    return u[:, 1]
end

_parallel3(u::AbstractVector, v::AbstractVector) = norm(cross(u, v)) < 1e-8 * norm(u) * norm(v)

# counter-clockwise angle of the proper part of R about the oriented axis u
function _rotation_angle(R::AbstractMatrix, u::AbstractVector) :: Float64
    P = _proper(R)
    e = abs(u[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    v = normalize(cross(u, e))
    return atan(dot(cross(v, P * v), u), dot(v, P * v))
end

# Context to orient the Mulliken labels (Bilbao conventions)
struct _NamingContext3D
    conventional::Matrix{Float64}            # columns: conventional lattice vectors a, b, c (Cartesian)
    principal::Vector{Float64}               # principal direction of the lattice: c (b for monoclinic lattices)
    secondary::Vector{Vector{Float64}}       # secondary directions: a, b, c (and a + b for hexagonal lattices)
    hexagonal::Bool                          # hexagonal or trigonal lattice
    cubic::Bool                              # cubic lattice: the principal direction is not distinguished
    k::Union{Nothing, Vector{Float64}}       # momentum (Cartesian, first Brillouin zone); nothing at Γ
end

# orientation of an axis: along the principal direction of the lattice, else along c, b, a
function _orient(u::AbstractVector, ctx::_NamingContext3D) :: Vector{Float64}
    for ref in (ctx.principal, ctx.conventional[:, 3], ctx.conventional[:, 2], ctx.conventional[:, 1])
        s = dot(u, ref)
        abs(s) > 1e-8 * norm(ref) && return s > 0 ? Vector{Float64}(u) : -Vector{Float64}(u)
    end
    return Vector{Float64}(u)
end

_along_secondary(v::AbstractVector, ctx::_NamingContext3D) = any(s -> _parallel3(v, s), ctx.secondary)

# Deterministic choice among elements whose axes (or mirror normals) are not distinguished by
# the rules: the one whose axis, in the conventional basis and with its first nonzero component
# positive, is lexicographically largest. Does not depend on the order of the elements.
function _lexicographic(Rs, list::Vector{Int}, ctx::_NamingContext3D) :: Int
    function key(i)
        x = ctx.conventional \ _rotation_axis(Rs[i])
        x = round.(x / maximum(abs, x); digits=8)
        j = findfirst(v -> v != 0, x)
        return x[j] < 0 ? -x : x
    end
    return list[argmax([Tuple(key(i)) for i in list])]
end

# Orientation data of a point group, needed for the Mulliken names
struct _Orientation
    u::Union{Nothing, Vector{Float64}}   # oriented principal axis
    n::Int                               # order of the principal generator (0 if none)
    generator::Int                       # rotation by +2π/n about u (for S4 and D2d: the -4 operation with proper part +90°)
    inversion::Int                       # index of the inversion (0 if none)
    sigma_h::Int                         # mirror perpendicular to u (0 if none)
    c2p::Int                             # twofold axis ⟂ u, secondary class (0 if none)
    sigma_v::Int                         # mirror containing u, secondary class (0 if none)
    axes::Vector{Int}                    # D2, D2h: twofold rotations about z, y, x
end

function _element(Rs, pred)
    i = findfirst(pred, Rs)
    return isnothing(i) ? 0 : i
end

# rotation with proper part of order n about u, by +2π/n (proper or improper as requested)
function _generator(Rs, u, n, proper::Bool)
    return _element(Rs, R -> (det(R) > 0) == proper && !(_proper(R) ≈ I) && _parallel3(_rotation_axis(R), u) &&
                             abs(rem2pi(_rotation_angle(R, u) - 2π / n, RoundNearest)) < 1e-8)
end

# element of `list` whose rotation axis (or mirror normal) is along the first possible secondary
# direction (in the order of `ctx.secondary`), else the lexicographic choice
function _pick_secondary(Rs, list::Vector{Int}, ctx::_NamingContext3D) :: Int
    isempty(list) && return 0
    for sdir in ctx.secondary
        j = findfirst(i -> _parallel3(_rotation_axis(Rs[i]), sdir), list)
        isnothing(j) || return list[j]
    end
    return _lexicographic(Rs, list, ctx)
end

# C2v: the mirror σ(xz) under which B1 is even, from the two mirrors `mvs` containing the
# twofold axis u. The choice only depends on k, u and the lattice, so that it is the same for
# all momenta of a star (related by symmetries of the cluster).
function _c2v_mirror(Rs, mvs::Vector{Int}, u::Vector{Float64}, ctx::_NamingContext3D) :: Int
    normal(i) = _rotation_axis(Rs[i])
    # k ≠ Γ not along the twofold axis: the mirror containing k
    if !isnothing(ctx.k) && !_parallel3(ctx.k, u)
        containing = [i for i in mvs if abs(dot(normal(i), ctx.k)) < 1e-8 * norm(ctx.k)]
        length(containing) == 1 && return containing[1]
    end
    if !ctx.cubic && _parallel3(u, ctx.principal)
        # twofold axis along the principal direction: for hexagonal lattices the mirror whose normal
        # is a secondary direction (as in 6mm), else the mirror containing a (x along a, as m_y in mm2)
        ctx.hexagonal && return _pick_secondary(Rs, mvs, ctx)
    elseif !ctx.cubic
        # twofold axis perpendicular to the principal direction (k along u, or Γ): the mirror
        # perpendicular to the principal direction
        j = findfirst(i -> _parallel3(normal(i), ctx.principal), mvs)
        isnothing(j) || return mvs[j]
    else
        # cubic lattices: the mirror whose normal is a cubic axis, if only one is
        along = [i for i in mvs if _along_secondary(normal(i), ctx)]
        length(along) == 1 && return along[1]
    end
    # the mirror containing the first secondary direction perpendicular to u (x along it)
    for sdir in ctx.secondary
        _parallel3(sdir, u) && continue
        j = findfirst(i -> abs(dot(normal(i), sdir)) < 1e-8 * norm(sdir), mvs)
        isnothing(j) || return mvs[j]
    end
    return _lexicographic(Rs, mvs, ctx)
end

function _orientation(gname::String, Rs::Vector{<:AbstractMatrix}, ctx::_NamingContext3D) :: _Orientation
    inversion = _element(Rs, R -> R ≈ -I)
    of_type(t) = [i for i in eachindex(Rs) if _rotation_type_3d(Rs[i]) == t]
    u, n, generator, axes = nothing, 0, 0, Int[]

    if gname in ("T", "Th", "O", "Td", "Oh")
        # threefold axis for the complex irreps of T and Th: preferably along a + b + c
        threefold = [_rotation_axis(Rs[i]) for i in of_type("3")]
        d = sum(ctx.conventional[:, j] for j in 1:3)
        i = findfirst(v -> _parallel3(v, d), threefold)
        u, n = _orient(threefold[isnothing(i) ? 1 : i], ctx), 3
        generator = _generator(Rs, u, 3, true)
    elseif gname in ("D2", "D2h")
        twofold = of_type("2")
        dirs = [_rotation_axis(Rs[i]) for i in twofold]
        # z: along the principal direction of the lattice, else along c, b, a, else the lexicographic choice
        z = nothing
        for ref in (ctx.principal, ctx.conventional[:, 3], ctx.conventional[:, 2], ctx.conventional[:, 1])
            z = findfirst(v -> _parallel3(v, ref), dirs)
            isnothing(z) || break
        end
        isnothing(z) && (z = findfirst(==(_lexicographic(Rs, twofold, ctx)), twofold))
        rest = [j for j in 1:3 if j != z]
        # x: in the plane of k and z (k ≠ Γ), i.e. along k or perpendicular to the axis y that is
        # perpendicular to k; else along a secondary direction (in the order a, b, c), else the
        # lexicographic choice
        x = nothing
        if !isnothing(ctx.k)
            along_k = [j for j in rest if _parallel3(dirs[j], ctx.k)]
            perpendicular_k = [j for j in rest if abs(dot(dirs[j], ctx.k)) < 1e-8 * norm(ctx.k)]
            if length(along_k) == 1
                x = along_k[1]
            elseif length(perpendicular_k) == 1
                x = only(j for j in rest if j != perpendicular_k[1])
            end
        end
        if isnothing(x)
            for sdir in ctx.secondary
                j = findfirst(j -> _parallel3(dirs[j], sdir), rest)
                isnothing(j) || (x = rest[j]; break)
            end
        end
        isnothing(x) && (x = findfirst(==(_lexicographic(Rs, twofold[rest], ctx)), twofold))
        y = only(j for j in rest if j != x)
        axes = [twofold[z], twofold[y], twofold[x]]
        u = _orient(dirs[z], ctx)
    elseif gname in ("S4", "D2d")
        u, n = _orient(_rotation_axis(Rs[first(of_type("-4"))]), ctx), 4
        generator = _generator(Rs, u, 4, false)
    elseif gname == "Cs"
        u = _orient(_rotation_axis(Rs[first(of_type("m"))]), ctx)
    elseif !(gname in ("C1", "Ci"))
        # principal axis: the proper rotation of highest order
        n = maximum(parse(Int, _rotation_type_3d(R)) for R in Rs if det(R) > 0)
        u = _orient(_rotation_axis(Rs[first(of_type(string(n)))]), ctx)
        generator = _generator(Rs, u, n, true)
    end

    sigma_h = isnothing(u) ? 0 : _element(Rs, R -> _rotation_type_3d(R) == "m" && _parallel3(_rotation_axis(R), u))
    c2p, sigma_v = 0, 0
    if !isnothing(u) && gname != "Cs"
        perpendicular(i) = abs(dot(_rotation_axis(Rs[i]), u)) < 1e-8
        c2s = [i for i in of_type("2") if perpendicular(i)]
        mvs = [i for i in of_type("m") if perpendicular(i)]     # mirrors containing u
        c2p = _pick_secondary(Rs, c2s, ctx)
        if gname in ("D6", "D6h")
            # Bilbao (622, 6/mmm): B1 is even under the twofold axes 2_120, i.e. along the tertiary
            # directions (perpendicular to a secondary direction), unlike 4/mmm (2_100)
            j = findfirst(i -> !_along_secondary(_rotation_axis(Rs[i]), ctx), c2s)
            isnothing(j) || (c2p = c2s[j])
        end
        if gname == "C2v"
            sigma_v = _c2v_mirror(Rs, mvs, u, ctx)
        else
            sigma_v = _pick_secondary(Rs, mvs, ctx)
        end
    end
    return _Orientation(u, n, generator, inversion, sigma_h, c2p, sigma_v, axes)
end

_is_real(a::Vector{Int}) = all(v -> v == 0 || v == 6, a)
_gu(a, o::_Orientation) = o.inversion == 0 ? "" : (a[o.inversion] == 0 ? "g" : "u")
_prime(a, o::_Orientation) = o.sigma_h == 0 || o.inversion != 0 ? "" : (a[o.sigma_h] == 0 ? "p" : "pp")

# name of a complex one-dimensional irrep: E (with index 1/2 for sixfold groups) and a/b
function _complex_name(a::Vector{Int}, o::_Orientation, suffix::String) :: String
    m = mod(a[o.generator] * o.n ÷ 12, o.n)
    index = o.n == 6 ? string(min(m, o.n - m)) : ""
    return "E" * index * suffix * (2m < o.n ? "a" : "b")
end

# Mulliken name of a one-dimensional irrep (exponents `a`) of a 3D point group
function _name_1d_3d(gname::String, a::Vector{Int}, Rs, o::_Orientation) :: String
    suffix = _gu(a, o) * _prime(a, o)
    _is_real(a) || return _complex_name(a, o, suffix)
    gname in ("C1", "Ci", "Cs", "C3", "C3i", "C3h", "T", "Th") && return "A" * suffix
    if gname in ("D2", "D2h")
        positive = [a[i] == 0 for i in o.axes]
        all(positive) && return "A" * suffix
        return "B" * string(findfirst(positive)) * suffix
    end
    if gname in ("O", "Oh")
        c4 = findfirst(R -> _rotation_type_3d(R) == "4", Rs)
        return "A" * (a[c4] == 0 ? "1" : "2") * suffix
    end
    if gname == "Td"
        m = findfirst(R -> _rotation_type_3d(R) == "m", Rs)
        return "A" * (a[m] == 0 ? "1" : "2") * suffix
    end
    letter = a[o.generator] == 0 ? "A" : "B"
    if gname in ("D3", "D4", "D6", "D2d", "D3d", "D3h", "D4h", "D6h")
        return letter * (a[o.c2p] == 0 ? "1" : "2") * suffix
    elseif gname in ("C2v", "C3v", "C4v", "C6v")
        return letter * (a[o.sigma_v] == 0 ? "1" : "2") * suffix
    end
    return letter * suffix   # C2, C2h, C4, C4h, S4, C6, C6h
end

# Mulliken name of a two-dimensional irrep of a 3D point group
function _name_2d_3d(gname::String, E::_Irrep2D, Rs, o::_Orientation) :: String
    name = "E"
    if gname in ("C6v", "D6", "D6h")
        c2 = _element(Rs, R -> det(R) > 0 && _rotation_type_3d(R) == "2" && _parallel3(_rotation_axis(R), o.u))
        name *= real(E.character[c2]) < 0 ? "1" : "2"
    end
    o.inversion != 0 && (name *= real(E.character[o.inversion]) > 0 ? "g" : "u")
    o.inversion == 0 && o.sigma_h != 0 && (name *= real(E.character[o.sigma_h]) > 0 ? "p" : "pp")
    return name
end

# all sectors of a three-dimensional point group, and the names of its skipped 3D irreps
function _pointgroup_sectors(Ws::Vector{Matrix{Int}}, lattice::Lattice, ctx::_NamingContext3D)
    G = _FiniteGroup(Ws)
    Ws = G.elements
    Rs = [lattice.A' * W / lattice.A' for W in Ws]
    gname = _pointgroup_name_3d(Rs)
    irreps1 = _irreps_1d(G)
    irreps2 = _irreps_2d(G, irreps1)
    n1, n2 = length(irreps1), length(irreps2)
    n3, remainder = divrem(length(G) - n1 - 4 * n2, 9)
    if remainder != 0 || n1 + n2 + n3 != _nclasses(G) || n3 != length(get(_IRREPS_3D, gname, String[]))
        error("Incomplete set of irreducible representations for $gname. This is a bug, please report!")
    end
    o = _orientation(gname, Rs, ctx)

    sectors = _PointGroupSector[]
    for a in irreps1
        name = _name_1d_3d(gname, a, Rs, o)
        push!(sectors, _PointGroupSector(name, name, 1, collect(1:length(G)), a))
    end
    for E in irreps2
        parent = _name_2d_3d(gname, E, Rs, o)
        Hidx = [i for i in 1:length(G) if E.det[i] == 0]
        H = _FiniteGroup(Ws[Hidx])
        RH = Rs[Hidx]
        hname = _pointgroup_name_3d(RH)
        oH = _orientation(hname, RH, ctx)
        partners = 0
        for ρ in _irreps_1d(H)
            multiplicity = sum(conj(_root(ρ[k])) * E.character[Hidx[k]] for k in eachindex(Hidx)) / length(Hidx)
            isapprox(multiplicity, 1; atol=1e-8) || continue
            name = _name_1d_3d(hname, ρ, RH, oH)
            startswith(name, parent) && length(name) == length(parent) + 1 ||
                error("Unexpected partner $name of $parent ($gname). This is a bug, please report!")
            push!(sectors, _PointGroupSector(name, parent, 2, Hidx, ρ))
            partners += 1
        end
        partners == 2 || error("Could not split the two-dimensional irrep $parent of $gname. This is a bug, please report!")
    end
    names = [s.name for s in sectors]
    allunique(names) || error("Ambiguous irrep names $names for $gname. This is a bug, please report!")
    return gname, Ws, sort(sectors; by=s -> s.name), get(_IRREPS_3D, gname, String[])
end
