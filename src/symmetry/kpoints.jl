using LinearAlgebra

# ----------------------------------------------------------------------
#      Holohedry (point group of the Bravais lattice) and labels of
#      momenta in two dimensions (Bilbao / CDML conventions)
# ----------------------------------------------------------------------

# Labels of momenta, in coordinates of the conventional reciprocal basis:
#  points: (label, k)
#  lines:  (label, k0, direction d, αmax) for k = k0 + α d with 0 < α < αmax
# For the centered rectangular lattice, lines are identified from the first-Brillouin-zone
# image instead (Σ: y = 0, Δ: x = 0, C: |y| = 1, F: |x| = 1), since their extent depends on
# the lattice parameters.
const _KPOINT_LABELS_2D = Dict(
    :hexagonal => (
        points = [("K", [1/3, 1/3]), ("M", [1/2, 0.0])],
        lines  = [("Sigma", [0.0, 0.0], [1.0, 0.0], 1/2), ("Lambda", [0.0, 0.0], [1.0, 1.0], 1/3),
                  ("T", [1/2, 0.0], [-1.0, 2.0], 1/6)]),
    :square => (
        points = [("M", [1/2, 1/2]), ("X", [0.0, 1/2])],
        lines  = [("Delta", [0.0, 0.0], [0.0, 1.0], 1/2), ("Sigma", [0.0, 0.0], [1.0, 1.0], 1/2),
                  ("Y", [0.0, 1/2], [1.0, 0.0], 1/2)]),
    :rectangular => (
        points = [("S", [1/2, 1/2]), ("X", [1/2, 0.0]), ("Y", [0.0, 1/2])],
        lines  = [("Sigma", [0.0, 0.0], [1.0, 0.0], 1/2), ("Delta", [0.0, 0.0], [0.0, 1.0], 1/2),
                  ("C", [0.0, 1/2], [1.0, 0.0], 1/2), ("D", [1/2, 0.0], [0.0, 1.0], 1/2)]),
    :centered_rectangular => (
        points = [("Y", [1.0, 0.0]), ("S", [1/2, 1/2])],
        lines  = Tuple{String, Vector{Float64}, Vector{Float64}, Float64}[]),
    :oblique => (
        points = [("Y", [0.0, 1/2]), ("B", [1/2, 0.0]), ("A", [1/2, -1/2])],
        lines  = Tuple{String, Vector{Float64}, Vector{Float64}, Float64}[]),
)

struct _Holohedry
    type::Symbol                          # :oblique, :rectangular, :centered_rectangular, :square, :hexagonal
    rotations::Vector{Matrix{Float64}}    # Cartesian point group of the Bravais lattice
    conventional::Matrix{Float64}         # columns: conventional lattice vectors (Cartesian)
    shortest::Vector{Vector{Float64}}     # shortest lattice vectors (Cartesian)
end

_rotate2d(v::AbstractVector, θ::Real) = [cos(θ) -sin(θ); sin(θ) cos(θ)] * v

# Reduced basis (columns) of the lattice spanned by the columns of B, so that small integer
# combinations reach all short lattice vectors also for a skewed input basis (see `_LatticeReduction`)
_reduced_basis(B::AbstractMatrix) :: Matrix{Float64} = _LatticeReduction(B).basis

# lattice vectors (Cartesian) with coordinates in -n:n in a reduced basis, sorted by length and angle
function _lattice_vectors_2d(lattice::Lattice, n::Int)
    R = _reduced_basis(lattice.A')
    vs = [R * [i, j] for i in -n:n for j in -n:n if (i, j) != (0, 0)]
    return sort(vs; by=v -> (round(norm(v); digits=8), mod(atan(v[2], v[1]), 2π)))
end

function _holohedry(lattice::Lattice) :: _Holohedry
    bravais = Lattice(lattice.A)
    Rs = [cartesian_rotation(op, bravais) for op in operations(spacegroup(bravais))]
    vs = _lattice_vectors_2d(lattice, 6)
    shortest = [v for v in vs if norm(v) < norm(vs[1]) + 1e-8]
    a = shortest[1]
    n = length(Rs)
    if n == 12
        return _Holohedry(:hexagonal, Rs, hcat(a, _rotate2d(a, 2π / 3)), shortest)
    elseif n == 8
        return _Holohedry(:square, Rs, hcat(a, _rotate2d(a, π / 2)), shortest)
    elseif n == 4
        # shortest lattice vectors along the two mirror lines span the conventional cell
        lines = [_mirror_line(R) for R in Rs if det(R) < 0]
        va, vb = [first(v for v in vs if _parallel(v, u)) for u in lines]
        va, vb = norm(va) <= norm(vb) ? (va, vb) : (vb, va)
        vb = va[1] * vb[2] - va[2] * vb[1] > 0 ? vb : -vb
        index = abs(det(hcat(lattice.A' \ va, lattice.A' \ vb)))
        type = isapprox(index, 1; atol=1e-6) ? :rectangular : :centered_rectangular
        return _Holohedry(type, Rs, hcat(va, vb), shortest)
    elseif n == 2
        # reduced basis: shortest vector a, shortest non-parallel vector b with a·b <= 0
        b = first(v for v in vs if !_parallel(v, a) && dot(a, v) <= 1e-8)
        b = a[1] * b[2] - a[2] * b[1] > 0 ? b : -b
        return _Holohedry(:oblique, Rs, hcat(a, b), shortest)
    else
        error("Unexpected point group of order $n of a two-dimensional Bravais lattice. This is a bug, please report!")
    end
end

# whether a Cartesian vector q is a reciprocal lattice vector
_isreciprocal(lattice::Lattice, q::AbstractVector; atol=1e-8) = all(is_whole.(lattice.A * q / (2π); atol=atol))

# all images of a Cartesian momentum k in the first Brillouin zone (several on the zone boundary)
function _bz_images(lattice::Lattice, k::AbstractVector) :: Vector{Vector{Float64}}
    reciprocal = _LatticeReduction(2π * inv(lattice.A))    # columns: reciprocal lattice vectors
    return [[abs(x) < 1e-12 ? 0.0 : x for x in k - G] for G in _closest_lattice_vectors(reciprocal, k)]
end

# image of a Cartesian momentum k in the first Brillouin zone with the largest components along the
# axes of the conventional basis `frame` (columns), lexicographically; on the zone boundary, this
# fixes one of the images independently of the Cartesian orientation
function _canonical_image(lattice::Lattice, frame::AbstractMatrix, k::AbstractVector) :: Vector{Float64}
    units = [normalize(frame[:, j]) for j in 1:size(frame, 2)]
    return argmax(q -> [round(dot(e, q); digits=8) + 0.0 for e in units], _bz_images(lattice, k))
end

# image of a Cartesian momentum k in the first Brillouin zone; on the zone boundary, the
# image with the largest (kx, ky) is chosen
_first_bz(lattice::Lattice, k::AbstractVector) :: Vector{Float64} =
    sort(_bz_images(lattice, k); by=q -> Tuple(round.(q; digits=9) .+ 0.0), rev=true)[1]

# Bilbao label (without numbering) of a Cartesian momentum k
function _kpoint_label_2d(lattice::Lattice, holo::_Holohedry, k::AbstractVector) :: String
    _isreciprocal(lattice, k) && return "Gamma"
    k = _first_bz(lattice, k)   # the shifts below are small
    table = _KPOINT_LABELS_2D[holo.type]
    Cstar = 2π * inv(holo.conventional)'   # columns: conventional reciprocal basis vectors
    for (label, kc) in table.points
        kp = Cstar * kc
        any(R -> _isreciprocal(lattice, R * k - kp), holo.rotations) && return label
    end
    if holo.type == :centered_rectangular
        x, y = holo.conventional' * _first_bz(lattice, k) / (2π)
        abs(y) < 1e-8 && return "Sigma"
        abs(x) < 1e-8 && return "Delta"
        abs(abs(y) - 1) < 1e-8 && return "C"
        abs(abs(x) - 1) < 1e-8 && return "F"
        return "GP"
    end
    Bstar = _reduced_basis(2π * inv(lattice.A))
    shifts = [Bstar * collect(n) for n in Iterators.product(-2:2, -2:2)]
    for (label, k0, d, αmax) in table.lines
        k0c, dc = Cstar * k0, Cstar * d
        for R in holo.rotations, G in shifts
            q = R * k - k0c - G
            α = dot(q, dc) / dot(dc, dc)
            if norm(q - α * dc) < 1e-8 && 1e-8 < α < αmax - 1e-8
                return label
            end
        end
    end
    return "GP"
end


_kpoint_label(lattice::Lattice, holo::_Holohedry, k::AbstractVector) = _kpoint_label_2d(lattice, holo, k)


# ----------------------------------------------------------------------
#      Holohedry and labels of momenta in three dimensions
#      (CDML labels of the Bilbao Crystallographic Server, see kpoint_tables.jl)
# ----------------------------------------------------------------------

# A family of momenta (special point, line or plane) k0 + Σ_i t_i d_i, prepared for exact
# membership tests modulo reciprocal lattice vectors in the coordinates κ = A k / 2π
struct _KFamily
    label::String
    k0::Vector{Float64}           # Cartesian
    directions::Matrix{Int}       # columns: primitive integer directions in κ coordinates
    order::Int                    # order of the little co-group in the holohedry
end

struct _Holohedry3D
    number::Int                           # space group of the Bravais lattice (one of the 14 holohedries)
    rotations::Vector{Matrix{Float64}}    # Cartesian point group of the Bravais lattice
    conventional::Matrix{Float64}         # columns: conventional lattice vectors (ITA standard setting, Cartesian)
    principal::Vector{Float64}            # c (b for monoclinic lattices)
    secondary::Vector{Vector{Float64}}    # unit vectors along a, b, c (and a + b for hexagonal lattices)
    hexagonal::Bool
    families::Vector{_KFamily}            # points, then lines, then planes
end

# primitive integer vector parallel to the rational vector v
function _integer_direction(v::AbstractVector) :: Vector{Int}
    r = rationalize.(v; tol=1e-8)
    n = [numerator(x) * (lcm(denominator.(r)) ÷ denominator(x)) for x in r]
    return n .÷ gcd(n)
end

function _holohedry3d(lattice::Lattice) :: _Holohedry3D
    bravais = Lattice(lattice.A)
    cell, latmat = _spglib_cell(bravais)
    dataset = Spglib.get_dataset(cell, 1e-5)
    number = dataset.spacegroup_number
    haskey(_KPOINT_LABELS_3D, number) || error("Unexpected space group $number of a Bravais lattice. This is a bug, please report!")
    C = latmat * inv(dataset.transformation_matrix)   # conventional basis in the Cartesian frame of the lattice
    Rs = [cartesian_rotation(op, bravais) for op in operations(spacegroup(bravais))]
    hexagonal = number in (166, 191)
    principal, secondary = _directions3d(number, C)
    Cstar = 2π * inv(C)'
    families = _KFamily[]
    for (label, k0, ds, order) in sort(_KPOINT_LABELS_3D[number]; by=f -> length(f[3]))
        label == "Gamma" && continue
        k0c = Cstar * k0
        directions = isempty(ds) ? zeros(Int, 3, 0) : hcat([_integer_direction(lattice.A * (Cstar * d) / (2π)) for d in ds]...)
        push!(families, _KFamily(label, k0c, directions, order))
    end
    return _Holohedry3D(number, Rs, C, principal, secondary, hexagonal, families)
end

# principal direction (c, or b for monoclinic lattices) and secondary directions (unit vectors along
# a, b, c, and a + b for hexagonal lattices) of a conventional basis C (columns) of holohedry `number`
function _directions3d(number::Integer, C::AbstractMatrix)
    principal = number in (10, 12) ? C[:, 2] : C[:, 3]
    hexagonal = number in (166, 191)
    secondary = normalize.(hexagonal ? [C[:, 1], C[:, 2], C[:, 1] + C[:, 2], C[:, 3]] : [C[:, 1], C[:, 2], C[:, 3]])
    return principal, secondary
end

# whether the momentum q (Cartesian) lies in the family f, modulo reciprocal lattice vectors
function _in_family(lattice::Lattice, f::_KFamily, q::AbstractVector) :: Bool
    x = lattice.A * (q - f.k0) / (2π)
    m = size(f.directions, 2)
    if m == 0
        return all(is_whole.(x; atol=1e-8))
    elseif m == 1
        # x - t e ∈ Z^3 for some t: t is fixed by one component up to |e_j| choices
        e = f.directions[:, 1]
        j = argmax(abs.(e))
        return any(n -> all(is_whole.(x - (x[j] - n) / e[j] * e; atol=1e-8)), 0:abs(e[j]) - 1)
    else
        # x ∈ plane + Z^3 iff w·x ∈ Z for the primitive integer normal w
        w = cross(f.directions[:, 1], f.directions[:, 2])
        w = w .÷ gcd(w)
        return is_whole(dot(w, x); atol=1e-8)
    end
end

function _kpoint_label(lattice::Lattice, holo::_Holohedry3D, k::AbstractVector) :: String
    _isreciprocal(lattice, k) && return "Gamma"
    for f in holo.families
        any(R -> _in_family(lattice, f, R * k), holo.rotations) && return f.label
    end
    return "GP"
end

# order of the little co-group of k in the holohedry (used to check the labels)
_holohedry_littlegroup_order(lattice::Lattice, rotations, k::AbstractVector) = count(R -> _isreciprocal(lattice, R * k - k), rotations)


# ----------------------------------------------------------------------
#                       Momenta of a finite lattice
# ----------------------------------------------------------------------

@doc raw"""
    ClusterMomentum

A momentum resolved by a periodic finite lattice, see [`momenta`](@ref).

# Fields
- `coords::Vector{Rational{Int}}`: coordinates ``\kappa`` in the basis of reciprocal lattice
  vectors, ``\mathbf{k} = \sum_j \kappa_j \mathbf{b}_j`` with ``\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\delta_{ij}``,
  reduced to ``[0, 1)``.
- `momentum::Vector{Float64}`: Cartesian coordinates of the image in the first Brillouin zone.
- `label::String`: label of the momentum following the conventions of the Bilbao
  Crystallographic Server for the Bravais lattice, e.g. `"Gamma"`, `"K"`, `"M"`, `"Sigma"`.
  Generic momenta are labelled `"GP0"`, `"GP1"`, …; other labels are numbered (`"Sigma0"`,
  `"Sigma1"`, …) only if they occur more than once among the representatives of the stars. The
  numbers do not depend on the orientation of the lattice or on the representatives.
- `star::Int`: index of the star (orbit under the point group of the cluster) of the momentum.
- `representative::Bool`: whether the momentum represents its star.
- `littlegroup::Vector{Matrix{Int}}`: little co-group, i.e. the rotations (lattice basis) of the
  cluster that leave the momentum invariant up to a reciprocal lattice vector.
- `littlegroup_name::String`: Schoenflies symbol of the little co-group.
"""
struct ClusterMomentum
    coords::Vector{Rational{Int}}
    momentum::Vector{Float64}
    label::String
    star::Int
    representative::Bool
    littlegroup::Vector{Matrix{Int}}
    littlegroup_name::String
end

# all momenta κ with B κ ∈ Z^D, modulo reciprocal lattice vectors
function _cluster_coords(boundary::Matrix{Int}) :: Vector{Vector{Rational{Int}}}
    Binv = inv(Rational{Int}.(boundary))
    gens = [mod.(Binv[:, j], 1) for j in 1:size(Binv, 2)]
    coords = [zeros(Rational{Int}, size(Binv, 1))]
    seen = Set(coords)
    frontier = copy(coords)
    while !isempty(frontier)
        new = Vector{Rational{Int}}[]
        for κ in frontier, g in gens
            p = mod.(κ + g, 1)
            if !(p in seen)
                push!(seen, p)
                push!(new, p)
                push!(coords, p)
            end
        end
        frontier = new
    end
    length(coords) == abs(round(Int, det(boundary))) || error("Wrong number of momenta. This is a bug, please report!")
    return coords
end

# Momenta of a finite lattice with point group `pointgroup` (lattice basis)
# `frame`: conventional basis (columns, Cartesian) that orders repeated labels, see `_canonical_frame`
function _momenta(flattice::FiniteLattice, pointgroup::Vector{Matrix{Int}}, frame::AbstractMatrix) :: Vector{ClusterMomentum}
    lattice = flattice.lattice
    holo = dim(lattice) == 2 ? _holohedry(lattice) : _holohedry3d(lattice)
    coords = _cluster_coords(flattice.boundary)
    cartesian(κ) = 2π * (lattice.A \ Float64.(κ))
    allimages = [_bz_images(lattice, cartesian(κ)) for κ in coords]
    images = [sort(q; by=v -> Tuple(round.(v; digits=9) .+ 0.0), rev=true)[1] for q in allimages]

    # stars: orbits under κ -> W^{-T} κ
    index = Dict(κ => i for (i, κ) in enumerate(coords))
    Winvt = [Matrix{Int}(round.(Int, inv(W))') for W in pointgroup]
    star = zeros(Int, length(coords))
    nstars = 0
    for i in eachindex(coords)
        star[i] == 0 || continue
        nstars += 1
        for M in Winvt
            star[index[mod.(M * coords[i], 1)]] = nstars
        end
    end

    # representative of each star: first one when sorted by (kx, ky) of the first-BZ image, descending
    order = sortperm([Tuple(round.(q; digits=9) .+ 0.0) for q in images]; rev=true)
    representative = falses(length(coords))
    seen = Set{Int}()
    for i in order
        if !(star[i] in seen)
            push!(seen, star[i])
            representative[i] = true
        end
    end

    # labels; numbered among representatives if repeated, generic momenta always numbered. The
    # numbers follow the stars sorted (descending) by the largest components of their momenta along
    # the axes of `frame`, which does not depend on the orientation of the lattice or on the representative
    base = [_kpoint_label(lattice, holo, cartesian(κ)) for κ in coords]
    units = [normalize(frame[:, j]) for j in 1:size(frame, 2)]
    starkey = Dict{Int, Vector{Float64}}()
    for i in eachindex(coords), q in allimages[i]
        key = [round(dot(e, q); digits=8) + 0.0 for e in units]
        (!haskey(starkey, star[i]) || key > starkey[star[i]]) && (starkey[star[i]] = key)
    end
    reps = sort([i for i in eachindex(coords) if representative[i]]; by=i -> starkey[star[i]], rev=true)
    counts = Dict{String, Int}()
    for i in reps
        counts[base[i]] = get(counts, base[i], 0) + 1
    end
    number = Dict{Int, Int}()   # star -> number
    used = Dict{String, Int}()
    for i in reps
        if base[i] == "GP" || counts[base[i]] > 1
            number[star[i]] = get(used, base[i], 0)
            used[base[i]] = number[star[i]] + 1
        end
    end
    labels = [haskey(number, star[i]) ? base[i] * string(number[star[i]]) : base[i] for i in eachindex(coords)]

    # little co-groups
    result = ClusterMomentum[]
    for (i, κ) in enumerate(coords)
        little = [W for W in pointgroup if all(iszero, mod.(W' * κ - κ, 1))]
        Rs = [lattice.A' * W / lattice.A' for W in little]
        push!(result, ClusterMomentum(κ, images[i], labels[i], star[i], representative[i], little, _pointgroup_name(Rs)))
    end
    return result
end
