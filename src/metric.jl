using LinearAlgebra
using Printf

abstract type Metric end

"""
    EuclideanMetric()

The standard Euclidean metric. Calling an instance on two `EuclideanVector`s
returns their Euclidean distance.

```julia
d = EuclideanMetric()
d(EuclideanVector([0, 0]), EuclideanVector([3, 4]))   # 5.0
```
"""
struct EuclideanMetric <: Metric end

"""
    PeriodicEuclideanMetric(flattice::FiniteLattice)

The Euclidean metric modulo the periodic boundary vectors of a [`FiniteLattice`](@ref).
Calling an instance on two `EuclideanVector`s (or `LatticeVector`s) returns the
shortest distance between them, taking into account all periodic images.
Boundary directions with open boundary conditions are treated as non-periodic.

```julia
fl = FiniteLattice(square, [4 0; 0 4], true)
d = PeriodicEuclideanMetric(fl)
d(EuclideanVector([0, 0]), EuclideanVector([3, 0]))   # 1.0 (across the boundary)
```
"""
struct PeriodicEuclideanMetric <: Metric
    flattice::FiniteLattice 

    # (all necessary consistency checks inside FiniteLattice constructor)
    PeriodicEuclideanMetric(flattice::FiniteLattice) = new(flattice)
end

function (d::EuclideanMetric)(x::EuclideanVector, y::EuclideanVector) :: Float64
    return LinearAlgebra.norm(y.coords - x.coords)
end

"""
    distance_vector(x1::EuclideanVector, x2::EuclideanVector; flattice=nothing) -> EuclideanVector

Returns the shortest vector pointing from `x1` to `x2`. If `flattice` is a
[`FiniteLattice`](@ref), the vector is taken modulo the periodic boundary vectors
of the finite lattice, i.e., the returned vector connects `x1` to the closest periodic
image of `x2`. Without `flattice`, simply `x2 - x1` is returned.
"""
function distance_vector(x1::EuclideanVector, x2::EuclideanVector; flattice=nothing) :: EuclideanVector
    r_euc = x2 - x1
    if isnothing(flattice)
        return r_euc
    else

        # get the boundary vectors along which periodicity is assumed, as Vector{LatticeVector}
        periodic_boundary_vecs = periodic_boundary(flattice)
        if length(periodic_boundary_vecs) == 0
            return r_euc # no periodicity -> standard euclidean treatment
        end

        # construct "torus matrix" in euclidean coordinates (boundary vecs in columns)
        T = hcat([to_euclidean_basis(vec).coords for vec in periodic_boundary_vecs]...)

        # subtract the closest vector of the lattice spanned by the periodic boundary vectors (for
        # fewer periodic directions than dimensions, the component of r_euc perpendicular to them
        # does not change which one is closest)
        return r_euc - EuclideanVector(_closest_lattice_vector(_LatticeReduction(T), r_euc.coords))
    end
end


# ----------------------------------------------------------------------
#        Closest lattice vector (exact for lattices of rank <= 3)
# ----------------------------------------------------------------------
#
# The vector of a lattice (spanned by the columns of B, rank n <= 3, possibly embedded in more
# dimensions) closest to a point x is found with the Voronoi-relevant vectors of the lattice. Every
# lattice of rank n <= 3 has an obtuse superbase b_1, ..., b_n, b_0 = -Σ_i b_i with all b_i·b_j <= 0,
# obtained by Lagrange-Gauss (n = 2) or Selling (n = 3) reduction, and the sums over its nonempty
# proper subsets contain all Voronoi-relevant vectors (6 for n = 2, 14 for n = 3; Conway and Sloane,
# Proc. R. Soc. A 436, 55 (1992)). Starting from rounding in the reduced basis, the candidate y is
# moved along a relevant vector as long as this brings it closer to x. This stops exactly when x - y
# lies in the Voronoi cell of the origin, which proves that y is a closest lattice vector. The
# result does not depend on how skewed the given basis is.

struct _LatticeReduction
    basis::Matrix{Float64}               # reduced basis (columns) of the lattice
    pinv::Matrix{Float64}                # its pseudo-inverse
    relevant::Vector{Vector{Float64}}    # contains all Voronoi-relevant vectors
    tol::Float64                         # tolerance for squared distances
end

# obtuse superbase [b_1, ..., b_n, b_0] of the lattice spanned by the columns of B (n <= 3)
function _obtuse_superbase(B::Matrix{Float64}) :: Vector{Vector{Float64}}
    n = size(B, 2)
    bs = [B[:, i] for i in 1:n]
    n == 1 && return [bs[1], -bs[1]]
    # pairwise size reduction with rounded multipliers: Lagrange-Gauss reduction for n = 2, and a
    # fast first step for strongly skewed bases for n = 3
    changed = true
    while changed
        changed = false
        for i in 1:n, j in 1:n
            i == j && continue
            m = round(dot(bs[i], bs[j]) / dot(bs[j], bs[j]))
            if m != 0 && norm(bs[i] - m * bs[j]) < (1 - 1e-12) * norm(bs[i])
                bs[i] -= m * bs[j]
                changed = true
            end
        end
    end
    if n == 2
        # |b_1·b_2| <= |b_i|^2 / 2, so the superbase is obtuse once b_1·b_2 <= 0
        dot(bs[1], bs[2]) > 0 && (bs[2] = -bs[2])
        return [bs[1], bs[2], -bs[1] - bs[2]]
    end
    # Selling reduction: while b_i·b_j > 0, replace b_i -> -b_i and b_k -> b_k + b_i for the other
    # two vectors; this strictly decreases the sum of the squared lengths
    sb = [bs[1], bs[2], bs[3], -bs[1] - bs[2] - bs[3]]
    tol = 1e-12 * maximum(b -> dot(b, b), sb)
    pairs = [(i, j) for i in 1:4 for j in i+1:4]
    for _ in 1:10_000
        p = findfirst(((i, j),) -> dot(sb[i], sb[j]) > tol, pairs)
        isnothing(p) && return sb
        i, j = pairs[p]
        for k in 1:4
            (k == i || k == j) || (sb[k] += sb[i])
        end
        sb[i] = -sb[i]
    end
    error("Selling reduction of the lattice spanned by $B did not converge. This is a bug, please report!")
end

function _LatticeReduction(B::AbstractMatrix)
    B = Matrix{Float64}(B)
    n = size(B, 2)
    1 <= n <= 3 || throw(ArgumentError("Closest lattice vectors are implemented for lattices of rank 1 to 3, got rank $n."))
    superbase = _obtuse_superbase(B)
    R = hcat(superbase[1:n]...)
    if cond(R) > 1e8
        error("The lattice spanned by the columns of $B is nearly degenerate (condition number $(cond(R)) of its reduced basis), so closest lattice vectors cannot be determined reliably.")
    end
    relevant = [sum(superbase[i] for i in 1:n+1 if isodd(mask >> (i - 1))) for mask in 1:2^(n+1)-2]
    return _LatticeReduction(R, pinv(R), relevant, 1e-10 * maximum(v -> dot(v, v), relevant))
end

# a lattice vector closest to x
function _closest_lattice_vector(red::_LatticeReduction, x::AbstractVector) :: Vector{Float64}
    y = red.basis * round.(red.pinv * x)
    r = x - y
    for _ in 1:1000
        gain, i = findmax(v -> 2 * dot(r, v) - dot(v, v), red.relevant)
        gain <= red.tol && return y      # x - y lies in the Voronoi cell of the origin
        y += red.relevant[i]
        r -= red.relevant[i]
    end
    error("No closest lattice vector found for $x. This is a bug, please report!")
end

# all lattice vectors closest to x (several if x lies on the boundary of a Voronoi cell). They are
# the vertices of a face of the Delaunay tiling, connected by its edges, which are relevant vectors.
function _closest_lattice_vectors(red::_LatticeReduction, x::AbstractVector) :: Vector{Vector{Float64}}
    y = _closest_lattice_vector(red, x)
    d2 = sum(abs2, x - y)
    found = [y]
    queue = [y]
    while !isempty(queue)
        z = pop!(queue)
        for v in red.relevant
            w = z + v
            if sum(abs2, x - w) <= d2 + red.tol && !any(u -> sum(abs2, u - w) < red.tol, found)
                push!(found, w)
                push!(queue, w)
            end
        end
    end
    return found
end

function (d::PeriodicEuclideanMetric)(x::EuclideanVector, y::EuclideanVector) :: Float64
    return norm(distance_vector(x, y; flattice=d.flattice))
end

function (d::PeriodicEuclideanMetric)(x::LatticeVector, y::LatticeVector) :: Float64
    return d(to_euclidean_basis(x), to_euclidean_basis(y))
end


@doc raw"""
    distance(x1::EuclideanVector, x2::EuclideanVector; flattice=nothing)

Computes the distance between two points

# Arguments
- `x1::EuclideanVector`: first point
- `x2::EuclideanVector`: second point

# Keyword arguments
- `flattice=nothing`: `FiniteLattice` instance defining the periodicity vectors.

If flattice=nothing, the standard euclidean distance is computed

`` d(\mathbf{x}_1, \mathbf{x}_2) = \lVert \mathbf{x}_1 - \mathbf{x}_2 \rVert_2``.

If flattice is defined (FiniteLattice), then compute the distance as the
minimum distance between x1 and x2 assuming full periodicity along the periodicity vectors.
In other words, this function returns

``\min_{n_1, \ldots, n_p \in \mathbb{Z}} \lVert \mathbf{x}_1 - \mathbf{x}_2 + \sum_{i=1}^p n_i \mathbf{p}_i\rVert,``

where ``\mathbf{p}_i`` are the periodicity vectors.
"""
function distance(x1::EuclideanVector, x2::EuclideanVector; flattice=nothing) :: Float64
    if isnothing(flattice)
        metric = EuclideanMetric()
    else
        metric = PeriodicEuclideanMetric(flattice)
    end
    return metric(x1, x2)
end

@doc """
    distance_matrix(points::Vector{EuclideanVector}; flattice=nothing)

Computes the symmetric matrix of pairwise distances between points.

# Arguments
- `points::Vector{EuclideanVector}`: the points between which distances are computed

# Keyword arguments
- `flattice=nothing`: if defined, the `FiniteLattice` instance that defines the periodicity vectors.
"""
function distance_matrix(xs::Vector{EuclideanVector}; flattice=nothing)
    if isnothing(flattice)
        metric = EuclideanMetric()
    else
        metric = PeriodicEuclideanMetric(flattice)
    end

    N = length(xs)
    matrix = zeros(Float64, N, N)
    for i in 1:N
        for j in (i+1):N
            matrix[i, j] = metric(xs[i], xs[j])
            matrix[j, i] = matrix[i, j] # symmetric
        end
    end
    return matrix
end

@doc """
    distances(points::Vector{EuclideanVector}; flattice=nothing)

Computes the sorted unique values of distances present between the points.
The first entry is always `0.0` (the self-distance).

# Arguments
- `points::Vector{EuclideanVector}`: vector of points of which the distance is computed

# Keyword arguments
- `flattice=nothing`: `FiniteLattice` defining the `PeriodicEuclideanMetric`. If `nothing`, the standard `EuclideanMetric` is used.
"""
function distances(points::Vector{EuclideanVector}; flattice=nothing)
    return sort(unique(x -> round(x, digits=12), distance_matrix(points; flattice=flattice)))
end


@doc """
    neighbors(points::Vector{EuclideanVector}; num_distance::Integer=1, flattice=nothing)
    neighbors(flattice::FiniteLattice; num_distance::Integer=1)

Computes which pairs of the input points are k-th nearest neighbors where k = num_distance.
For a `FiniteLattice`, the points are its sites `atoms(flattice)` and distances are computed
with the periodic metric of the finite lattice.

# Arguments
- `points::Vector{EuclideanVector}`: Vectors taken into account for distance computation.

# Keyword arguments
- `num_distance::Integer=1`: At which distance neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.
- `flattice=nothing`: `FiniteLattice` defining the `PeriodicEuclideanMetric`. If `nothing`, the standard `EuclideanMetric` is used.

# Returns
- neighbors::Vector{Tuple{Int64, Int64}}: Vectors of index pairs (i, j) with i<j of points that are num_distance-th nearest neighbors.
"""
function neighbors(points::Vector{EuclideanVector}; num_distance::Integer=1, flattice=nothing) :: Vector{Tuple{Int64, Int64}}

    N = length(points)
    dists = distances(points; flattice=flattice)

    if num_distance < 1
        error("Invalid num_distance < 1")
    elseif num_distance > length(dists) - 1  # -1 because the first distance is always 0 (self-distance)
        error(@sprintf "Num_distance (%s) is larger than the number of available distances (%s)." num_distance length(dists))
    end
    
    num_distance_val = dists[num_distance+1] # +1 because the first distance is always 0 (self-distance)

    if isnothing(flattice)
        metric = EuclideanMetric()
    else
        metric = PeriodicEuclideanMetric(flattice)
    end

    dist(x, y) = metric(x, y)
    
    result = Vector{Tuple{Int64, Int64}}()
    for i in 1:N
        for j in (i+1):N
            if isapprox(dist(points[i], points[j]), num_distance_val) && i < j
                push!(result, (i, j))
            end
        end
    end

    return result
end

function neighbors(flattice::FiniteLattice; num_distance::Integer=1)
    return neighbors(atoms(flattice); num_distance=num_distance, flattice=flattice)
end
