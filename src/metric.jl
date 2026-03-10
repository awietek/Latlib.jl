using LinearAlgebra
using Printf

abstract type Metric end

struct EuclideanMetric <: Metric end

struct PeriodicEuclideanMetric <: Metric
    flattice::FiniteLattice 

    # (all necessary consistency checks inside FiniteLattice constructor)
    PeriodicEuclideanMetric(flattice::FiniteLattice) = new(flattice)
end

function (d::EuclideanMetric)(x::EuclideanVector, y::EuclideanVector) :: Float64
    return LinearAlgebra.norm(y.coords - x.coords)
end

# return shortest path from x1 to x2 in euclidean space modulo finite lattice
# TO-DO: fix the problem explained below!
function distance_vector(x1::EuclideanVector, x2::EuclideanVector; flattice=nothing) :: EuclideanVector
    r_euc = x2 - x1
    if isnothing(flattice)
        return r_euc
    else

        # ----- ATTENTION ----- 
        # Naively, one would now expect that round(r_fl.coords) gives the coordinates
        # of the FiniteLatticeVector that is closest to r_fl in euclidean real space distance.
        # For strongly elongated/ non-orthogonal unit cells this is not always the case, however.
        # To properly circumvent this issue, there are well known methods like Niggli/Delaunay 
        # reduction etc.
        # We use a dirty work-around here that hopefully works in most cases but might still fail for extreme cases!!!
        
        r_fl = FiniteLatticeVector(flattice, r_euc)
        n0 = round.(r_fl.coords)  # initial guess via component-wise rounding
        ndim = length(n0)

        # Search all 2^d neighboring integer vectors around n0 to find the true minimum image
        best_dist = Inf
        best_shift = nothing
        for offset in Iterators.product(ntuple(_ -> -1:1, ndim)...)
            n_trial = n0 .+ collect(offset)
            shift = FiniteLatticeVector(flattice, n_trial)
            dist = norm(to_euclidean_basis(shift) - r_euc)
            if dist < best_dist
                best_dist = dist
                best_shift = shift
            end
        end
        return r_euc - to_euclidean_basis(best_shift)
    end
end

function (d::PeriodicEuclideanMetric)(x::EuclideanVector, y::EuclideanVector) :: Float64
    return norm(distance_vector(x, y; flattice=d.flattice))
end

function (d::PeriodicEuclideanMetric)(x::LatticeVector, y::LatticeVector) :: Float64
    return d(to_euclidean_basis(x), to_euclidean_basis(y))
end


@doc raw"""
    distance(x1::AbstractVector, x2::AbstractVector; periodicity_vectors=nothing)

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

`` \min{n_1, \ldots, n_p \in \mathbf{Z}} \lVert \mathbf{x}_1 - \mathbf{x}_2 + \sum_{i=1}^p n_i \mathbf{p}_i\rVert,`` 

where $\mathbf{p}_i$ are the periodicity vectors.
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
    distance_matrix(points::AbstractMatrix; flattice=nothing)

Computes the pairwise distances between points

# Arguments
- `points::AbstractMatrix`: matrix whose columns are the points of which the distance is computed

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

Computes which unique values of distances are present between the points

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

Computes which pairs of the input points are k-th nearest neighbors where k = num_distance.

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
