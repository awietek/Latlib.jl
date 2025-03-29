using LinearAlgebra

abstract type Metric end

struct EuclideanMetric <: Metric end

struct PeriodicMetric{T<:AbstractArray} <: Metric
    periodicity_vectors::T
end

function (dist::PeriodicMetric)(x, y)
    r = dist.periodicity_vectors \ (x - y)
    r = round.(r)
    return norm(dist.periodicity_vectors * r - x + y)
end

function (dist::EuclideanMetric)(x, y)
    return norm(x - y)
end

@doc raw"""
    distance(x1::AbstractVector, x2::AbstractVector; periodicity_vectors=nothing)

Computes the distance between two points

# Arguments
- `x1::AbstractVector`: first point
- `x2::AbstractVector`: second point

# Keyword arguments
- `periodicity_vectors=nothing`: if defined, the vectors that define a periodic direction

If periodicity_vectors is not defined, the euclidian distance is computed

`` d(\mathbf{x}_1, \mathbf{x}_2) = \lVert \mathbf{x}_1 - \mathbf{x}_2 \rVert``.

The periodicity vectors are a matrix ``\mathcal{P}= (\mathbf{p}_1 | \ldots | \mathbf{p}_p)`` 
of dimension ``D \times p`` where ``D`` denotes the dimension of the points, and ``p`` 
denotes the number of periodic directions.

If periodicity_vectors is defined the distance is defined as,

`` d(\mathbf{x}_1, \mathbf{x}_2; \mathcal{P}) \equiv \min\limits_{n_1, \ldots, n_p \in \mathbf{Z}} \lVert \mathbf{x}_1 + \sum_{i=1}^p n_i \mathbf{p}_i - \mathbf{x}_2 \rVert`` 
"""
function distance(x1::AbstractVector, x2::AbstractVector; periodicity_vectors=nothing)
    if periodicity_vectors === nothing
        metric = EuclideanMetric()
    else
        metric = PeriodicMetric(periodicity_vectors)
    end
    return metric(x1, x2)
end


function distance_vector(x1::AbstractVector, x2::AbstractVector; periodicity_vectors=nothing)
    if periodicity_vectors === nothing
        return x2 - x1
    else 
        r = periodicity_vectors \ (x1 - x2)
        r = round.(r)
        return periodicity_vectors * r - x1 + x2
    end
end

"""
    distance_matrix(points::AbstractMatrix; periodicity_vectors=nothing)

Computes the pairwise distances between points

# Arguments
- `points::AbstractMatrix`: matrix whose columns are the points of which the distance is computed

# Keyword arguments
- `periodicity_vectors=nothing`: if defined, the vectors that define a periodic direction
"""
function distance_matrix(points::AbstractMatrix; periodicity_vectors=nothing)
    if periodicity_vectors === nothing
        metric = EuclideanMetric()
    else
        metric = PeriodicMetric(periodicity_vectors)
    end
    dist(x, y) = metric(x, y)

    N = size(points)[2]
    matrix = zeros(Float64, N, N)
    for (i, j) in product(1:N, 1:N)
        matrix[i, j] = dist(points[:,i], points[:,j])
    end
    return matrix
end

"""
    distances(points::AbstractMatrix; periodicity_vectors=nothing)

Computes which unique values of distances are present between the points

# Arguments
- `points::AbstractMatrix`: matrix whose columns are the points of which the distance is computed

# Keyword arguments
- `periodicity_vectors=nothing`: if defined, the vectors that define a periodic direction
"""
function distances(points::AbstractMatrix; periodicity_vectors=nothing)
    return sort(unique(x -> round(x, digits=12), distance_matrix(points; periodicity_vectors)))
end


"""
    neighbors(points::AbstractMatrix; num_distance::Integer=1, periodicity_vectors=nothing)

Computes which pairs of the input points are neighbors

# Arguments
- `points::AbstractMatrix`: matrix whose columns are the points of which the distance is computed

# Keyword arguments
- `num_distance::Integer=1`: at which distance neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.
- `periodicity_vectors=nothing`: if defined, the vectors that define a periodic direction
"""
function neighbors(points::AbstractMatrix; num_distance::Integer=1, periodicity_vectors=nothing)

    N = size(points)[2]
    dists = distances(points; periodicity_vectors=periodicity_vectors)

    if num_distance < 0
        error("Invalid num_distance < 0")
    elseif num_distance > length(dists) - 1
        error("num_distance larger than available distances")
    end
    
    println(dists)
    distance = dists[num_distance+1]

    if periodicity_vectors === nothing
        metric = EuclideanMetric()
    else
        metric = PeriodicMetric(periodicity_vectors)
    end

    dist(x, y) = metric(x, y)
    
    nbors = Vector{Int64}[]
    for (i, j) in product(1:N, 1:N)
        if isapprox(dist(points[:,i], points[:,j]), distance) && i <= j
            push!(nbors, [i, j])
        end
    end

    sort!(nbors)
    return hcat(nbors...)
end
