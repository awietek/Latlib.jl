using LinearAlgebra

abstract type Metric end

struct EuclideanMetric <: Metric end
struct PeriodicMetric{T<:AbstractArray} <: Metric
    periodicity::T
end

function (dist::PeriodicMetric)(x, y)
    r = dist.periodicity \ (x - y)
    r = round.(r)
    return norm(dist.periodicity * r - x + y)
end

function (dist::EuclideanMetric)(x, y)
    return norm(x - y)
end


function distances(points; periodicity=nothing)
    N = size(points)[2]

    # Choose the correct metric
    if periodicity === nothing
        metric = EuclideanMetric()
    else
        metric = PeriodicMetric(periodicity)
    end
    dist(x, y) = metric(x, y)

    dists = Float64[]
    for (i, j) in product(1:N, 1:N)
        d = dist(points[:,i], points[:,j])
        
        found = false
        for dd in dists
            if isapprox(dd, d)
                found = true
            end
        end
        if !found
            append!(dists, d)
        end
    end
    return sort(dists)
end

function neighbors(coords, num_distance; periodicity=nothing)
    N = size(coords)[2]

    dists = distances(coords; periodicity=periodicity)
    distance = dists[num_distance+1]
    
    nbors_all = Vector{Int64}[]
    for i in 1:N
        nbors = Int64[]
        p = coords[:,i]
        for j in 1:N
            p2 = coords[:,j]
            dist = 0.0
            if periodicity === nothing
                metric = EuclideanMetric()
            else
                metric = PeriodicMetric(periodicity)
            end
            dist(x, y) = metric(x, y)

            
            if isapprox(dist(p, p2), distance)
                append!(nbors, j)
            end
        end
        push!(nbors_all, nbors)
    end
    return nbors_all        
end


function periodicity(lattice::Lattice, periodic_dims::Vector{<:Integer}=[])
    periodicity = nothing
    coord_dim = dimension(lattice, true)
    n_periodic_dims = length(periodic_dims)
    if n_periodic_dims > dimension(lattice)
        error("More periodic dimensions than dimensions of the Lattice")
    end
    if n_periodic_dims > 0
        periodicity = zeros(Float64, coord_dim, n_periodic_dims)
        for dim in 1:n_periodic_dims
            periodicity[:,dim] = lattice.lattice * lattice.boundary[:,periodic_dims[dim]] 
        end
    end
    return periodicity
end
