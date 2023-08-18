using Spglib
using IterTools

struct Lattice
    lattice::Matrix{Real}
    basis::Matrix{Real}
    boundary::Matrix{Integer}
    types::Union{Vector{AbstractString}, Nothing}

    function Lattice(lattice, basis, boundary, types)
        dims_lattice = size(lattice)
        dims_basis = size(basis)
        dims_boundary = size(boundary)

        # Basic consistency checks
        if dims_lattice[2] > 3
            error("Invalid lattice dimension. Must be <= 3")
        elseif dims_lattice[1] > 3 || dims_lattice[1] < dims_lattice[2]
            error("Invalid lattice coordinate dimension. Must be <= 3 and >= dim")
        elseif dims_basis[1] != dims_lattice[1]
            error("Incompatible lattice and basis dimensions")
        elseif dims_boundary[1] != dims_lattice[2] || dims_boundary[2] != dims_lattice[2]
            error("Incompatible lattice dimension and boundary")
        elseif types != nothing
            if length(types) != dims_basis[2]
                error("Incompatible length of types and basis")
            end                
        end
        new(lattice, basis, boundary, types)
    end
end


Lattice(lattice, basis, boundary) = Lattice(lattice, basis, boundary, nothing)

function dimension(lattice::Lattice, coord_dimension=false)
    if coord_dimension
        return size(lattice.lattice)[1]
    else
        return size(lattice.lattice)[2]
    end
end

function bravais_coordinates(lattice::Lattice)
    boundary_lattice = lattice.lattice * lattice.boundary
    dim = dimension(lattice)
    dim_coord = dimension(lattice, true)

    # compute Bravais coordinates within the boundary
    bravais_coords = Vector{Float64}[]
    max_mult = sum(abs.(lattice.boundary), dims=1)
    ranges = [-mm:mm for mm in max_mult]
    for c in product(ranges...)
        bravais_coord = convert.(Float64, lattice.lattice) * convert.(Float64, collect(c))
        frac_coord = boundary_lattice \ bravais_coord
        incell = all(frac_coord .> -1e-6) && all(frac_coord .< 1 - 1e-6)
        if incell
            push!(bravais_coords, bravais_coord)
        end
    end
    sort!(bravais_coords)
    return hcat(bravais_coords...)
end

function coordinates(lattice::Lattice)
    bravais_coords = bravais_coordinates(lattice)

    # generate (unsorted) list of all coordinates
    all_coords = Vector{Float64}[]
    for bravais_coord in eachcol(bravais_coords)
        for basis_coord in eachcol(lattice.basis)
            push!(all_coords, bravais_coord + basis_coord)
        end
    end

    # sort the coordinates w.r.t. given ordering
    sort!(all_coords)
    return hcat(all_coords...)
end
