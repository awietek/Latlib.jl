using Printf
using GLMakie

function plot2d(flattice::FiniteLattice;
                annotate_sites::Bool=true,
                show_boundary::Bool=true,
                show_neighbors::Bool=true)
    f = Figure()
    ax = Axis(f[1,1]) 
    
    coords = coordinates(flattice)
    n_coords = size(coords)[2]

    if show_boundary
        boundary_cartesian = boundary_vectors(flattice)
        A = Tuple([0, 0])
        B = Tuple(boundary_cartesian[:,1])
        C = Tuple(boundary_cartesian[:,1] + boundary_cartesian[:,2])
        D = Tuple(boundary_cartesian[:,2])
        poly!(ax, Point2f[A, B, C, D], color=1, colormap=:tab10, colorrange = (1, 10), alpha=0.2)
    end

    if show_neighbors
        nbors = neighbors(flattice)
        ndist = distances(flattice)[2]

        for n in eachcol(nbors)
            p1 = coords[:,n[1]]
            p2 = coords[:,n[2]]

            # two points differ by periodic direction
            if !isapprox(distance(p1, p2, flattice), distance(p1, p2))
                linestyle = :dash
                d = distance_vector(p2, p1, flattice)
                p1 = p2 + d
                scatter!(ax, [p1[1]], [p1[2]], color=:grey, markersize=12)
                text!(ax, [p1[1]], [p1[2]], text = string.([n[1]]), color=:grey)
            else
                linestyle = :solid
            end
            
            c1 = Tuple(p1)
            c2 = Tuple(p2)
            lines!(ax, [c1, c2], color=2, colormap=:tab10, colorrange = (1, 10),
                   linestyle=linestyle)
        end
    end

    scatter!(ax, coords[1,:], coords[2,:],
             color=1, colormap=:tab10, colorrange = (1, 10),
             markersize=12)
    
    if annotate_sites
        text!(ax, coords[1,:], coords[2,:], text = string.(1:n_coords))
    end

    display(f)
end

"""
    plot

Plots a finite lattice.

# Arguments
- `flattice::FiniteLattice`: finite lattice to plot

# Keyword arguments
- `annotate_sites::Bool=true`: flag to label the lattice sites with a number
- `show_boundary::Bool=true`: flag to show the boundary box of the lattice
- `show_neighbors::Bool=true`: flag to show the nearest neighbors of the lattice
"""
function plot(flattice::FiniteLattice;
              annotate_sites::Bool=true,
              show_boundary::Bool=true,
              show_neighbors::Bool=true)
    dim = dimension(flattice)
    if dim == 2
        plot2d(flattice; annotate_sites, show_boundary, show_neighbors)
    else
        error(@sprintf "Plotting of FiniteLattice not implemented for dimension %d" dim)
    end
end
