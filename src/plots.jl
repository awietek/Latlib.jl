using Printf
using GLMakie

function plot2d(flattice::FiniteLattice, ax::Makie.Axis;
    annotate_sites::Bool=true,
    show_boundary::Bool=true,
    show_neighbors::Bool=true)


    coords = coordinates(flattice)
    n_coords = size(coords)[2]

    if show_boundary
        boundary_cartesian = boundary_vectors(flattice)
        A = Tuple([0, 0])
        B = Tuple(boundary_cartesian[:, 1])
        C = Tuple(boundary_cartesian[:, 1] + boundary_cartesian[:, 2])
        D = Tuple(boundary_cartesian[:, 2])
        poly!(ax, Point2f[A, B, C, D], color=1, colormap=:tab10, colorrange=(1, 10), alpha=0.2)
    end

    if show_neighbors
        nbors = neighbors(flattice)
        for n in eachcol(nbors)
            p1 = coords[:, n[1]]
            p2 = coords[:, n[2]]

            # two points differ by periodic direction
            if !isapprox(distance(p1, p2, flattice), distance(p1, p2))
                linestyle = :dash
                d = distance_vector(p2, p1, flattice)
                p1 = p2 + d
                scatter!(ax, [p1[1]], [p1[2]], color=:grey, markersize=12)
                if annotate_sites
                    text!(ax, [p1[1]], [p1[2]], text=string.([n[1]]), color=:grey)
                end
            else
                linestyle = :solid
            end

            c1 = Tuple(p1)
            c2 = Tuple(p2)
            lines!(ax, [c1, c2], color=2, colormap=:tab10, colorrange=(1, 10),
                linestyle=linestyle)
        end
    end

    scatter!(ax, coords[1, :], coords[2, :],
        color=1, colormap=:tab10, colorrange=(1, 10),
        markersize=12)

    if annotate_sites
        text!(ax, coords[1, :], coords[2, :], text=string.(1:n_coords))
    end

end


function plot_ops(opsum::OpSum, flattice::FiniteLattice, ax::Makie.Axis;)
    # Assign unique numbers starting from 2 to each bond type
    bond_types = unique(op.coupling for op in opsum.ops)
    color_map = Dict(bond_type => i + 2 for (i, bond_type) in enumerate(bond_types))

    # Get coordinates of lattice points
    coords = coordinates(flattice)

    plotted_bonds = []
    for op in opsum.ops
        # Get the two sites connected by the bond
        site1, site2 = op.sites
        p1 = coords[:, site1]
        p2 = coords[:, site2]
        if !isapprox(distance(p1, p2, flattice), distance(p1, p2))
            linestyle = :dash
            d = distance_vector(p2, p1, flattice)
            p1 = p2 + d
        else
            linestyle = :solid
        end
        # Get the color for the bond type
        bond_color = color_map[op.coupling]

        c1 = Tuple(p1)
        c2 = Tuple(p2)

        label = (op.coupling ∉ plotted_bonds) ? string(op.coupling) : ""
        push!(plotted_bonds, op.coupling)

        if label == op.coupling
            # Plot the line connecting the two points
            lines!(ax, [c1, c2], color = bond_color, colormap=:tab10, colorrange=(1, 10), linewidth = 2, linestyle=linestyle,label=label)
        else
            lines!(ax, [c1, c2], color = bond_color, colormap=:tab10, colorrange=(1, 10), linewidth = 2, linestyle=linestyle)
        end
        #text!(ax, midpoint[1], midpoint[2]*1.1, text = "$(op.coupling)", align = (:center, :center), fontsize = 12,color = bond_color, colormap=:tab10, colorrange=(1, 10))
    end
    axislegend(ax)

end


@doc raw"""
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
    ax=nothing,
    annotate_sites::Bool=true,
    show_boundary::Bool=true,
    show_neighbors::Bool=true)
    dim = dimension(flattice)

    if dim == 2
        if ax == nothing
            f = Figure()
            ax = Axis(f[1, 1])
            show = true
        else
            show = false
        end

        plot2d(flattice, ax; annotate_sites, show_boundary, show_neighbors)
        if show == true
            display(f)
        end
    else
        error(@sprintf "Plotting of FiniteLattice not implemented for dimension %d" dim)
    end
end



@doc raw"""
    plot_opsum

Plots a finite lattice with the interaction bonds contained in the OpSum. The interaction with differnt couplings are represented in distintic colors.

# Arguments
- `opsum::OpSum`: finite lattice to plot
- `flattice::FiniteLattice`: finite lattice to plot

# Keyword arguments
- `annotate_sites::Bool=true`: flag to label the lattice sites with a number
- `show_boundary::Bool=true`: flag to show the boundary box of the lattice
- `show_neighbors::Bool=true`: flag to show the nearest neighbors of the lattice
"""
function plot_opsum(opsum::OpSum, flattice::FiniteLattice;
    ax=nothing,
    annotate_sites::Bool=true,
    show_boundary::Bool=true,
    show_neighbors::Bool=true)

    dim = dimension(flattice)

    if dim == 2
        if ax == nothing
            f = Figure()
            ax = Axis(f[1, 1])
            show = true
        else
            show = false
        end

        plot2d(flattice, ax; annotate_sites, show_boundary, show_neighbors)
        plot_ops(opsum, flattice, ax)

        if show == true
            display(f)
        end
    else
        error(@sprintf "Plotting of FiniteLattice not implemented for dimension %d" dim)
    end
end