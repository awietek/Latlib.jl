using Printf
using GLMakie

function plot2d(flattice::FiniteLattice,
    ax::Makie.Axis;
    annotate_sites::Bool=true,
    show_boundary::Bool=true,
    show_neighbors::Bool=true)

    if dim(flattice) != 2
        error(@sprintf "plot2d only supports 2D lattices, but got dimension %d" dim(flattice))
    end

    coords = atoms(flattice)  # Vector{EuclideanVector}
    n_coords = length(coords)

    if show_boundary
        bvecs = [to_euclidean_basis(v) for v in boundary(flattice)]
        A = (0.0, 0.0)
        B = Tuple(bvecs[1].coords)
        C = Tuple(bvecs[1].coords + bvecs[2].coords)
        D = Tuple(bvecs[2].coords)
        poly!(ax, Point2f[A, B, C, D], color=1, colormap=:tab10, colorrange=(1, 10), alpha=0.2)
    end

    if show_neighbors
        nbors = neighbors(flattice)  # Vector{Tuple{Int64, Int64}}
        metric_pbc = PeriodicEuclideanMetric(flattice)
        metric_euc = EuclideanMetric()
        for (i, j) in nbors
            p1 = coords[i]
            p2 = coords[j]

            # two points differ by periodic direction
            if !isapprox(metric_pbc(p1, p2), metric_euc(p1, p2))
                linestyle = :dash
                d = distance_vector(p2, p1; flattice=flattice)
                p1_shifted = p2 + d
                scatter!(ax, [p1_shifted.coords[1]], [p1_shifted.coords[2]], color=:grey, markersize=12)
                if annotate_sites
                    text!(ax, [p1_shifted.coords[1]], [p1_shifted.coords[2]], text=string.([i]), color=:grey)
                end
                c1 = Tuple(p1_shifted.coords)
            else
                linestyle = :solid
                c1 = Tuple(p1.coords)
            end

            c2 = Tuple(p2.coords)
            lines!(ax, [c1, c2], color=2, colormap=:tab10, colorrange=(1, 10),
                linestyle=linestyle)
        end
    end

    xs = [c.coords[1] for c in coords]
    ys = [c.coords[2] for c in coords]
    scatter!(ax, xs, ys,
        color=1, colormap=:tab10, colorrange=(1, 10),
        markersize=12)

    if annotate_sites
        text!(ax, xs, ys, text=string.(1:n_coords))
    end

end


function plot_ops(opsum::OpSum, flattice::FiniteLattice, ax::Makie.Axis;)
    # Assign unique numbers starting from 2 to each bond type
    bond_types = unique(op.cpl for op in opsum.ops)
    color_map = Dict(bond_type => i + 2 for (i, bond_type) in enumerate(bond_types))

    # Get coordinates of lattice points
    coords = atoms(flattice)  # Vector{EuclideanVector}

    metric_pbc = PeriodicEuclideanMetric(flattice)
    metric_euc = EuclideanMetric()

    plotted_bonds = []
    for op in opsum.ops
        # Get the two sites connected by the bond
        site1, site2 = op.sites
        p1 = coords[site1]
        p2 = coords[site2]
        if !isapprox(metric_pbc(p1, p2), metric_euc(p1, p2))
            linestyle = :dash
            d = distance_vector(p2, p1; flattice=flattice)
            p1 = p2 + d
        else
            linestyle = :solid
        end
        # Get the color for the bond type
        bond_color = color_map[op.cpl]

        c1 = Tuple(p1.coords)
        c2 = Tuple(p2.coords)

        label = (op.cpl ∉ plotted_bonds) ? string(op.cpl) : ""
        push!(plotted_bonds, op.cpl)

        if label == string(op.cpl)
            lines!(ax, [c1, c2], color = bond_color, colormap=:tab10, colorrange=(1, 10), linewidth = 2, linestyle=linestyle, label=label)
        else
            lines!(ax, [c1, c2], color = bond_color, colormap=:tab10, colorrange=(1, 10), linewidth = 2, linestyle=linestyle)
        end
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
    d = dim(flattice)

    if d == 2
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
        error(@sprintf "Plotting of FiniteLattice not implemented for dimension %d" d)
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

    d = dim(flattice)

    if d == 2
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
        error(@sprintf "Plotting of FiniteLattice not implemented for dimension %d" d)
    end
end