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



@doc raw"""
    plot_3d

Plots a 3D `FiniteLattice` interactively using GLMakie.

Draws:
- A grid of Bravais lattice vectors (arrows from each unit cell origin)
- The unit cell parallelepiped of the underlying `Lattice`
- All atom sites in the finite lattice, colored by atom type
- Nearest-neighbor bonds

# Arguments
- `flattice::FiniteLattice`: a 3D finite lattice to plot

# Keyword arguments
- `annotate_sites::Bool=false`: label each site with its index
- `show_boundary::Bool=true`: show the boundary parallelepiped of the finite lattice
- `show_unit_cell::Bool=true`: show the unit cell parallelepiped at the origin
- `show_bravais_grid::Bool=true`: show Bravais lattice vectors as arrows at each cell origin
- `show_neighbors::Bool=true`: show nearest-neighbor bonds
"""
function plot_3d(flattice::FiniteLattice;
    show_boundary::Bool=true,
    show_unit_cell::Bool=true,
    show_bravais_grid::Bool=true,
    show_neighbors::Bool=true)

    if dim(flattice) != 3
        error(@sprintf "plot_3d only supports 3D lattices, but got dimension %d" dim(flattice))
    end

    lat = flattice.lattice
    As = lattice_vecs(flattice)
    a1 = As[1].coords # converts EuclideanVector to Vector{Float64}
    a2 = As[2].coords
    a3 = As[3].coords

    f = Figure(size=(1000, 800))
    ax = LScene(f[1, 1]; show_axis=true)

    bravais_vecs = [a1, a2, a3]
    bravais_arrow_color = :black

    # --- Bravais vectors: one arrow per lattice vector, drawn from the origin ---
    if show_bravais_grid
        ox = Float64[]; oy = Float64[]; oz = Float64[]
        dx = Float64[]; dy = Float64[]; dz = Float64[]
        for avec in bravais_vecs
            push!(ox, 0.0); push!(oy, 0.0); push!(oz, 0.0)
            push!(dx, avec[1]); push!(dy, avec[2]); push!(dz, avec[3])
        end
        arrows!(ax, ox, oy, oz, dx, dy, dz;
            color=bravais_arrow_color,
            arrowsize=Vec3f(0.15, 0.15, 0.2),
            linewidth=0.04)
    end

    # --- Unit cell parallelepiped ---
    if show_unit_cell
        _draw_parallelepiped!(ax, [0.0, 0.0, 0.0], a1, a2, a3;
            color=(:dodgerblue, 0.12), linecolor=:dodgerblue, linewidth=2)
    end

    # --- Boundary parallelepiped ---
    if show_boundary
        bvecs = [to_euclidean_basis(v).coords for v in boundary(flattice)]
        _draw_parallelepiped!(ax, [0.0, 0.0, 0.0], bvecs[1], bvecs[2], bvecs[3];
            color=(:orange, 0.06), linecolor=:orange, linewidth=1.5)
    end

    # --- Neighbor bonds ---
    coords = atoms(flattice)
    if show_neighbors
        nbors = neighbors(flattice)
        metric_pbc = PeriodicEuclideanMetric(flattice)
        metric_euc = EuclideanMetric()
        for (i, j) in nbors
            p1 = coords[i]
            p2 = coords[j]
            is_periodic = !isapprox(metric_pbc(p1, p2), metric_euc(p1, p2))
            linestyle = is_periodic ? :dash : :solid
            if is_periodic
                d = distance_vector(p2, p1; flattice=flattice)
                p1_shifted = p2 + d
                c1 = p1_shifted.coords
            else
                c1 = p1.coords
            end
            c2 = p2.coords
            lines!(ax, [Point3f(c1...), Point3f(c2...)];
                color=(:gray, 0.6), linewidth=1.5, linestyle=linestyle)
        end
    end

    # --- Atom sites ---
    n_types = length(unique(lat.types))
    type_colors = Makie.wong_colors()
    n_atoms_per_cell = natoms(flattice)
    n_cells = length(bravais_cells(flattice))

    for t in 1:n_types
        # atoms() orders: for each atom type, all cells, then next type
        # Actually: for each position, all bravais cells
        idxs = Int[]
        for p in 1:n_atoms_per_cell
            if lat.types[p] == t
                offset = (p - 1) * n_cells
                append!(idxs, (offset + 1):(offset + n_cells))
            end
        end
        xs = [coords[i].coords[1] for i in idxs]
        ys = [coords[i].coords[2] for i in idxs]
        zs = [coords[i].coords[3] for i in idxs]
        meshscatter!(ax, xs, ys, zs;
            color=type_colors[mod1(t, length(type_colors))],
            markersize=0.15,
            label="type $t")
    end

    display(f)
    return f
end


"""Draw a parallelepiped defined by origin + three edge vectors."""
function _draw_parallelepiped!(ax, origin, v1, v2, v3;
    color=(:blue, 0.1), linecolor=:blue, linewidth=2)

    o = origin
    # 8 corners
    c000 = o
    c100 = o .+ v1
    c010 = o .+ v2
    c001 = o .+ v3
    c110 = o .+ v1 .+ v2
    c101 = o .+ v1 .+ v3
    c011 = o .+ v2 .+ v3
    c111 = o .+ v1 .+ v2 .+ v3

    corners = [c000, c100, c010, c001, c110, c101, c011, c111]
    pts = [Point3f(c...) for c in corners]

    # 12 edges of a parallelepiped
    edges = [
        (1,2), (1,3), (1,4),    # from origin
        (2,5), (2,6),           # from c100
        (3,5), (3,7),           # from c010
        (4,6), (4,7),           # from c001
        (5,8), (6,8), (7,8)    # to c111
    ]
    for (i, j) in edges
        lines!(ax, [pts[i], pts[j]]; color=linecolor, linewidth=linewidth)
    end

    # 6 faces as triangulated quads (2 triangles per face)
    face_indices = [
        (1,2,5,3), (4,6,8,7),  # bottom, top
        (1,2,6,4), (3,5,8,7),  # front, back
        (1,3,7,4), (2,5,8,6)   # left, right
    ]
    for (a, b, c, d) in face_indices
        mesh!(ax,
            [pts[a], pts[b], pts[c], pts[d]],
            [1 2 3; 1 3 4],
            color=color, transparency=true)
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
    elseif d == 3
        plot_3d(flattice; annotate_sites, show_boundary, show_neighbors)
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