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
- `annotate_sites::Bool=false`: label each site with its number index in atoms(flattice).
- `annotate_sites_zero_based::Bool=true`: Determines if site annotations are 0-based or 1-based.
- `show_boundary::Bool=true`: show the boundary parallelepiped of the finite lattice
- `show_unit_cell::Bool=false`: show the unit cell parallelepiped at the origin
- `show_bravais_grid::Bool=false`: show Bravais lattice vectors as arrows at each cell origin
- `show_neighbors::Bool=true`: show nearest-neighbor bonds (usually correspond to the lattice edges).
- `site_marksize::Float64=0.15`: Controls the size of the spheres representing lattice sites.
- `draw_periodic_flattice::Bool=false`: show grey copies of finite lattice for each face of the finite lattice boundary.
- `draw_periodic_flattice_shifts::Vector{Tuple{Int, Int, Int}}=nothing`: when `draw_periodic_flattice=true`, determines 
    which periodic images to draw. Each tuple corresponds to a shift along the three boundary vectors of the finite lattice.
    By default, all 6 face-sharing neighbors are drawn: `[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]`.
- `scale_factor::Number=1.5`: controls, among other things, the scale of labels
    """
function plot_3d(flattice::FiniteLattice;
    annotate_sites::Bool=false,
    annotate_sites_zero_based::Bool=true,
    show_boundary::Bool=true,
    show_unit_cell::Bool=false,
    show_bravais_grid::Bool=false,
    show_neighbors::Bool=true,
    site_marksize::Float64=0.15,
    draw_periodic_flattice::Bool=false,
    draw_periodic_flattice_shifts::Union{Nothing, Vector{Tuple{Int, Int, Int}}}=nothing,
    scale_factor::Number=1.5,
    )

    if dim(flattice) != 3
        error(@sprintf "plot_3d only supports 3D lattices, but got dimension %d" dim(flattice))
    end

    GLMakie.activate!(ssao=true, fxaa=true, scalefactor=scale_factor)

    lat = flattice.lattice
    As = lattice_vecs(flattice)
    a1 = As[1].coords # converts EuclideanVector to Vector{Float64}
    a2 = As[2].coords
    a3 = As[3].coords

    f = Figure(size=(1000, 800), figure_padding=20)
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
        arrows3d!(ax, ox, oy, oz, dx, dy, dz;
            color=bravais_arrow_color,
            tipradius=0.15,
            tiplength=0.2,
            shaftradius=0.04)
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
            color=(:orange, 0.06), linecolor=:orange, linewidth=1.5,
            extend_flattice_cell_edges=0.2)
    end

    # --- Neighbor bonds ---
    coords = atoms(flattice)
    if show_neighbors
        nbors = neighbors(flattice)
        metric_pbc = PeriodicEuclideanMetric(flattice)
        metric_euc = EuclideanMetric()

        # Precompute per-bond: is_periodic flag and (for periodic bonds) the shift
        # that places p1 next to p2
        bond_data = map(nbors) do (i, j)
            p1 = coords[i]; p2 = coords[j]
            is_periodic = !isapprox(metric_pbc(p1, p2), metric_euc(p1, p2))
            if is_periodic
                d = distance_vector(p2, p1; flattice=flattice)
                p1_shifted = p2 + d
                (i, j, true, p1_shifted.coords, p2.coords)
            else
                (i, j, false, p1.coords, p2.coords)
            end
        end

        # Determine which cell offsets to draw bonds for
        if draw_periodic_flattice
            bvecs_bond = [to_euclidean_basis(v).coords for v in boundary(flattice)]
            resolved_shifts = if isnothing(draw_periodic_flattice_shifts)
                [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
            else
                draw_periodic_flattice_shifts
            end
            cell_offsets = [[0.0, 0.0, 0.0]]
            for (n1, n2, n3) in resolved_shifts
                push!(cell_offsets, n1 .* bvecs_bond[1] .+ n2 .* bvecs_bond[2] .+ n3 .* bvecs_bond[3])
            end
        else
            cell_offsets = [[0.0, 0.0, 0.0]]
        end

        for offset in cell_offsets
            is_main = all(offset .== 0.0)
            for (i, j, is_periodic, c1, c2) in bond_data
                linestyle = (is_periodic && !draw_periodic_flattice) ? :dash : :solid
                bond_color = is_main ? (:gray, 0.8) : (:gray, 0.6)
                pa = Point3f((c1 .+ offset)...)
                pb = Point3f((c2 .+ offset)...)
                lines!(ax, [pa, pb]; color=bond_color, linewidth=1.5, linestyle=linestyle)
            end
        end
    end

    # --- Atom sites ---
    n_types = length(unique(lat.types))
    type_colors = Makie.wong_colors()
    n_atoms_per_cell = natoms(flattice)
    n_cells = length(bravais_cells(flattice))

    for t in 1:n_types
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
            markersize=site_marksize,
            label="type $t")
    end

    # --- Site index labels ---
    if annotate_sites
        positions = [Point3f(c.coords...) for c in coords]
        if annotate_sites_zero_based
            labels = string.(0:(length(coords)-1))
        else
            labels = string.(1:length(coords))
        end
        text!(ax, positions; text=labels, fontsize=14, color=:black, offset=(5, 5))
    end

    # --- Shadow cells: grey copies in all 26 neighbouring periodic images ---
    if draw_periodic_flattice
        bvecs = [to_euclidean_basis(v).coords for v in boundary(flattice)]
        base_positions = [c.coords for c in coords]
        shadow_color = (:grey, 0.35)

        if annotate_sites
            if annotate_sites_zero_based
                shadow_labels = string.(0:(length(coords)-1))
            else
                shadow_labels = string.(1:length(coords))
            end
        end

        
        if isnothing(draw_periodic_flattice_shifts)
            # 6 face-sharing neighbours: ±1 along each boundary vector
            draw_periodic_flattice_shifts = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
        end
        for (n1, n2, n3) in draw_periodic_flattice_shifts
            shift = n1 .* bvecs[1] .+ n2 .* bvecs[2] .+ n3 .* bvecs[3]
            sx = [p[1] + shift[1] for p in base_positions]
            sy = [p[2] + shift[2] for p in base_positions]
            sz = [p[3] + shift[3] for p in base_positions]
            meshscatter!(ax, sx, sy, sz;
                color=shadow_color, markersize=0.12)
            if annotate_sites
                spts = [Point3f(sx[k], sy[k], sz[k]) for k in eachindex(sx)]
                text!(ax, spts; text=shadow_labels,
                    fontsize=12, color=(:grey, 0.8), offset=(5, 5))
            end
        end
    end

    display(f)
    return f, ax
end


@doc raw"""
    plot_3d(flattice, opsum; kwargs...)

Plots a 3D `FiniteLattice` with colored interaction edges from an `OpSum`.

# Arguments:
- `flattice::FiniteLattice`: a 3D finite lattice to plot;
- `opsum::OpSum`: an operator sum containing interactions to plot as edges;

# Keyword arguments:
- `cpl_dict::Dict`: a dictionary mapping coupling strings to colors;
- `draw_periodic_flattice::Bool=false`: see plot_3d(flattice)
- `draw_periodic_flattice_shifts::Vector{Tuple{Int, Int, Int}}=nothing`: see plot_3d(flattice)
- all other arguments are forwarded to `plot_3d(flattice; ...)`.
"""
function plot_3d(flattice::FiniteLattice, opsum::OpSum;
    cpl_dict::Dict=nothing,
    draw_periodic_flattice::Bool=false,
    draw_periodic_flattice_shifts::Union{Nothing, Vector{Tuple{Int, Int, Int}}}=nothing,
    kwargs...)

    f, ax = plot_3d(flattice;
        show_neighbors=false,
        draw_periodic_flattice=draw_periodic_flattice,
        draw_periodic_flattice_shifts=draw_periodic_flattice_shifts,
        kwargs...)

    coords = atoms(flattice)
    metric_pbc = PeriodicEuclideanMetric(flattice)
    metric_euc = EuclideanMetric()

    # if cpl_dict is unspecified, assign a distinct color to each unique coupling string
    if isnothing(cpl_dict)
        couplings = unique(op.cpl for op in opsum.ops)
        cpl_colors = Dict(cpl => Makie.wong_colors()[mod1(i, length(Makie.wong_colors()))] # TO-DO: handle more than 10 unique couplings here!!!
                           for (i, cpl) in enumerate(couplings))
    else
        cpl_colors = cpl_dict
    end

    # Determine cell offsets for interaction edges
    if draw_periodic_flattice
        bvecs = [to_euclidean_basis(v).coords for v in boundary(flattice)]
        resolved_shifts = if isnothing(draw_periodic_flattice_shifts)
            [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
        else
            draw_periodic_flattice_shifts
        end
        cell_offsets = [[0.0, 0.0, 0.0]]
        for (n1, n2, n3) in resolved_shifts
            push!(cell_offsets, n1 .* bvecs[1] .+ n2 .* bvecs[2] .+ n3 .* bvecs[3])
        end
    else
        cell_offsets = [[0.0, 0.0, 0.0]]
    end

    # Precompute bond endpoints (with periodic wrapping resolved)
    op_data = map(opsum.ops) do op
        site1, site2 = op.sites
        p1 = coords[site1]; p2 = coords[site2]
        is_periodic = !isapprox(metric_pbc(p1, p2), metric_euc(p1, p2))
        if is_periodic
            d = distance_vector(p2, p1; flattice=flattice)
            c1 = (p2 + d).coords
        else
            c1 = p1.coords
        end
        (op.cpl, is_periodic, c1, p2.coords)
    end

    # draw Ops as colored edges with color given by coupling
    for offset in cell_offsets
        # draw interactions in main cell with full opacity, and lower it otherwise
        edge_alpha = 0.6
        if isapprox(offset, [0.0, 0.0, 0.0])
            edge_alpha = 1.0
        end
        # start printing loop
        for (cpl, is_periodic, c1, c2) in op_data
            linestyle = (is_periodic && !draw_periodic_flattice) ? :dash : :solid
            pa = Point3f((c1 .+ offset)...)
            pb = Point3f((c2 .+ offset)...)
            if cpl ∉ keys(cpl_colors)
                # couplings that are not specified in cpl_dict are not drawn at all
                #println("Warning: coupling string \"$cpl\" not found in cpl_dict, skipping drawing these bonds.")
                #println("Available couplings in cpl_dict: ", keys(cpl_colors))
                continue
            end
            lines!(ax, [pa, pb]; color=(cpl_colors[cpl], edge_alpha), linewidth=2.5,
                linestyle=linestyle, label=string(cpl))
        end
    end

    # --- Legend panel on the right side ---
    legend_entries = [LineElement(color=cpl_colors[cpl], linewidth=3)
                      for cpl in keys(cpl_colors)]
    legend_labels  = [string(cpl) for cpl in keys(cpl_colors)]
    Legend(f[1, 2], legend_entries, legend_labels;
        framevisible=true, labelsize=14, patchsize=(20, 3))
    # make sure the 3D scene takes most of the width
    colsize!(f.layout, 1, Relative(0.85))

    display(f)
    return f, ax
end


@doc raw"""
    plot_3d(flattice, opsum, spin_data; kwargs...)

Plots a 3D `FiniteLattice` with colored interaction edges from an `OpSum`
and 3D spin arrows at each lattice site.

# Arguments:
- `flattice::FiniteLattice`: a 3D finite lattice to plot;
- `opsum::OpSum`: an operator sum containing interactions to plot as edges;
- `spin_data::AbstractVector`: a length-N vector of 3-component vectors giving
    the spin at each lattice site. The k-th entry corresponds to the k-th site
    in `atoms(flattice)`.

# Keyword arguments:
- `spin_length_multiplier::Float64=1.0`: multiplier for the length of the spin arrows, to match with lattice spacings.;
- `spin_color:gray`: color of the spin arrows;
- `spin_tipradius::Float64=0.12`: radius of the arrowhead cone;
- `spin_tiplength::Float64=0.16`: length of the arrowhead cone;
- `spin_shaftradius::Float64=0.06`: radius of the arrow shaft;
- all other keyword arguments are forwarded to `plot_3d(flattice, opsum; ...)`.
"""
function plot_3d(flattice::FiniteLattice, opsum::OpSum, spin_data::AbstractVector;
    spin_length_multiplier::Float64 = 1.0,
    spin_color=:gray,
    spin_tipradius::Float64=0.1,
    spin_tiplength::Float64=0.2,
    spin_shaftradius::Float64=0.04,
    draw_periodic_flattice::Bool=false,
    draw_periodic_flattice_shifts::Union{Nothing, Vector{Tuple{Int, Int, Int}}}=nothing,
    kwargs...)

    f, ax = plot_3d(flattice, opsum;
        draw_periodic_flattice=draw_periodic_flattice,
        draw_periodic_flattice_shifts=draw_periodic_flattice_shifts,
        kwargs...)

    coords = atoms(flattice)
    N = length(coords)

    if length(spin_data) != N
        error("spin_data has length $(length(spin_data)) but the finite lattice has $N sites")
    end

    # Arrow origins (centered at lattice site) and directions
    ox = Float64[]; oy = Float64[]; oz = Float64[]
    dx = Float64[]; dy = Float64[]; dz = Float64[]
    for k in 1:N
        s = spin_data[k] * spin_length_multiplier
        c = coords[k].coords
        # center the arrow at the site: origin = site - spin/2
        push!(ox, c[1] - s[1] / 2)
        push!(oy, c[2] - s[2] / 2)
        push!(oz, c[3] - s[3] / 2)
        push!(dx, s[1]); push!(dy, s[2]); push!(dz, s[3])
    end

    arrows3d!(ax, ox, oy, oz, dx, dy, dz;
        color=spin_color,
        tipradius=spin_tipradius*spin_length_multiplier,
        tiplength=spin_tiplength*spin_length_multiplier,
        shaftradius=spin_shaftradius*spin_length_multiplier)

    # Draw spin arrows on periodic copies as well
    if draw_periodic_flattice
        bvecs = [to_euclidean_basis(v).coords for v in boundary(flattice)]
        resolved_shifts = if isnothing(draw_periodic_flattice_shifts)
            [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
        else
            draw_periodic_flattice_shifts
        end
        for (n1, n2, n3) in resolved_shifts
            shift = n1 .* bvecs[1] .+ n2 .* bvecs[2] .+ n3 .* bvecs[3]
            sox = ox .+ shift[1]
            soy = oy .+ shift[2]
            soz = oz .+ shift[3]
            arrows3d!(ax, sox, soy, soz, dx, dy, dz;
                color=(spin_color, 0.35),
                tipradius=spin_tipradius*spin_length_multiplier,
                tiplength=spin_tiplength*spin_length_multiplier,
                shaftradius=spin_shaftradius*spin_length_multiplier)
        end
    end

    display(f)
    return f, ax
end


"""Draw a parallelepiped defined by origin + three edge vectors.

When `extend_flattice_cell_edges > 0`, each edge is drawn longer than the actual
parallelepiped by that fraction of its length in both directions,
giving a visual guide for translating the cell."""
function _draw_parallelepiped!(ax, origin, v1, v2, v3;
    color=(:blue, 0.1), linecolor=:blue, linewidth=2,
    extend_flattice_cell_edges::Real=0.0)

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
        p1 = pts[i]; p2 = pts[j]
        if extend_flattice_cell_edges > 0
            d = p2 - p1
            p1 = p1 - extend_flattice_cell_edges * d
            p2 = p2 + extend_flattice_cell_edges * d
        end
        lines!(ax, [p1, p2]; color=linecolor, linewidth=linewidth)
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