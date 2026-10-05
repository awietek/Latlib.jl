using Dates

@doc raw"""
    write_toml(flattice::FiniteLattice, opsum::OpSum, filename::String; zero_based::Bool=false, return_string::Bool=false)

Writes `OpSum` defined on a `FiniteLattice` into a TOML file
containing site coordinates and interactions.

The output format is:
```TOML
Coordinates = [
  [3.0, -3.0, 2.0],
  [2.0, -2.0, 0.0],
  ...
]

Interactions = [
  ['coupling_1', 'type_1', i_1, j_1],
  ['coupling_2', 'type_2', i_2, j_2],
  ...
]
```

The file starts with comments listing the lattice vectors, atom positions, and
boundary (torus) vectors of the finite lattice. Coordinates are obtained from
[`atoms`](@ref)`(flattice)` and rounded according to `flattice.tol`. Interactions
are listed as `[coupling, type, site1, site2]` for each operator in `opsum`.

If `symmetries` is given, a `Symmetries` section with the site permutations of the
symmetry operations of the cluster is appended (see [`toml_symmetries`](@ref)). For symmorphic
space groups, it is followed by the irreducible representations (see [`toml_irreps`](@ref) and
[`irreps`](@ref)), unless `irreps=false`. Irreducible representations are not implemented for
non-symmorphic space groups yet; in this case only the symmetry operations are written and a
warning is issued. Skipped three-dimensional irreps are reported in a banner at the top of the file.

# Arguments
- `flattice::FiniteLattice`: FiniteLattice object containing the lattice;
- `opsum::OpSum`: OpSum object containing the operators;
- `filename::String`: Full path to the output TOML file.

# Keyword arguments
- `zero_based::Bool=false`: If true, site indices are 0-based instead of 1-based.
- `return_string::Bool=false`: If true, the function returns the TOML string instead of writing to file and ignores `filename`.
- `symmetries=nothing`: `true` to determine the symmetries of `flattice` (fully periodic finite
  lattices only), or the symmetries of `flattice` as a [`FiniteSpaceGroup`](@ref) (see
  [`spacegroup`](@ref)) or a [`ClusterSymmetries`](@ref) (see [`symmetries`](@ref)).
- `irreps::Bool=true`: whether to write the irreducible representations (if `symmetries` is given).
- `origin=nothing`: symmetry center used for the irreducible representations, see
  [`symmetries`](@ref). It has to be given if the lattice has several inequivalent points with
  its full point-group symmetry, e.g. for the square lattice. Not needed for `irreps=false`.

# Examples
```julia
fl = FiniteLattice(maple_leaf, [1 1; 1 -2], true)
H = neighbor_interaction("SdotS", "J", fl)
write_toml(fl, H, "maple_leaf.toml"; zero_based=true, symmetries=true)   # with irreps

fl = FiniteLattice(square, [4 0; 0 4], true)
H = neighbor_interaction("SdotS", "J", fl)
write_toml(fl, H, "square.toml"; zero_based=true, symmetries=true, origin=LatticeVector(square, [0.0, 0.0]))
write_toml(fl, H, "square.toml"; zero_based=true, symmetries=true, irreps=false)   # permutations only
```
"""
function write_toml(flattice::FiniteLattice, opsum::OpSum, filename::String; zero_based::Bool=false, return_string::Bool=false,
                    symmetries=nothing, irreps::Bool=true, origin=nothing)
    # determine the symmetries first, so that errors are raised before anything is written
    if symmetries === true
        symmetries = spacegroup(flattice)
    elseif symmetries isa Union{FiniteSpaceGroup, ClusterSymmetries}
        group = symmetries isa FiniteSpaceGroup ? symmetries : symmetries.spacegroup
        if atoms(group.flattice) != atoms(flattice)
            throw(ArgumentError("The symmetries passed as `symmetries` belong to a different finite lattice (sites differ)."))
        end
        if symmetries isa ClusterSymmetries && !isnothing(origin)
            throw(ArgumentError("`origin` cannot be combined with a `ClusterSymmetries`, which already contains its symmetry center."))
        end
    elseif !(isnothing(symmetries) || symmetries === false)
        throw(ArgumentError("`symmetries` must be `nothing`, `true`, `false`, a `FiniteSpaceGroup` or a `ClusterSymmetries`."))
    end

    # irreducible representations: only for symmorphic space groups so far
    if symmetries isa FiniteSpaceGroup && irreps
        sg = symmetries.spacegroup
        if !sg.symmorphic
            @warn "Irreducible representations are not implemented for non-symmorphic space groups yet (the lattice has the group $(sg.symbol)). Only the symmetry operations are written."
        else
            symmetries = Latlib.symmetries(symmetries; origin=origin)
        end
    end
    sectors = symmetries isa ClusterSymmetries && irreps ? _sectors(symmetries) : nothing

    # write meta-data to TOML file
    out_str = toml_metadata()
    out_str *= "\n"
    if !isnothing(sectors) && !isempty(_skipped_labels(sectors))
        out_str *= _skipped_banner(_skipped_labels(sectors)) * "\n"
    end

    # write lattice data to TOML file
    out_str = toml_lattice(flattice; out_str=out_str)
    out_str *= "\n"
    
    # get string for `Coordinates` section of TOML file
    out_str *= toml_coordinates(flattice)
    out_str *= "\n"

    # get string for `Interactions` section of TOML file
    out_str = toml_interactions(opsum; zero_based=zero_based, out_str=out_str)

    # get string for `Symmetries` section and the irreducible representations
    if symmetries isa FiniteSpaceGroup
        out_str *= "\n"
        out_str = toml_symmetries(symmetries; zero_based=zero_based, out_str=out_str)
    elseif symmetries isa ClusterSymmetries
        out_str *= "\n"
        out_str = toml_symmetries(symmetries; zero_based=zero_based, out_str=out_str)
        if irreps
            out_str *= "\n"
            out_str = _toml_irreps(sectors; zero_based=zero_based, out_str=out_str)
        end
    end

    if return_string
        return out_str
    else
        open(filename, "w") do f
            write(f, out_str)
        end
    end
end

"""
    toml_metadata() -> String

Returns the header comment written to TOML files by [`write_toml`](@ref),
containing the date and the Latlib version.
"""
function toml_metadata() :: String
    date = Dates.today()
    latlib_version = get_latlib_version() # defined in utils.jl
    meta_str = "# This file was generated by Latlib.jl on " * string(date)
    meta_str *= " under version " * latlib_version * ".\n"
    return meta_str
end

"""
    toml_lattice(flattice::FiniteLattice; out_str::String="") -> String

Appends comments describing the lattice vectors, atom positions, and boundary vectors
of `flattice` to `out_str` and returns the result. Used by [`write_toml`](@ref).
"""
function toml_lattice(flattice; out_str::String="") :: String

    # determine rounding precision from tolerance
    digits = max(0, -floor(Int, log10(flattice.tol)))


    # lattice vectors
    out_str *= "# Lattice vectors: "
    for (i, v) in enumerate(lattice_vecs(flattice))
        rounded = [round(x; digits=digits) for x in v.coords]
        trailing = i < length(lattice_vecs(flattice)) ? ", " : ""
        out_str *= "a$i=(" * join(rounded, ", ") * ")" * trailing
    end
    out_str *= "\n"
    
    # atom positions (Cartesian basis)
    out_str *= "# Atom positions (Cartesian basis): "
    for (i, v) in enumerate(positions(flattice))
        rounded = [round(x; digits=digits) for x in to_euclidean_basis(v).coords]
        trailing = i < length(positions(flattice)) ? ", " : ""
        out_str *= "(" * join(rounded, ", ") * ")" * trailing
    end
    out_str *= "\n"

    # atom positions (lattice basis)
    out_str *= "# Atom positions (lattice basis): "
    for (i, v) in enumerate(positions(flattice))
        rounded = [round(x; digits=digits) for x in v.coords]
        trailing = i < length(positions(flattice)) ? ", " : ""
        out_str *= "(" * join(rounded, ", ") * ")" * trailing
    end
    out_str *= "\n"

    # torus vectors (Cartesian basis)
    out_str *= "# Torus vectors (Cartesian basis): "
    for (i, v) in enumerate(boundary(flattice))
        rounded = [round(x; digits=digits) for x in to_euclidean_basis(v).coords]
        trailing = i < length(boundary(flattice)) ? ", " : ""
        out_str *= "t$i=(" * join(rounded, ", ") * ")" * trailing
    end
    out_str *= "\n"

    # torus vectors (lattice basis)
    out_str *= "# Torus vectors (lattice basis): "
    for (i, v) in enumerate(boundary(flattice))
        rounded = [round(x; digits=digits) for x in v.coords] # TO-DO: rounding may be suboptimal here
        trailing = i < length(boundary(flattice)) ? ", " : ""
        out_str *= "t$i=(" * join(rounded, ", ") * ")" * trailing
    end
    out_str *= "\n"

    return out_str
end

"""
    toml_coordinates(flattice::FiniteLattice; out_str::String="") -> String

Appends the `Coordinates` section (Cartesian coordinates of all sites of `flattice`)
to `out_str` and returns the result. Used by [`write_toml`](@ref).
"""
function toml_coordinates(flattice::FiniteLattice; out_str::String="") :: String
    
    # determine rounding precision from tolerance
    digits = max(0, -floor(Int, log10(flattice.tol)))

    # get site coordinates as Vector{EuclideanVector}
    euc_coords = atoms(flattice)

    # write coordinates
    out_str *= "Coordinates = [\n"
    for (i, v) in enumerate(euc_coords)
        rounded = [round(x; digits=digits) for x in v.coords]
        trailing = ","
        out_str *= "  [" * join(rounded, ", ") * "]" * trailing * "\n"
    end
    out_str *= "]\n"
    return out_str
end


@doc raw"""
    toml_symmetries(cs::ClusterSymmetries; zero_based::Bool=false, out_str::String="") -> String
    toml_symmetries(g::FiniteSpaceGroup; zero_based::Bool=false, out_str::String="") -> String

Appends the `Symmetries` section to `out_str` and returns the result. Used by [`write_toml`](@ref).

The section lists the site permutations of the symmetry operations of the finite lattice:
the entry `[p_1, p_2, ...]` maps site `i` onto site `p_i`. Comments before the section state
the symmetry group of the cluster and of the infinite lattice. For a [`ClusterSymmetries`](@ref),
they also state the symmetry center and list all momenta with their labels and little
co-groups, marking the representatives of the stars with `*`; these are the momenta of the
irreducible representations written by [`toml_irreps`](@ref).

```TOML
Symmetries = [
  [0, 1, 2, 3, ...],
  [1, 0, 3, 2, ...],
  ...
]
```

Only one operation per distinct site permutation is written (see [`distinct_operations`](@ref)).
If operations other than the identity act trivially on the sites (see [`trivial_operations`](@ref)),
which happens on small or thin clusters, the omitted operations are reported in a comment block.
The format of this report is a placeholder and may change.

Site indices follow `zero_based` as in [`toml_interactions`](@ref). Codes such as
[XDiag](https://github.com/awietek/xdiag) expect 0-based indices, so a warning is issued for
`zero_based=false`.
"""
function toml_symmetries(g::FiniteSpaceGroup; zero_based::Bool=false, out_str::String="") :: String
    zero_based || @warn "Symmetries are written with 1-based site indices, but codes like XDiag expect 0-based indices. Use `zero_based=true`."
    offset = zero_based ? 1 : 0
    sg = g.spacegroup

    out_str *= @sprintf "# Symmetry group of the cluster: %s (#%d), point group %s (%s), %d operations\n" g.symbol g.number g.pointgroup g.schoenflies length(g)
    out_str *= @sprintf "# Symmetry group of the infinite lattice: %s (#%d), point group %s (%s)\n" sg.symbol sg.number sg.pointgroup sg.schoenflies
    if length(trivial_operations(g)) > 1
        out_str *= toml_omitted_symmetries_placeholder(g)
    end
    return _toml_permutations(g.permutations[distinct_operations(g)], offset, out_str)
end

function toml_symmetries(cs::ClusterSymmetries; zero_based::Bool=false, out_str::String="") :: String
    zero_based || @warn "Symmetries are written with 1-based site indices, but codes like XDiag expect 0-based indices. Use `zero_based=true`."
    offset = zero_based ? 1 : 0
    g = cs.spacegroup
    sg = g.spacegroup
    lattice = g.flattice.lattice

    out_str *= "# Symmetry center: " * _describe_point(lattice, cs.origin) * "\n"
    out_str *= @sprintf "# Symmetry group of the cluster: %s (#%d), point group %s (%s), %d operations\n" g.symbol g.number g.pointgroup g.schoenflies length(g)
    if length(cs.operations) < length(g)
        out_str *= @sprintf "# Only the %d operations with translations by lattice vectors are used (the unit cell is not primitive).\n" length(cs.operations)
    end
    out_str *= @sprintf "# Symmetry group of the infinite lattice: %s (#%d), point group %s (%s)\n" sg.symbol sg.number sg.pointgroup sg.schoenflies

    # all momenta, grouped by stars, representatives first
    ks = momenta(cs)
    reps = sort([k for k in ks if k.representative]; by=k -> _kpoint_sortkey(k.label))
    out_str *= "# Momenta (Cartesian, first Brillouin zone; representatives of the stars marked with *):\n"
    for r in reps, k in sort([k for k in ks if k.star == r.star]; by=k -> !k.representative)
        coords = join([@sprintf("%.10f", x) for x in k.momentum], ", ")
        out_str *= @sprintf "#   %-8s %-5s (%s)%s\n" k.label k.littlegroup_name coords (k.representative ? " *" : "")
    end

    if length(cs.trivial) > 1
        ops = operations(g)
        trivial = [ops[cs.operations[j]] for j in cs.trivial[2:end]]
        out_str *= toml_omitted_symmetries_placeholder(trivial, length(cs.operations) - length(cs.permutations), length(cs.operations))
    end
    return _toml_permutations(cs.permutations, offset, out_str)
end

function _toml_permutations(permutations::Vector{Vector{Int}}, offset::Int, out_str::String) :: String
    out_str *= "Symmetries = [\n"
    for p in permutations
        out_str *= "  [" * join(p .- offset, ", ") * "],\n"
    end
    out_str *= "]\n"
    return out_str
end

# PLACEHOLDER for omitted symmetry operations.
# Operations that act on the sites like another operation (trivially acting operations
# on small or thin clusters) are left out of `Symmetries`, since all irreducible
# representations that are non-trivial on them vanish. How this should be stated in the
# TOML file is to be decided (to be agreed upon with XDiag). Until then, a comment block
# with the marker "PLACEHOLDER" reports the omission.
function toml_omitted_symmetries_placeholder(trivial::Vector{SymmetryOperation}, nomitted::Int, ntotal::Int) :: String
    s  = "# PLACEHOLDER (TOML format to be decided): omitted symmetry operations\n"
    s *= @sprintf "# %d of the %d symmetry operations of this cluster are not listed in `Symmetries`, because\n" nomitted ntotal
    s *= "# they permute the sites in the same way as another operation. The following operations\n"
    s *= "# act trivially on the sites (coordinates in the lattice basis):\n"
    for op in trivial
        s *= "#   " * _xyz_string(op) * "\n"
    end
    s *= "# Irreducible representations that are not trivial on these operations vanish on this cluster.\n"
    return s
end

function toml_omitted_symmetries_placeholder(g::FiniteSpaceGroup) :: String
    trivial = g.operations[trivial_operations(g)[2:end]]
    return toml_omitted_symmetries_placeholder(trivial, length(g) - length(distinct_operations(g)), length(g))
end

# PLACEHOLDER for omitted irreducible representations (see above).
function toml_omitted_irreps_placeholder(labels::Vector{String}) :: String
    s  = "# PLACEHOLDER (TOML format to be decided): omitted irreducible representations\n"
    s *= "# The following irreducible representations vanish on this cluster, because they are not\n"
    s *= "# trivial on operations that act trivially on the sites:\n"
    for label in labels
        s *= "#   " * label * "\n"
    end
    return s
end

@doc raw"""
    toml_irreps(cs::ClusterSymmetries; zero_based::Bool=false, out_str::String="") -> String

Appends the irreducible representations of the symmetry group of a finite lattice (see
[`irreps`](@ref)) to `out_str` and returns the result. Used by [`write_toml`](@ref).

Every irreducible representation is written as a table named after its label, containing the
characters (as `[real, imag]` pairs) of the allowed symmetries, the indices of the allowed
symmetries in the `Symmetries` section, and the Cartesian momentum:

```TOML
[K.C3.Ea]
characters = [
  [1.0000000000000000, 0.0000000000000000],
  [-0.4999999999999998, -0.8660254037844388],
  ...
]
allowed_symmetries = [0, 1, 2, 3, 4, 5, 6, 7, 8]
momentum = [1.5832138822983011, 0.0000000000000000]
```

The indices of the allowed symmetries follow `zero_based`. Irreducible representations that
vanish on the cluster are left out and listed in a comment block; its format is a
placeholder and may change. Skipped three-dimensional irreps are listed in a banner, and a
warning is issued.
"""
function toml_irreps(cs::ClusterSymmetries; zero_based::Bool=false, out_str::String="") :: String
    return _toml_irreps(_sectors(cs); zero_based=zero_based, out_str=out_str)
end

function _toml_irreps(sectors; zero_based::Bool=false, out_str::String="") :: String
    offset = zero_based ? 1 : 0
    out_str *= "# Irreducible representations\n"
    skipped = _skipped_labels(sectors)
    if !isempty(skipped)
        _warn_skipped(skipped)
        out_str *= _skipped_banner(skipped)
    end
    vanishing = [irrep.label for (irrep, status) in sectors if status == :vanishing]
    isempty(vanishing) || (out_str *= toml_omitted_irreps_placeholder(vanishing))
    for (irrep, status) in sectors
        status == :ok || continue
        if irrep.dimension == 2 && endswith(irrep.name, "a")
            out_str *= "# $(irrep.kpoint).$(irrep.littlegroup).$(irrep.parent) is two-dimensional: its partners $(irrep.parent)a and $(irrep.parent)b\n"
            out_str *= "# are sectors of a subgroup of the little group and are exactly degenerate.\n"
        end
        out_str *= "[" * irrep.label * "]\n"
        out_str *= "characters = [\n"
        for χ in irrep.characters
            out_str *= @sprintf "  [%.16f, %.16f],\n" real(χ) imag(χ)
        end
        out_str *= "]\n"
        out_str *= "allowed_symmetries = [" * join(irrep.allowed_symmetries .- offset, ", ") * "]\n"
        out_str *= "momentum = [" * join([@sprintf("%.16f", x) for x in irrep.momentum], ", ") * "]\n"
        out_str *= "\n"
    end
    return out_str
end

# banner for skipped three-dimensional irreps
function _skipped_banner(skipped::Vector{String}) :: String
    line = "# " * "!"^90 * "\n"
    s  = line
    s *= "# !!! WARNING: three-dimensional irreducible representations were SKIPPED, since they cannot be\n"
    s *= "# !!! written as one-dimensional characters. The irreducible representations in this file\n"
    s *= "# !!! do NOT span the full Hilbert space. Skipped:\n"
    for label in skipped
        s *= "# !!!   " * label * "\n"
    end
    return s * line
end

"""
    toml_interactions(opsum::OpSum; zero_based::Bool=false, out_str::String="") -> String

Appends the `Interactions` section (one entry `[coupling, type, site1, site2]` per operator
in `opsum`) to `out_str` and returns the result. Used by [`write_toml`](@ref).
"""
function toml_interactions(opsum::OpSum; zero_based::Bool=false, out_str::String="") :: String
    
    # 0 or 1-based indexing?
    offset = zero_based ? 1 : 0 
    
    # write interactions
    out_str *= "Interactions = [\n"
    for (i, op) in enumerate(opsum.ops)
        s1 = op.sites[1] - offset
        s2 = op.sites[2] - offset
        trailing = ","
        out_str *= "  ['" * op.cpl * "', '" * op.type * "', " * string(s1) * ", " * string(s2) * "]" * trailing * "\n"
    end
    out_str *= "]\n"
    return out_str
end


