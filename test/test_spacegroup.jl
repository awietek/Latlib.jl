using LinearAlgebra
using TOML

# This file tests the spglib interface in src/symmetry/spacegroup.jl.
#
# The expected groups of the infinite lattices are the textbook ones. The expected
# groups of the finite clusters were derived by hand: a rotation of the infinite
# lattice survives on the cluster iff it maps the torus vectors onto torus vectors
# (comments below give the argument). In addition, all symmetry operations are
# cross-checked against an independent brute-force search, and the site permutations
# are verified geometrically.


# ================================================================
# Independent brute-force reference (test only)
# ================================================================

# whether d (lattice basis) is a Bravais lattice vector, up to a Cartesian distance tol
bf_isbravais(lattice, d; tol=1e-6) = norm(lattice.A' * (d - round.(d))) < tol

# all integer matrices W (entries in -3:3) that are isometries of the Bravais lattice
function bf_isometries(lattice; tol=1e-8)
    D = dim(lattice)
    G = lattice.A * lattice.A'    # Gram matrix of the lattice vectors
    box = [collect(v) for v in Iterators.product(ntuple(_ -> -3:3, D)...)]
    columns = [[v for v in box if abs(dot(v, G * v) - G[j, j]) < tol] for j in 1:D]
    Ws = Matrix{Int}[]
    for cols in Iterators.product(columns...)
        W = reduce(hcat, cols)
        norm(W' * G * W - G) < tol && push!(Ws, W)
    end
    return Ws
end

# all space-group operations (W, w) of the infinite lattice, w modulo Bravais translations
function bf_operations(lattice)
    P = natoms(lattice)
    X = [lattice.positions[a, :] for a in 1:P]
    is_symmetry(W, w) = all(any(lattice.types[b] == lattice.types[a] &&
                                bf_isbravais(lattice, W * X[a] + w - X[b]) for b in 1:P) for a in 1:P)
    ops = Tuple{Matrix{Int}, Vector{Float64}}[]
    for W in bf_isometries(lattice), b in 1:P
        lattice.types[b] == lattice.types[1] || continue
        w = X[b] - W * X[1]
        if is_symmetry(W, w) && !any(op -> op[1] == W && bf_isbravais(lattice, op[2] - w), ops)
            push!(ops, (W, w))
        end
    end
    return ops
end

# whether the rotational part W maps the torus of the finite lattice onto itself
function bf_torus_compatible(flattice, W)
    T = Matrix{Float64}(flattice.boundary')
    M = T \ (W * T)
    return norm(M - round.(M)) < 1e-8
end

# whether the operations of sg coincide with the brute-force reference ref
function same_operations(sg::SpaceGroup, ref)
    length(operations(sg)) == length(ref) || return false
    return all(any(op.W == W && bf_isbravais(sg.lattice, op.w - w) for (W, w) in ref) for op in operations(sg))
end

# fractional coordinates and atom types of all sites of a finite lattice
function site_data(flattice)
    lattice = flattice.lattice
    X = [lattice.A' \ x.coords for x in atoms(flattice)]
    types = [lattice.types[findfirst(a -> bf_isbravais(lattice, x - lattice.positions[a, :]), 1:natoms(lattice))] for x in X]
    return X, types
end

# whether two points (lattice basis) coincide modulo the torus of the finite lattice
function same_site(flattice, x, y)
    T = Matrix{Float64}(flattice.boundary')
    u = T \ (x - y)
    return norm(flattice.lattice.A' * T * (u - round.(u))) < 1e-6
end

# every operation maps every site i onto the site perm[i] of the same type
function permutations_are_correct(g::FiniteSpaceGroup)
    X, types = site_data(g.flattice)
    for (op, p) in zip(operations(g), site_permutations(g))
        isperm(p) || return false
        for i in eachindex(X)
            types[p[i]] == types[i] || return false
            same_site(g.flattice, op(X[i]), X[p[i]]) || return false
        end
    end
    return true
end

# The permutations form a group: identity first, closed under composition, and every
# permutation appears equally often. (On small clusters several operations can act
# identically on the sites; they form a normal subgroup.)
function permutations_form_group(g::FiniteSpaceGroup; maxpairs=40_000)
    perms = site_permutations(g)
    perms[1] == collect(1:length(perms[1])) || return false
    counts = Dict{Vector{Int}, Int}()
    for p in perms
        counts[p] = get(counts, p, 0) + 1
    end
    length(unique(values(counts))) == 1 || return false
    distinct = collect(keys(counts))
    n = length(distinct)
    for k in 0:max(1, (n * n) ÷ maxpairs):(n * n - 1)
        p, q = distinct[k ÷ n + 1], distinct[k % n + 1]
        haskey(counts, p[q]) || return false   # composition: first q, then p
    end
    return true
end


# ================================================================
# Test lattices
# ================================================================

# (name, lattice, symbol, ITA number, point group, |point group|, symmorphic, pure translations per cell)
const LATTICE_CASES = [
    ("square",             square,             "p4mm",     11,  "C4v", 8,  true,  1),
    ("triangular",         triangular,         "p6mm",     17,  "C6v", 12, true,  1),
    ("honeycomb",          honeycomb,          "p6mm",     17,  "C6v", 12, true,  1),
    ("kagome",             kagome,             "p6mm",     17,  "C6v", 12, true,  1),
    ("lieb",               lieb,               "p4mm",     11,  "C4v", 8,  true,  1),
    ("trellis",            trellis,            "c2mm",     9,   "C2v", 4,  true,  1),
    ("maple_leaf",         maple_leaf,         "p6",       16,  "C6",  6,  true,  1),
    # sites only form a square lattice with half the lattice constant: 4 translations per cell
    ("shastry_sutherland", shastry_sutherland, "p4mm",     11,  "C4v", 8,  true,  4),
    # orthogonal-dimer geometry: the sites alone have the non-symmorphic symmetry of the model
    ("shastry_sutherland_non_symmorphic", shastry_sutherland_non_symmorphic, "p4gm", 12, "C4v", 8, false, 1),
    ("simple_cubic",       simple_cubic,       "Pm-3m",    221, "Oh",  48, true,  1),
    ("bcc",                bcc,                "Im-3m",    229, "Oh",  48, true,  1),
    ("fcc",                fcc,                "Fm-3m",    225, "Oh",  48, true,  1),
    ("simple_hexagonal",   simple_hexagonal,   "P6/mmm",   191, "D6h", 24, true,  1),
    ("hcp",                hcp,                "P6_3/mmc", 194, "D6h", 24, false, 1),
    ("diamond",            diamond,            "Fd-3m",    227, "Oh",  48, false, 1),
    ("pyrochlore",         pyrochlore,         "Fd-3m",    227, "Oh",  48, false, 1),
    ("hyperhoneycomb",     hyperhoneycomb,     "Fddd",     70,  "D2h", 8,  false, 1),
]

# (name, boundary, symbol, |point group|, |group|), with |group| = |point group| x N_cells x translations per cell
const CLUSTER_CASES = Dict(
    "square" => [
        ([4 0; 0 4],  "p4mm", 8, 128),   # torus invariant under the full 4mm
        ([4 0; 0 2],  "p2mm", 4, 32),    # rectangular torus: no 4-fold axis, no diagonal mirrors
        ([3 1; -1 3], "p4",   4, 40),    # tilted torus: C4 maps (3,1) -> (-1,3), mirrors do not preserve it
    ],
    "triangular" => [
        ([3 0; 0 3],  "p6mm", 12, 108),  # 3x the Bravais lattice: all operations survive
        ([2 0; 0 1],  "c2mm", 4, 8),     # only C2 and the mirrors along/perpendicular to a2
        ([3 1; -1 4], "p6",   6, 78),    # C6 maps (3,1) -> (-1,4); chiral torus, no mirrors
    ],
    "honeycomb" => [
        ([3 0; 0 3],  "p6mm", 12, 108),
        ([2 -2; 1 1], "c2mm", 4, 16),    # rectangular torus (0,2), (√3,0): mirrors along x and y
        ([3 1; 4 -3], "p6",   6, 78),    # C6 maps (3,1) -> (4,-3) in the honeycomb basis
    ],
    "kagome" => [
        ([2 0; 0 2],  "p6mm", 12, 48),
        ([2 0; 0 1],  "c2mm", 4, 8),
        ([3 1; -1 4], "p6",   6, 78),
    ],
    "lieb" => [
        ([2 0; 0 2], "p4mm", 8, 32),
        ([2 0; 0 1], "p2mm", 4, 8),
    ],
    "trellis" => [
        ([4 0; -2 4], "c2mm", 4, 64),    # torus along the rectangular axes: (4,0) and (0, 4 + 2√3)
        ([3 0; 1 2],  "p2",   2, 12),    # oblique torus: only C2 survives
    ],
    "maple_leaf" => [
        ([1 1; 1 -2], "p6", 6, 18),      # legacy cluster 18.v1
        ([1 0; 0 2],  "p2", 2, 4),       # legacy cluster 12.v1
    ],
    "shastry_sutherland" => [
        ([2 0; 0 2], "p4mm", 8, 128),
        ([2 0; 0 1], "p2mm", 4, 32),
    ],
    "shastry_sutherland_non_symmorphic" => [
        ([2 0; 0 2],  "p4gm", 8, 32),
        ([2 0; 0 1],  "p2gg", 4, 8),     # axial reflections of p4gm are glides
        ([1 1; 2 -2], "c2mm", 4, 16),    # diagonal reflections of p4gm are mirrors: symmorphic subgroup
    ],
    "simple_cubic" => [
        ([2 0 0; 0 2 0; 0 0 2], "Pm-3m",  48, 384),
        ([2 0 0; 0 2 0; 0 0 3], "P4/mmm", 16, 192),
        ([1 0 0; 0 2 0; 0 0 3], "Pmmm",   8,  48),
    ],
    "bcc" => [
        ([0 1 1; 1 0 1; 1 1 0], "Im-3m",  48, 96),   # conventional cubic cell
        ([0 1 1; 1 0 1; 2 2 0], "I4/mmm", 16, 64),   # conventional cell doubled along z
    ],
    "fcc" => [
        ([-1 1 1; 1 -1 1; 1 1 -1], "Fm-3m",  48, 192),
        ([-1 1 1; 1 -1 1; 2 2 -2], "I4/mmm", 16, 128),
    ],
    "simple_hexagonal" => [
        ([3 0 0; 0 3 0; 0 0 2], "P6/mmm", 24, 432),
        ([2 0 0; 0 1 0; 0 0 1], "Cmmm",   8,  16),   # in-plane: mirrors along/perpendicular to a2 only
        ([3 1 0; 2 3 0; 0 0 1], "P6/m",   12, 84),   # chiral in-plane torus
    ],
    "hcp" => [
        ([3 0 0; 0 3 0; 0 0 2], "P6_3/mmc", 24, 432),
        ([2 0 0; 0 1 0; 0 0 1], "Cmcm",     8,  16),
    ],
    "diamond" => [
        ([-1 1 1; 1 -1 1; 1 1 -1], "Fd-3m",    48, 192),   # conventional cubic cell, 8 sites
        ([-1 1 1; 1 -1 1; 2 2 -2], "I4_1/amd", 16, 128),   # tetragonal subgroup
        ([-1 1 1; 2 -2 2; 3 3 -3], "Fddd",     8,  192),   # orthorhombic subgroup
    ],
    "pyrochlore" => [
        ([-1 1 1; 1 -1 1; 1 1 -1], "Fd-3m",    48, 192),
        ([-1 1 1; 1 -1 1; 2 2 -2], "I4_1/amd", 16, 128),
    ],
    "hyperhoneycomb" => [
        ([-1 1 1; 0 1 -1; -1 1 -1], "Fddd", 8, 16),   # 8-site cluster of examples/Hyperhoneycomb
        ([-1 1 1; 0 1 -1; -2 2 -2], "C2/c", 4, 16),   # 16-site cluster: monoclinic subgroup
    ],
)

const NONSYMMORPHIC_SYMBOLS = Set(["p4gm", "p2gg", "P6_3/mmc", "Cmcm", "Fd-3m", "I4_1/amd", "Fddd", "C2/c"])


# ================================================================
# Tests
# ================================================================

@testset "Space groups" begin

    @testset "SymmetryOperation" begin
        sg = spacegroup(triangular)
        ops = operations(sg)
        # composition and application are consistent
        for op1 in ops, op2 in ops
            x = [0.1, 0.27]
            @test (op1 * op2)(x) ≈ op1(op2(x))
        end
        # rotations are orthogonal in Cartesian coordinates, and C6 is a 60 degree rotation
        @test all(op -> cartesian_rotation(op, triangular)' * cartesian_rotation(op, triangular) ≈ I, ops)
        angles = [atan(R[2, 1], R[1, 1]) for R in (cartesian_rotation(op, triangular) for op in ops) if det(R) > 0]
        @test any(θ -> isapprox(θ, π / 3; atol=1e-10), angles)
    end

    reference = Dict(name => bf_operations(lattice) for (name, lattice) in LATTICE_CASES)

    @testset "infinite lattices" begin
        for (name, lattice, symbol, number, schoenflies, npg, symmorphic, ntranslations) in LATTICE_CASES
            @testset "$name" begin
                sg = spacegroup(lattice)
                @test sg.symbol == symbol
                @test sg.number == number
                @test sg.schoenflies == schoenflies
                @test length(pointgroup_operations(sg)) == npg
                @test issymmorphic(sg) == symmorphic
                @test sg.ntranslations == ntranslations
                @test length(sg) == npg * ntranslations
                @test operations(sg)[1].W == I && iszero(operations(sg)[1].w)
                # independent brute-force search finds exactly the same operations
                @test same_operations(sg, reference[name])
                if symmorphic
                    # the pure rotations about the origin are symmetries
                    @test all(W -> any(op -> op[1] == W && bf_isbravais(lattice, (I - W) * sg.origin - op[2]), reference[name]),
                              pointgroup_operations(sg))
                else
                    @test isnothing(sg.origin)
                end
            end
        end

        # origins at the expected high-symmetry points
        @test spacegroup(square).origin ≈ [0.0, 0.0]
        @test spacegroup(triangular).origin ≈ [0.0, 0.0]
        @test spacegroup(honeycomb).origin ≈ [2/3, 2/3]        # center of a hexagon, not a site
        @test spacegroup(kagome).origin ≈ [1/2, 1/2]           # center of a hexagon
        @test maple_leaf.A' * spacegroup(maple_leaf).origin ≈ [2.456769074559977, 0.9819805060619656]  # legacy "Symmetry center"

        # atom types are taken into account: two inequivalent sublattices break the sixfold axis
        honeycomb_ab = Lattice(honeycomb.A, honeycomb.positions; types=[1, 2])
        @test spacegroup(honeycomb_ab).symbol == "p3m1"
    end

    @testset "finite lattices" begin
        for (name, lattice) in ((c[1], c[2]) for c in LATTICE_CASES)
            @testset "$name" begin
                cases = CLUSTER_CASES[name]
                # the clusters of each lattice realize different symmetry groups
                @test length(cases) >= 2
                @test length(unique(c[2] for c in cases)) == length(cases)

                for (boundary, symbol, npg, order) in cases
                    fl = FiniteLattice(lattice, boundary, true)
                    g = spacegroup(fl)
                    @test g.symbol == symbol
                    @test length(pointgroup_operations(g)) == npg
                    @test length(g) == order
                    @test issymmorphic(g) == !(symbol in NONSYMMORPHIC_SYMBOLS)
                    @test g.spacegroup.symbol == spacegroup(lattice).symbol

                    # brute force: rotations preserving the torus, times all translations
                    ref = reference[name]
                    @test Set(pointgroup_operations(g)) == Set(op[1] for op in ref if bf_torus_compatible(fl, op[1]))
                    @test order == count(op -> bf_torus_compatible(fl, op[1]), ref) * length(bravais_cells(fl))

                    # site permutations
                    @test length(site_permutations(g)) == length(g)
                    @test all(p -> length(p) == length(atoms(fl)), site_permutations(g))
                    @test permutations_are_correct(g)
                    @test permutations_form_group(g)
                end
            end
        end

        # only fully periodic clusters are supported
        @test_throws ArgumentError spacegroup(FiniteLattice(square, [4 0; 0 4], [true, false]))
    end

    @testset "legacy maple-leaf clusters" begin
        files = ["maple.leaf.JhexagonJtriangleJdimer.12.v1.2sl.toml",
                 "maple.leaf.JhexagonJtriangleJdimer.18.v1.2sl.toml",
                 "maple.leaf.JhexagonJtriangleJdimer.24.v2.2sl.toml",
                 "maple.leaf.JhexagonJtriangleJdimer.36.v3.2sl.toml",
                 "maple.leaf.JhexagonJtriangleJdimer.54.v1.2sl.toml"]
        parse_vector(text, name) = parse.(Float64, split(match(Regex(name * raw"=\(([^)]*)\)"), text)[1], ","))
        for file in files
            @testset "$file" begin
                path = joinpath(@__DIR__, "data", "legacy", file)
                text = read(path, String)
                data = TOML.parsefile(path)

                # torus from the header (the "Simulation torus matrix" line of legacy files is unreliable)
                t = hcat(parse_vector(text, "t1"), parse_vector(text, "t2"))
                boundary = round.(Int, (maple_leaf.A' \ t)')
                fl = FiniteLattice(maple_leaf, boundary, true)
                g = spacegroup(fl)

                @test g.spacegroup.schoenflies == match(r"Space Group \(infinite Lattice\): (\w+)", text)[1]
                @test g.schoenflies == match(r"Lattice Point Group: (\w+)", text)[1]
                @test length(g) == length(data["Symmetries"])

                # map legacy site indices onto ours by their positions modulo the torus
                X, _ = site_data(fl)
                σ = [findfirst(y -> same_site(fl, maple_leaf.A' \ Float64.(c), y), X) for c in data["Coordinates"]]
                @test !any(isnothing, σ) && isperm(σ)
                legacy = Set{Vector{Int}}()
                for q in data["Symmetries"]   # 0-based: legacy site i is mapped onto legacy site q[i]
                    p = zeros(Int, length(q))
                    for i in eachindex(q)
                        p[σ[i]] = σ[q[i] + 1]
                    end
                    push!(legacy, p)
                end
                @test legacy == Set(site_permutations(g))
            end
        end
    end
end
