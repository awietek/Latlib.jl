# Checks TOML files with symmetries written by Latlib with XDiag (https://github.com/awietek/xdiag).
# Not part of the test suite. Run from the repository root with
#
#     julia test/xdiag/check_with_xdiag.jl
#
# It creates a temporary environment with this version of Latlib and XDiag, writes TOML files
# for a few clusters, and checks for each of them that
#  - XDiag reads every irreducible representation (XDiag verifies that the characters form a
#    one-dimensional representation of the group of allowed symmetries),
#  - the dimensions of all sectors, weighted with the sizes of the stars, add up to 2^N,
# and, for the legacy maple-leaf cluster, that the ground-state energies of all sectors agree
# with those obtained from the legacy file.

using Pkg
Pkg.activate(; temp=true)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..", "..")))
Pkg.add("XDiag")

using Latlib, XDiag, TOML

function sector_labels(file)
    data = TOML.parsefile(file)
    return sort(["$k.$g.$i" for (k, v) in data if v isa Dict for (g, v2) in v for (i, _) in v2])
end

function check(name, fl, cs)
    file = joinpath(mktempdir(), "$name.toml")
    write_toml(fl, neighbor_interaction("SdotS", "J", fl), file; zero_based=true, symmetries=cs)
    ks = momenta(cs)
    starsize = Dict(k.label => count(q -> q.star == k.star, ks) for k in ks if k.representative)
    N = length(atoms(fl))
    f = FileToml(file)
    total = 0
    for label in sector_labels(file)
        total += starsize[split(label, ".")[1]] * size(Spinhalf(N, read_representation(f, label)))
    end
    println(rpad(name, 22), length(sector_labels(file)), " sectors read by XDiag, Σ|star|·dim = ", total,
            total == 2^N ? " = 2^$N" : " ≠ 2^$N  FAILED")
    return file
end

site(lattice) = LatticeVector(lattice, [0.0, 0.0])
clusters = [
    ("maple_leaf_18", FiniteLattice(maple_leaf, [1 1; 1 -2], true), nothing),
    ("triangular_9", FiniteLattice(triangular, [3 0; 0 3], true), nothing),
    ("honeycomb_8", FiniteLattice(honeycomb, [2 -2; 1 1], true), nothing),
    ("kagome_12", FiniteLattice(kagome, [2 0; 0 2], true), nothing),
    ("square_16", FiniteLattice(square, [4 0; 0 4], true), site(square)),
    ("square_8_thin", FiniteLattice(square, [4 0; 0 2], true), site(square)),
]
files = Dict(name => check(name, fl, symmetries(fl; origin=origin)) for (name, fl, origin) in clusters)

# ground-state energies of all sectors: new maple-leaf file vs legacy file (all couplings 1)
legacy = joinpath(@__DIR__, "..", "data", "legacy", "maple.leaf.JhexagonJtriangleJdimer.18.v1.2sl.toml")
new, old = FileToml(files["maple_leaf_18"]), FileToml(legacy)
Hnew = read_opsum(new, "Interactions"); Hnew["J"] = 1.0
Hold = read_opsum(old, "Interactions")
for c in ("Jhexagon", "Jtriangle", "Jdimer")
    Hold[c] = 1.0
end
for label in sector_labels(files["maple_leaf_18"])
    e_new = eigval0(Hnew, Spinhalf(18, 9, read_representation(new, label)))
    e_old = eigval0(Hold, Spinhalf(18, 9, read_representation(old, label)))
    println(rpad(label, 14), " E0 = ", round(e_new; digits=10), "   legacy: ", round(e_old; digits=10),
            abs(e_new - e_old) < 1e-8 ? "" : "  FAILED")
end
