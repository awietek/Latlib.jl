# Generates src/symmetry/kpoint_tables.jl: the k-vector types (CDML labels, as used by ISOTROPY and
# the Bilbao Crystallographic Server) of the 14 holohedries, extracted from Crystalline.jl.
# Run from the repository root with
#
#     julia dev/generate_kpoint_tables.jl
#
# Crystalline.jl is only needed here, not as a dependency of Latlib.

using Pkg
Pkg.activate(; temp=true)
Pkg.add("Crystalline")

using Crystalline, LinearAlgebra
const HOLOHEDRIES = [(2, "P-1"), (10, "P2/m"), (12, "C2/m"), (47, "Pmmm"), (65, "Cmmm"), (69, "Fmmm"), (71, "Immm"),
                     (123, "P4/mmm"), (139, "I4/mmm"), (166, "R-3m"), (191, "P6/mmm"), (221, "Pm-3m"), (225, "Fm-3m"), (229, "Im-3m")]
greek = Dict('Γ' => "Gamma", 'Δ' => "Delta", 'Λ' => "Lambda", 'Σ' => "Sigma", 'Ω' => "Omega", 'Π' => "Pi", 'Φ' => "Phi", 'Ψ' => "Psi", 'Ξ' => "Xi", 'Θ' => "Theta")
ascii(l) = join([haskey(greek, c) ? greek[c] : string(c) for c in l])
fmt(x) = (r = rationalize(x; tol=1e-9); denominator(r) == 1 ? string(numerator(r)) * ".0" : "$(numerator(r))/$(denominator(r))")
vec_str(v) = "[" * join(fmt.(v), ", ") * "]"
allnames = Set{String}()
out = IOBuffer()
println(out, "# k-vector types of the 14 holohedries (symmorphic space groups of the Bravais lattices),")
println(out, "# in the CDML notation of ISOTROPY and of the Bilbao Crystallographic Server (REPRES), taken from")
println(out, "# Crystalline.jl (lgirreps) with the script dev/generate_kpoint_tables.jl. Coordinates refer to the conventional")
println(out, "# reciprocal basis of the ITA standard setting: k = k0 + Σ_i t_i d_i.")
println(out, "#   (label, k0, directions d_i, order of the little co-group)")
println(out, "const _KPOINT_LABELS_3D = Dict{Int, Vector{Tuple{String, Vector{Float64}, Vector{Vector{Float64}}, Int}}}(")
for (sg, sym) in HOLOHEDRIES
    lgirs = lgirreps(sg, Val(3))
    entries = String[]
    rows = []
    for (label, irs) in lgirs
        label == "Ω" && continue
        kv = position(first(irs))
        dirs = [Vector{Float64}(kv.free[:, j]) for j in 1:3 if any(!iszero, kv.free[:, j])]
        push!(rows, (ascii(label), Vector{Float64}(kv.cnst), dirs, length(group(first(irs)))))
        push!(allnames, label)
    end
    sort!(rows; by=r -> (length(r[3]), r[1]))
    println(out, "    $sg => [   # $sym")
    for (l, c, d, n) in rows
        println(out, "        (\"$l\", $(vec_str(c)), [" * join(vec_str.(d), ", ") * "], $n),")
    end
    println(out, "    ],")
end
println(out, ")")
write(joinpath(@__DIR__, "..", "src", "symmetry", "kpoint_tables.jl"), take!(out))
println("labels used: ", sort(collect(allnames)))
