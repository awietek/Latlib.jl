# ----------------------------------------------------------------------
#        Dimensions of the sectors in the spin-1/2 Hilbert space
# ----------------------------------------------------------------------
#
# A spin configuration is invariant under a site permutation h if and only if the spins are
# equal along every cycle of h, so the trace of h on the spin-1/2 Hilbert space is
# 2^(number of cycles of h). The dimension of a sector with characters χ on a group H is
#
#     dim = (1/|H|) Σ_{h ∈ H} conj(χ(h)) 2^(cycles(h)).
#
# This is evaluated exactly. The elements of H with a given cycle type are closed under
# h ↦ h^m for m coprime to |H|, which maps χ(h) to its Galois conjugate χ(h)^m. The sum of
# conj(χ) over these elements is therefore a rational algebraic integer, i.e. an integer.
# Characters that do not form a one-dimensional representation generally violate this.

# cycle type (sorted cycle lengths) of a permutation
function _cycle_type(p::Vector{Int}) :: Vector{Int}
    seen = falses(length(p))
    lengths = Int[]
    for i in eachindex(p)
        seen[i] && continue
        len = 0
        j = i
        while !seen[j]
            seen[j] = true
            j = p[j]
            len += 1
        end
        push!(lengths, len)
    end
    return sort!(lengths)
end

# cycle type of each permutation, as an index into the list of distinct cycle types
function _cycle_types(permutations::Vector{Vector{Int}})
    index = Dict{Vector{Int}, Int}()
    types = Vector{Int}[]
    ids = Int[]
    for p in permutations
        type = _cycle_type(p)
        if !haskey(index, type)
            push!(types, type)
            index[type] = length(types)
        end
        push!(ids, index[type])
    end
    return ids, types
end

# number of spin configurations (with `nup` up spins, or all) invariant under a permutation of
# cycle type `type`: coefficient of x^nup in Π_cycles (1 + x^length)
function _fixed_configurations(type::Vector{Int}, nup::Union{Nothing, Int}) :: BigInt
    isnothing(nup) && return big(2)^length(type)
    c = zeros(BigInt, nup + 1)
    c[1] = 1
    for len in type, n in nup+1:-1:len+1
        c[n] += c[n - len]
    end
    return c[nup + 1]
end

# number of spin configurations invariant under a permutation of cycle type `type` followed by
# the global spin flip: the spins alternate along every cycle (all with N/2 up spins)
_fixed_configurations_flipped(type::Vector{Int}) :: BigInt = all(iseven, type) ? big(2)^length(type) : big(0)

# Dimension of the sector with `characters` on the permutations with cycle types `ids` (into
# `types`), which form a group. Throws if the characters fail the integrality conditions.
function _sector_dimension(ids::Vector{Int}, types::Vector{Vector{Int}}, characters::Vector{ComplexF64},
                           label::String; nup::Union{Nothing, Int}=nothing, spinflip::Union{Nothing, Int}=nothing) :: BigInt
    sums = Dict{Int, ComplexF64}()
    for (id, χ) in zip(ids, characters)
        sums[id] = get(sums, id, 0.0im) + conj(χ)
    end
    total = big(0)
    for (id, s) in sums
        n = round(Int, real(s))
        if abs(real(s) - n) > 1e-6 || abs(imag(s)) > 1e-6
            error("Consistency check failed for $label: the sum of the characters over the operations with cycle type $(types[id]) is $s, not an integer, so the characters are not a one-dimensional representation. This is a bug in Latlib.")
        end
        trace = _fixed_configurations(types[id], nup)
        isnothing(spinflip) || (trace += spinflip * _fixed_configurations_flipped(types[id]))
        total += n * trace
    end
    order = length(ids) * (isnothing(spinflip) ? 1 : 2)
    dimension, remainder = divrem(total, order)
    if remainder != 0 || dimension < 0
        error("Consistency check failed for $label: the dimension $total/$order is not a non-negative integer. This is a bug in Latlib.")
    end
    return dimension
end

@doc raw"""
    sector_dimension(cs::ClusterSymmetries, irrep::Irrep; nup=nothing, spinflip=nothing) -> BigInt

Dimension of the sector `irrep` (see [`irreps`](@ref)) in the Hilbert space of spin-1/2 on the
sites of the finite lattice of `cs`, i.e. the size of the block an exact diagonalization code
works with.

It is computed exactly from the cycles of the site permutations and does not enumerate any states:
a spin configuration is invariant under a permutation ``h`` if and only if the spins are equal along
every cycle of ``h``, so the dimension is
```math
\dim = \frac{1}{|H|} \sum_{h \in H} \overline{\chi(h)}\, 2^{\,\mathrm{cycles}(h)},
```
where ``H`` are the `allowed_symmetries` of the sector and ``\chi`` its characters.

# Keyword arguments
- `nup`: number of up spins, i.e. the sector ``S^z = n_\uparrow - N/2``. If `nothing`, all
  configurations are counted.
- `spinflip`: `+1` or `-1` restricts the sector to states that are even or odd under the global
  spin flip. Requires `nup = N/2` or `nup = nothing` (the spin flip maps ``n_\uparrow`` to
  ``N - n_\uparrow``).

The dimensions of all sectors, weighted with the sizes of the stars of their momenta, add up to
``2^N`` (or ``\binom{N}{n_\uparrow}``); this is checked whenever the irreps are computed.

# Examples
```julia
cs = symmetries(FiniteLattice(square, [6 0; 0 6], true); origin=LatticeVector(square, [0.0, 0.0]))
irrep = only(i for i in irreps(cs) if i.label == "Gamma.C4v.A1")
sector_dimension(cs, irrep; nup=18, spinflip=1)   # 15804956
```
"""
function sector_dimension(cs::ClusterSymmetries, irrep::Irrep; nup::Union{Nothing, Integer}=nothing,
                          spinflip::Union{Nothing, Integer}=nothing) :: BigInt
    N = length(cs.permutations[1])
    if !isnothing(nup) && !(0 <= nup <= N)
        throw(ArgumentError("nup = $nup must be between 0 and the number of sites $N."))
    end
    if !isnothing(spinflip)
        spinflip in (1, -1) || throw(ArgumentError("spinflip must be +1 or -1, got $spinflip."))
        if !isnothing(nup) && 2 * nup != N
            throw(ArgumentError("The spin flip maps nup = $nup to $(N - nup) up spins; it is a symmetry only for nup = N/2."))
        end
    end
    if isempty(irrep.allowed_symmetries)
        throw(ArgumentError("$(irrep.label) is a skipped three-dimensional irrep without characters."))
    end
    ids, types = _cycle_types(cs.permutations[irrep.allowed_symmetries])
    return _sector_dimension(ids, types, irrep.characters, irrep.label;
                             nup=isnothing(nup) ? nothing : Int(nup), spinflip=isnothing(spinflip) ? nothing : Int(spinflip))
end

# Checks the sectors of the representatives of the stars `ks` (see `_sectors`) against the
# decomposition of the spin-1/2 Hilbert space by momentum, which follows from the translations
# alone. At every representative momentum k, the dimensions of the sectors must add up to the
# dimension D_k of the subspace with momentum k, apart from skipped three-dimensional irreps
# (which contribute 3 × their multiplicity), and Σ_k |star(k)| D_k = 2^N. The dimension of
# every sector is computed exactly, which also checks the characters (see above).
function _check_sectors(cs::ClusterSymmetries, ks::Vector{ClusterMomentum}, sectors::Vector{Tuple{Irrep, Symbol}})
    ops = operations(cs.spacegroup)
    N = length(cs.permutations[1])
    ids, types = _cycle_types(cs.permutations)
    translations = [j for j in eachindex(cs.operations) if ops[cs.operations[j]].W == I]
    tids = ids[cs.permutation_index[translations]]
    starsize = Dict{Int, Int}()
    for k in ks
        starsize[k.star] = get(starsize, k.star, 0) + 1
    end
    total = big(0)
    for k in ks
        k.representative || continue
        χ = [cispi(2 * Float64(mod(sum(k.coords .* cs.translations[j]), 1))) for j in translations]
        Dk = _sector_dimension(tids, types, χ, "momentum $(k.label)")
        dims = big(0)
        skipped = false
        for (irrep, status) in sectors
            irrep.kpoint == k.label || continue
            skipped |= status == :skipped
            status == :ok || continue
            dims += _sector_dimension(ids[irrep.allowed_symmetries], types, irrep.characters, irrep.label)
        end
        rest = Dk - dims
        if rest < 0 || (skipped ? rem(rest, 3) != 0 : rest != 0)
            error("Consistency check failed at momentum $(k.label): the sectors have $dims states, the subspace of this momentum has $Dk" *
                  (skipped ? " (the difference must be a non-negative multiple of 3 because of skipped three-dimensional irreps)." : ".") *
                  " This is a bug in Latlib.")
        end
        total += starsize[k.star] * Dk
    end
    total == big(2)^N || error("Consistency check failed: the momenta of the stars span $total states instead of 2^$N. This is a bug in Latlib.")
    return nothing
end
