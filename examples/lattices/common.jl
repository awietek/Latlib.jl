# Helper functions shared by the example scripts in this directory (included by every script).

using Latlib
using LinearAlgebra

# All nonzero torus vectors (Cartesian) of a fully periodic finite lattice with length <= ρ. The
# search box is provably large enough: every torus vector T c with |T c| <= ρ satisfies
# |c_i| <= ρ |row i of pinv(T)|.
function torus_vectors(fl::FiniteLattice, ρ::Real) :: Vector{Vector{Float64}}
    T = Matrix{Float64}((fl.boundary * fl.lattice.A)')   # columns: boundary vectors (Cartesian)
    P = pinv(T)
    ranges = [-ceil(Int, ρ * norm(P[i, :])):ceil(Int, ρ * norm(P[i, :])) for i in axes(T, 2)]
    candidates = (T * collect(c) for c in Iterators.product(ranges...) if any(!=(0), c))
    return [t for t in candidates if norm(t) <= ρ + 1e-8]
end

# Checks that the cluster represents the couplings of H properly. `couplings` maps every coupling
# name to (bonds per unit cell, bond length). For every coupling,
#  - the number of bonds is the number per unit cell times the number of unit cells,
#  - every bond has the expected length under the periodic metric,
#  - no bond connects a site to itself and no pair of sites is coupled twice,
#  - every bond has a unique shortest periodic image, i.e. it does not wrap around the torus.
# On a cluster that is too small for a coupling, one of these fails.
function check_couplings(fl::FiniteLattice, H::OpSum, couplings::Dict)
    ncells = length(bravais_cells(fl))
    xs = atoms(fl)
    for (name, (bonds, len)) in couplings
        ops = [op for op in H.ops if op.cpl == name]
        if length(ops) != bonds * ncells
            error("Coupling $name has $(length(ops)) bonds, expected $bonds per unit cell × $ncells unit cells.")
        end
        ts = torus_vectors(fl, 2 * len)   # a second image as short as r requires |t| <= 2|r|
        pairs = Set{Tuple{Int, Int}}()
        for op in ops
            i, j = op.sites
            i == j && error("Coupling $name connects site $i to itself.")
            pair = (min(i, j), max(i, j))
            pair in pairs && error("Coupling $name connects sites $i and $j twice.")
            push!(pairs, pair)
            r = distance_vector(xs[i], xs[j]; flattice=fl).coords
            isapprox(norm(r), len; atol=1e-6) || error("Bond ($i, $j) of coupling $name has length $(norm(r)), expected $len.")
            if any(t -> norm(r - t) <= norm(r) + 1e-8, ts)
                error("Bond ($i, $j) of coupling $name wraps around the torus: it has two shortest periodic images.")
            end
        end
    end
    unchecked = setdiff(Set(op.cpl for op in H.ops), keys(couplings))
    isempty(unchecked) || error("Couplings without expected bonds: $(collect(unchecked)).")
    return nothing
end

# File name `<lattice>-<model>-N-<sites>-ver-<version>.toml` in the directory `dir`
toml_name(dir::String, lattice::String, model::String, fl::FiniteLattice, ver::Int) =
    joinpath(dir, "$lattice-$model-N-$(length(atoms(fl)))-ver-$ver.toml")

# Shows a figure and waits until its window is closed
function show_figure(f)
    println("    (close the plot window to continue)")
    wait(display(f))
end
