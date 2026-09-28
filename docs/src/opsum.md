# Operators and interactions

A Hamiltonian is represented as an [`OpSum`](@ref), a sum of operators
[`Op`](@ref). Each `Op` carries a type string (e.g. `"SdotS"`), a coupling
(e.g. `"J"`), and the 1-based indices of the sites it acts on.

There are two ways to generate interactions on a [`FiniteLattice`](@ref):

- [`neighbor_interaction`](@ref) couples all k-th nearest neighbors, using the
  periodic metric of the finite lattice.
- [`lattice_interaction`](@ref) couples two specific atoms of the unit cell,
  possibly in different Bravais cells, and repeats the bond over all cells of
  the finite lattice. This is needed for bond-dependent models such as the
  Kitaev model.

```@example opsum
using Latlib

fl = FiniteLattice(honeycomb, [2 -2; 1 1], true)

# Heisenberg-Kitaev model
H = OpSum()
H += neighbor_interaction("SdotS", "J", fl; num_distance=1)
H += lattice_interaction("SxSx", "KX", fl, 1, 2, [0, 0])   # atom 1 to atom 2 in the same cell
H += lattice_interaction("SySy", "KY", fl, 1, 2, [0, -1])  # atom 1 to atom 2 in cell [0, -1]
H += lattice_interaction("SzSz", "KZ", fl, 2, 1, [1, 0])   # atom 2 to atom 1 in cell [1, 0]
H.ops
```

```@docs
Op
OpSum
neighbor_interaction
lattice_interaction
unique_ops!
```
