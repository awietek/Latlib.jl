# Example clusters of the predefined lattices

One directory per predefined lattice of Latlib. Each script builds a few periodic clusters, defines
a Hamiltonian, checks that every cluster is large enough for its couplings, and writes one TOML
file per cluster and model (`<lattice>-<model>-N-<sites>-ver-<version>.toml`, at most 1 MB each)
into its directory. The TOML files are not part of the repository (see `.gitignore`); run the
scripts to generate them.
Run the scripts from the root of the repository, e.g.

```
julia --project=. examples/lattices/square/square.jl
```

Every cluster is plotted as well, unless `show_plot = false` is set at the top of the script.

| Lattice | Model | Sites | Symmetries in the TOML files |
|:--------|:------|:------|:-----------------------------|
| `square` | J1-J2 Heisenberg | 16, 18, 20, 32, 36 | operations and irreps |
| `triangular` | J1-J2 Heisenberg | 16, 21, 27, 36 | operations and irreps |
| `honeycomb` | J1-J2 Heisenberg; Kitaev-Heisenberg | 18, 24, 32 | operations and irreps¹ |
| `kagome` | J1-J2 Heisenberg | 12, 27, 36 | operations and irreps |
| `shastry_sutherland` | Shastry-Sutherland (J, Jd) | 16, 20, 32, 36 | operations and irreps of the sites² |
| `shastry_sutherland_non_symmorphic` | Shastry-Sutherland (J, Jd) | 16, 20, 32, 36 | operations |
| `lieb` | J1-J2 Heisenberg | 12, 24, 27, 48 | operations and irreps |
| `trellis` | trellis model (Jleg, Jrung, Jinter) | 16, 24, 32, 48 | operations and irreps |
| `maple_leaf` | Jhexagon, Jtriangle, Jdimer | 12, 18, 24, 36, 42, 54 | operations and irreps |
| `hyperhoneycomb` | Heisenberg; Kitaev-Heisenberg | 8, 16, 32, 64 | operations¹ |
| `simple_cubic` | J1-J2 Heisenberg | 27 | operations and irreps |
| `bcc` | J1-J2 Heisenberg | 27 | operations and irreps |
| `fcc` | J1-J2 Heisenberg | 27 | operations and irreps |
| `diamond` | J1-J2 Heisenberg | 32, 64 | operations |
| `pyrochlore` | J1-J2 Heisenberg | 32, 64 | operations |
| `simple_hexagonal` | J1-J2 Heisenberg | 27 | operations and irreps |
| `hcp` | J1-J2 Heisenberg | 36, 48, 54 | operations |

Irreducible representations are written for the symmorphic space groups only; three-dimensional
irreps of cubic little co-groups are skipped (with a banner in the file).

The symmetries are those of the lattice, which can be larger than those of the Hamiltonian:

1. The Kitaev couplings are invariant under translations, but the point-group operations exchange
   the bonds KX, KY and KZ; they are symmetries only together with a transformation of the spin
   components.
2. The sites of `shastry_sutherland` form a square lattice (p4mm); only part of its operations are
   symmetries of the Shastry-Sutherland model (p4gm). `shastry_sutherland_non_symmorphic` has the
   symmetry of the model built into its geometry.

The checks (`check_couplings` in `common.jl`) require, for every coupling, the expected number of
bonds per unit cell, all of the expected length under the periodic metric, no self-interactions,
no pair of sites coupled twice, and a unique shortest periodic image of every bond, so that no bond
wraps around the torus.
