# FiniteLattice

```@docs
FiniteLattice
FiniteLattice(lattice::Lattice, boundary::AbstractMatrix; periodic=true)
coordinates(flattice::FiniteLattice)
bravais_coordinates(flattice::FiniteLattice)
boundary_vectors(flattice::FiniteLattice)
periodicity_vectors(flattice::FiniteLattice)
neighbors(flattice::FiniteLattice)
distance(x1::AbstractVector, x2::AbstractVector, flattice::FiniteLattice)
distances(flattice::FiniteLattice)
dimension(flattice::FiniteLattice)
natoms(flattice::FiniteLattice)
```