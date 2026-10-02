# Distances and neighbors

Distances between sites can be computed either in plain Euclidean space or
modulo the periodic boundary vectors of a [`FiniteLattice`](@ref). All
functions below take an optional keyword argument `flattice`; if given, the
periodic metric of that finite lattice is used.

```@example metric
using Latlib

fl = FiniteLattice(square, [4 0; 0 4], true)
x = EuclideanVector([0, 0])
y = EuclideanVector([3, 0])

(distance(x, y), distance(x, y; flattice=fl))
```

The nearest neighbors of a finite lattice are pairs of site indices:

```@example metric
neighbors(fl)
```

```@docs
distance
distance_vector
distance_matrix
distances
neighbors
EuclideanMetric
PeriodicEuclideanMetric
```
