# OpSum

```@docs
Op
OpSum
neighbor_bonds(type::AbstractString, coupling::AbstractString, lattice::FiniteLattice; num_distance::Integer=1)
lattice_bonds(type::AbstractString, coupling::AbstractString, Flattice::FiniteLattice, b1::Vector{<:Integer}, b2::Vector{<:Integer})
write_opsum_to_toml!(opsum::OpSum, filename::String; index_zero::Bool=false)
```