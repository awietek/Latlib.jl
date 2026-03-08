using Revise
using Test
using Latlib


@testset "Latlib.jl" begin
    
    @testset "Lattice class" begin
        include("test_lattice.jl")
    end

    @testset "FiniteLattice class" begin
        include("test_flattice.jl")
    end

end
