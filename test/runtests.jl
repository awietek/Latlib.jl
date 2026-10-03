using Test
using Latlib
using LinearAlgebra


@testset "Latlib.jl" begin
    
    @testset "lattice.jl" begin
        include("test_lattice.jl")
    end

    @testset "flattice.jl" begin
        include("test_flattice.jl")
    end

    @testset "metric.jl" begin
        include("test_metric.jl")
    end

    @testset "opsum.jl" begin
        include("test_opsum.jl")
    end

    @testset "predefined_lattices.jl" begin
        include("test_predefined_lattices.jl")
    end

    @testset "spacegroup.jl" begin
        include("test_spacegroup.jl")
    end

    @testset "irreps.jl" begin
        include("test_irreps.jl")
    end


end
