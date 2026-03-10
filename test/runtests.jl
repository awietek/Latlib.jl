using Revise
using Test
using Latlib


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
    

end
