using LinearAlgebra

# This file tests the predefined lattice constants in src/lattice/predefined_lattices.jl.


@testset "Predefined Lattices" begin

    # ================================================================
    # Square lattice
    # ================================================================
    @testset "square" begin
        @test isa(square, Lattice)
        @test dim(square) == 2
        @test natoms(square) == 1

        # lattice vectors: identity matrix rows
        @test isapprox(square.A[1, :], [1.0, 0.0])
        @test isapprox(square.A[2, :], [0.0, 1.0])

        # single atom at the origin (in lattice coordinates)
        @test isapprox(square.positions[1, :], [0.0, 0.0])

        # nearest-neighbor distance = 1 (lattice constant)
        fl = FiniteLattice(square, [4 0; 0 4], true)
        coords = atoms(fl)
        dists = distances(coords; flattice=fl)
        @test isapprox(dists[2], 1.0)   # dists[1] == 0 (self-distance)

        # next-nearest-neighbor distance = sqrt(2) (diagonal)
        @test isapprox(dists[3], sqrt(2.0))
    end

    # ================================================================
    # Triangular lattice
    # ================================================================
    @testset "triangular" begin

        @test isa(triangular, Lattice)
        @test dim(triangular) == 2
        @test natoms(triangular) == 1

        # all lattice vectors have unit length
        @test isapprox(sum(triangular.A[1, :].^2), 1.0)
        @test isapprox(sum(triangular.A[2, :].^2), 1.0)

        # the angle between a1 and a2 is 120° (cos = –0.5)
        a1 = triangular.A[1, :]
        a2 = triangular.A[2, :]
        cos_theta = dot(a1, a2) / (LinearAlgebra.norm(a1) * LinearAlgebra.norm(a2))
        @test isapprox(cos_theta, -0.5; atol=1e-10)

        # single atom at origin
        @test isapprox(triangular.positions[1, :], [0.0, 0.0])

        # nearest-neighbor distance = 1 on a 4×4 periodic patch
        fl = FiniteLattice(triangular, [4 0; 0 4], true)
        coords = atoms(fl)
        dists = distances(coords; flattice=fl)
        @test isapprox(dists[2], 1.0)

        # each site has exactly 6 nearest neighbours
        nbs = neighbors(coords; num_distance=1, flattice=fl)
        nsites = length(coords)
        @test length(nbs) == 3 * nsites   # 6 per site / 2 (each bond counted once)
    end

    # ================================================================
    # Kagome lattice
    # ================================================================
    @testset "kagome" begin
        @test isa(kagome, Lattice)
        @test dim(kagome) == 2
        @test natoms(kagome) == 3

        # lattice vectors
        @test isapprox(kagome.A[1, :], [1.0, 0.0])
        @test isapprox(kagome.A[2, :], [0.5, sqrt(3)/2])

        # three distinct atom positions (in lattice coordinates)
        @test isapprox(kagome.positions[1, :], [0.0, 0.0])
        @test isapprox(kagome.positions[2, :], [0.5, 0.0])
        @test isapprox(kagome.positions[3, :], [0.0, 0.5])

        # nearest-neighbor distance = 0.5 (half the lattice vector length)
        fl = FiniteLattice(kagome, [4 0; 0 4], true)
        coords = atoms(fl)
        dists = distances(coords; flattice=fl)
        @test isapprox(dists[2], 0.5)

        # each site has exactly 4 nearest neighbours on the kagome lattice
        nbs = neighbors(coords; num_distance=1, flattice=fl)
        nsites = length(coords)
        @test length(nbs) == 2 * nsites   # 4 per site / 2
    end



    @testset "Hyperhoneycomb" begin

        @testset "Lattice Construction" begin
            @test isa(hyperhoneycomb, Lattice)
            @test hyperhoneycomb.dim == 3
            @test length(positions(hyperhoneycomb)) == 4
            @test hyperhoneycomb.A == [2.0 4.0 0.0; 3.0 3.0 2.0; -1.0 1.0 2.0]
        end

        # create N = 16 finite lattice with periodic boundaries (Alex's TOML file as reference)
        flat_16_vecs = [
            LatticeVector(hyperhoneycomb, [-1, 1, 1]),  # t1
            LatticeVector(hyperhoneycomb, [1, 1, -1]),  # t2
            LatticeVector(hyperhoneycomb, [-1, 1, -1]), # t3
        ]

        @test to_euclidean_basis(flat_16_vecs[1]).coords == [0.0, 0.0, 4.0]   # t1 euclidean
        @test to_euclidean_basis(flat_16_vecs[2]).coords == [6.0, 6.0, 0.0]   # t2 euclidean
        @test to_euclidean_basis(flat_16_vecs[3]).coords == [2.0, -2.0, 0.0]  # t3 euclidean

        flat_16 = FiniteLattice(flat_16_vecs, true)
        @test dim(flat_16) == 3
        @test natoms(flat_16) == 4
        @test periodicity(flat_16) == [true, true, true]
        @test length(bravais_cells(flat_16)) == 4
        @test length(atoms(flat_16)) == 16

        println("coordinates of 16-site hyperhoneycomb cluster")
        for (i, v_euc) in enumerate(atoms(flat_16))
            println("$i: ", v_euc)
        end

        # create N = 32 finite lattice with periodic boundaries
        flat_32_vecs = [
            LatticeVector(hyperhoneycomb, [-1, 1, 1]),
            LatticeVector(hyperhoneycomb, [1, 1, -1]),
            LatticeVector(hyperhoneycomb, [-2, 1, -2]),
        ]
        println("Simulation torus vectors for 32-site hyperhoneycomb cluster:")
        for (i, v) in enumerate(flat_32_vecs)
            @show v
        end


        flat_32 = FiniteLattice(flat_32_vecs, true)
        @test dim(flat_32) == 3
        @test natoms(flat_32) == 4
        @test periodicity(flat_32) == [true, true, true]
        @test length(bravais_cells(flat_32)) == 8
        @test length(atoms(flat_32)) == 32

        println("coordinates of 32-site hyperhoneycomb cluster")
        for (i, v_euc) in enumerate(atoms(flat_32))
            println("$i: ", v_euc)
        end








    end

end
