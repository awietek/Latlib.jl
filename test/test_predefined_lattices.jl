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


        @testset "N=16 cluster" begin
        
            # create N = 16 finite lattice with periodic boundaries (Alex's TOML file as reference)
            flat_16_vecs = [
                LatticeVector(hyperhoneycomb, [-1, 1, 1]),  # t1
                LatticeVector(hyperhoneycomb, [1, 1, -1]),  # t2
                LatticeVector(hyperhoneycomb, [-1, 1, -1]), # t3
            ]

            flat_16 = FiniteLattice(flat_16_vecs, true)
            @test dim(flat_16) == 3
            @test natoms(flat_16) == 4
            @test periodicity(flat_16) == [true, true, true]
            @test length(bravais_cells(flat_16)) == 4
            @test length(atoms(flat_16)) == 16

            # check TOML coordinates string
            coord_string = toml_coordinates(flat_16)
            coord_string_ref = "Coordinates = [
                [1.0, -1.0, 2.0],
                [0.0, 0.0, 0.0],
                [4.0, 2.0, 0.0],
                [3.0, 3.0, 2.0],
                [2.0, -0.0, 2.0],
                [1.0, 1.0, -0.0],
                [5.0, 3.0, -0.0],
                [4.0, 4.0, 2.0],
                [2.0, 1.0, 3.0],
                [1.0, 2.0, 1.0],
                [5.0, 4.0, 1.0],
                [4.0, 5.0, 3.0],
                [3.0, 2.0, 3.0],
                [2.0, 3.0, 1.0],
                [6.0, 5.0, 1.0],
                [5.0, 6.0, 3.0],
            ]"
            @test replace(coord_string, r"\s+" => "") == replace(coord_string_ref, r"\s+" => "") # compares the two strings after removing all whitespace or line breaks

            # check TOML interaction string
            opsum_16_HB = neighbor_interaction("SdotS", "J", flat_16)
            interaction_string = toml_interactions(opsum_16_HB; zero_based=true)
            interaction_string_ref = "Interactions = [
                ['J', 'SdotS', 0, 4],
                ['J', 'SdotS', 0, 14],
                ['J', 'SdotS', 0, 15],
                ['J', 'SdotS', 1, 5],
                ['J', 'SdotS', 1, 14],
                ['J', 'SdotS', 1, 15],
                ['J', 'SdotS', 2, 6],
                ['J', 'SdotS', 2, 12],
                ['J', 'SdotS', 2, 13],
                ['J', 'SdotS', 3, 7],
                ['J', 'SdotS', 3, 12],
                ['J', 'SdotS', 3, 13],
                ['J', 'SdotS', 4, 8],
                ['J', 'SdotS', 4, 9],
                ['J', 'SdotS', 5, 8],
                ['J', 'SdotS', 5, 9],
                ['J', 'SdotS', 6, 10],
                ['J', 'SdotS', 6, 11],
                ['J', 'SdotS', 7, 10],
                ['J', 'SdotS', 7, 11],
                ['J', 'SdotS', 8, 12],
                ['J', 'SdotS', 9, 13],
                ['J', 'SdotS', 10, 14],
                ['J', 'SdotS', 11, 15],
                ]"
            @test replace(interaction_string, r"\s+" => "") == replace(interaction_string_ref, r"\s+" => "")

            # check Kitaev interaction string
            opsum_16_kitaev = lattice_interaction("SxSx", "KX", flat_16, 1, 2, [0, 0, 0])
            opsum_16_kitaev += lattice_interaction("SxSx", "KX", flat_16, 3, 4, [0, 0, 0])
            opsum_16_kitaev += lattice_interaction("SySy", "KY", flat_16, 2, 3, [0, 0, 0])
            opsum_16_kitaev += lattice_interaction("SySy", "KY", flat_16, 4, 1, [0, 1, 0])
            opsum_16_kitaev += lattice_interaction("SzSz", "KZ", flat_16, 3, 2, [0, 0, 1])
            opsum_16_kitaev += lattice_interaction("SzSz", "KZ", flat_16, 4, 1, [1, 0, 0])
            interaction_string = toml_interactions(opsum_16_kitaev; zero_based=true)
            interaction_string_ref = "Interactions = [
                ['KX', 'SxSx', 0, 4],
                ['KX', 'SxSx', 1, 5],
                ['KX', 'SxSx', 2, 6],
                ['KX', 'SxSx', 3, 7],
                ['KX', 'SxSx', 8, 12],
                ['KX', 'SxSx', 9, 13],
                ['KX', 'SxSx', 10, 14],
                ['KX', 'SxSx', 11, 15],
                ['KY', 'SySy', 4, 8],
                ['KY', 'SySy', 5, 9],
                ['KY', 'SySy', 6, 10],
                ['KY', 'SySy', 7, 11],
                ['KY', 'SySy', 2, 12],
                ['KY', 'SySy', 3, 13],
                ['KY', 'SySy', 0, 14],
                ['KY', 'SySy', 1, 15],
                ['KZ', 'SzSz', 5, 8],
                ['KZ', 'SzSz', 4, 9],
                ['KZ', 'SzSz', 7, 10],
                ['KZ', 'SzSz', 6, 11],
                ['KZ', 'SzSz', 3, 12],
                ['KZ', 'SzSz', 2, 13],
                ['KZ', 'SzSz', 1, 14],
                ['KZ', 'SzSz', 0, 15],
                ]"
            @test replace(interaction_string, r"\s+" => "") == replace(interaction_string_ref, r"\s+" => "")
            # check if the interaction string is the same as the one obtained from reading the TOML file in the example
            
            toml_opsum_16 = read_toml_interaction(joinpath(@__DIR__, "..", "examples/Hyperhoneycomb/hyperhoneycomb-N-16-ver-1.toml"); zero_based=true)            
            @test (opsum_16_HB + opsum_16_kitaev) == toml_opsum_16
        end

        @testset "N=32 cluster" begin

            # create N = 32 finite lattice with periodic boundaries
            flat_32_vecs = [
                LatticeVector(hyperhoneycomb, [-1, 1, 1]),
                LatticeVector(hyperhoneycomb, [1, 1, -1]),
                LatticeVector(hyperhoneycomb, [-2, 1, -2]),
            ]
            flat_32 = FiniteLattice(flat_32_vecs, true)
            @test dim(flat_32) == 3
            @test natoms(flat_32) == 4
            @test periodicity(flat_32) == [true, true, true]
            @test length(bravais_cells(flat_32)) == 8
            @test length(atoms(flat_32)) == 32

            # check TOML coordinates string
            coord_string = toml_coordinates(flat_32)
            coord_string_ref = "Coordinates = [
                [3.0, -3.0, 2.0],
                [2.0, -2.0, 0.0],
                [1.0, -1.0, 2.0],
                [6.0, 0.0, 0.0],
                [5.0, 1.0, 2.0],
                [0.0, 0.0, 0.0],
                [4.0, 2.0, 0.0],
                [3.0, 3.0, 2.0],
                [4.0, -2.0, 2.0],
                [3.0, -1.0, -0.0],
                [2.0, -0.0, 2.0],
                [7.0, 1.0, 0.0],
                [6.0, 2.0, 2.0],
                [1.0, 1.0, -0.0],
                [5.0, 3.0, -0.0],
                [4.0, 4.0, 2.0],
                [4.0, -1.0, 3.0],
                [3.0, 0.0, 1.0],
                [2.0, 1.0, 3.0],
                [7.0, 2.0, 1.0],
                [6.0, 3.0, 3.0],
                [1.0, 2.0, 1.0],
                [5.0, 4.0, 1.0],
                [4.0, 5.0, 3.0],
                [5.0, -0.0, 3.0],
                [4.0, 1.0, 1.0],
                [3.0, 2.0, 3.0],
                [8.0, 3.0, 1.0],
                [7.0, 4.0, 3.0],
                [2.0, 3.0, 1.0],
                [6.0, 5.0, 1.0],
                [5.0, 6.0, 3.0],
                ]"
            @test replace(coord_string, r"\s+" => "") == replace(coord_string_ref, r"\s+" => "") # compares the two strings after removing all whitespace or line breaks

            # check TOML interaction string
            opsum_32 = neighbor_interaction("SdotS", "J", flat_32)
            interaction_string = toml_interactions(opsum_32; zero_based=true)
            interaction_string_ref = "Interactions = [
                ['J', 'SdotS', 0, 8],
                ['J', 'SdotS', 0, 27],
                ['J', 'SdotS', 0, 29],
                ['J', 'SdotS', 1, 9],
                ['J', 'SdotS', 1, 27],
                ['J', 'SdotS', 1, 28],
                ['J', 'SdotS', 2, 10],
                ['J', 'SdotS', 2, 28],
                ['J', 'SdotS', 2, 30],
                ['J', 'SdotS', 3, 11],
                ['J', 'SdotS', 3, 24],
                ['J', 'SdotS', 3, 31],
                ['J', 'SdotS', 4, 12],
                ['J', 'SdotS', 4, 24],
                ['J', 'SdotS', 4, 25],
                ['J', 'SdotS', 5, 13],
                ['J', 'SdotS', 5, 30],
                ['J', 'SdotS', 5, 31],
                ['J', 'SdotS', 6, 14],
                ['J', 'SdotS', 6, 25],
                ['J', 'SdotS', 6, 26],
                ['J', 'SdotS', 7, 15],
                ['J', 'SdotS', 7, 26],
                ['J', 'SdotS', 7, 29],
                ['J', 'SdotS', 8, 16],
                ['J', 'SdotS', 8, 23],
                ['J', 'SdotS', 9, 16],
                ['J', 'SdotS', 9, 17],
                ['J', 'SdotS', 10, 17],
                ['J', 'SdotS', 10, 18],
                ['J', 'SdotS', 11, 19],
                ['J', 'SdotS', 11, 21],
                ['J', 'SdotS', 12, 19],
                ['J', 'SdotS', 12, 20],
                ['J', 'SdotS', 13, 18],
                ['J', 'SdotS', 13, 21],
                ['J', 'SdotS', 14, 20],
                ['J', 'SdotS', 14, 22],
                ['J', 'SdotS', 15, 22],
                ['J', 'SdotS', 15, 23],
                ['J', 'SdotS', 16, 24],
                ['J', 'SdotS', 17, 25],
                ['J', 'SdotS', 18, 26],
                ['J', 'SdotS', 19, 27],
                ['J', 'SdotS', 20, 28],
                ['J', 'SdotS', 21, 29],
                ['J', 'SdotS', 22, 30],
                ['J', 'SdotS', 23, 31],
                ]"
            @test replace(interaction_string, r"\s+" => "") == replace(interaction_string_ref, r"\s+" => "")
        end

    end

    # ================================================================
    # Shastry-Sutherland lattice with one open one periodic boundary
    # ================================================================
    @testset "Shastry-Sutherland cylinder" begin
        L = 6
        W = 4
        boundary = [L÷2 0; 0 W÷2]
        fl = FiniteLattice(shastry_sutherland, boundary, [false, true]; atom_order=order_xy)
        
        # check if coordinates match
        coordinate_str = toml_coordinates(fl)
        coordinate_str_expected = "Coordinates = [
            [0.0, 0.0],
            [0.0, 0.5],
            [0.0, 1.0],
            [0.0, 1.5],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.5, 1.0],
            [0.5, 1.5],
            [1.0, 0.0],
            [1.0, 0.5],
            [1.0, 1.0],
            [1.0, 1.5],
            [1.5, 0.0],
            [1.5, 0.5],
            [1.5, 1.0],
            [1.5, 1.5],
            [2.0, 0.0],
            [2.0, 0.5],
            [2.0, 1.0],
            [2.0, 1.5],
            [2.5, 0.0],
            [2.5, 0.5],
            [2.5, 1.0],
            [2.5, 1.5],
            ]"
        @test replace(coordinate_str, r"\s+" => "") == replace(coordinate_str_expected, r"\s+" => "")

        # create standard J-Jd Hamiltonian for Shastry-Sutherland cylinder
        H = OpSum()
        H += neighbor_interaction("SdotS", "J", fl; num_distance=1)
        # Jd dimer bond 1: atom 1 and atom 4 within the same unit cell
        H += lattice_interaction("SdotS", "Jd", fl, 1, 4, [0, 0])
        # Jd dimer bond 2: atom 3 (origin cell) to atom 2 (cell offset [1,-1])
        H += lattice_interaction("SdotS", "Jd", fl, 3, 2, [1, -1])
        interaction_str = toml_interactions(H; zero_based=false)
        interaction_str_expected = "Interactions = [
            ['J', 'SdotS', 1, 2],
            ['J', 'SdotS', 1, 4],
            ['J', 'SdotS', 1, 5],
            ['J', 'SdotS', 2, 3],
            ['J', 'SdotS', 2, 6],
            ['J', 'SdotS', 3, 4],
            ['J', 'SdotS', 3, 7],
            ['J', 'SdotS', 4, 8],
            ['J', 'SdotS', 5, 6],
            ['J', 'SdotS', 5, 8],
            ['J', 'SdotS', 5, 9],
            ['J', 'SdotS', 6, 7],
            ['J', 'SdotS', 6, 10],
            ['J', 'SdotS', 7, 8],
            ['J', 'SdotS', 7, 11],
            ['J', 'SdotS', 8, 12],
            ['J', 'SdotS', 9, 10],
            ['J', 'SdotS', 9, 12],
            ['J', 'SdotS', 9, 13],
            ['J', 'SdotS', 10, 11],
            ['J', 'SdotS', 10, 14],
            ['J', 'SdotS', 11, 12],
            ['J', 'SdotS', 11, 15],
            ['J', 'SdotS', 12, 16],
            ['J', 'SdotS', 13, 14],
            ['J', 'SdotS', 13, 16],
            ['J', 'SdotS', 13, 17],
            ['J', 'SdotS', 14, 15],
            ['J', 'SdotS', 14, 18],
            ['J', 'SdotS', 15, 16],
            ['J', 'SdotS', 15, 19],
            ['J', 'SdotS', 16, 20],
            ['J', 'SdotS', 17, 18],
            ['J', 'SdotS', 17, 20],
            ['J', 'SdotS', 17, 21],
            ['J', 'SdotS', 18, 19],
            ['J', 'SdotS', 18, 22],
            ['J', 'SdotS', 19, 20],
            ['J', 'SdotS', 19, 23],
            ['J', 'SdotS', 20, 24],
            ['J', 'SdotS', 21, 22],
            ['J', 'SdotS', 21, 24],
            ['J', 'SdotS', 22, 23],
            ['J', 'SdotS', 23, 24],
            ['Jd', 'SdotS', 1, 6],
            ['Jd', 'SdotS', 3, 8],
            ['Jd', 'SdotS', 9, 14],
            ['Jd', 'SdotS', 11, 16],
            ['Jd', 'SdotS', 17, 22],
            ['Jd', 'SdotS', 19, 24],
            ['Jd', 'SdotS', 5, 12],
            ['Jd', 'SdotS', 7, 10],
            ['Jd', 'SdotS', 13, 20],
            ['Jd', 'SdotS', 15, 18],
        ]"
        @test replace(interaction_str, r"\s+" => "") == replace(interaction_str_expected, r"\s+" => "")
    end
end