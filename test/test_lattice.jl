# This file tests the functionalities of the `Lattice` class defined inside src/lattice/.


@testset "Tests for Lattice class" begin

    # Test of different constructors
    @testset "Constructors" begin
        
        # simple square lattice
        A = [
            1.0 0.0
            0.0 1.0
        ]
        @test isa(Lattice(A), Lattice)

        # simple cubic lattice
        A = [
            1.0 0.0 0.0;            # a1
            0.0 1.0 0.0;            # a2
            0.0 0.0 1.0             # a3
        ]
        @test isa(Lattice(A), Lattice)

        # 2D lattice with 3 vectors ()
        A = [
            1.0 0.0
            0.0 1.0
            1.0 1.0
        ]
        @test_throws ErrorException Lattice(A)

        # 3D lattice with 2 vectors (MUST THROW)
        A = [
            1.0 0.0 0.0
            0.0 1.0 0.0
        ]
        @test_throws ErrorException Lattice(A)

        # simple square lattice with 2 atoms
        A = [
            1.0 0.0
            0.0 1.0
        ]
        positions = [
            - 0.25 0.0;                # atom 1
            + 0.25 0.0                 # atom 2
        ]
        @test isa(Lattice(A, positions), Lattice)

        # square lattice with 3D positions (MUST THROW)
        positions = [
            0.0 0.0 0.0;             # atom 1
            0.5 0.0 0.0              # atom 2
        ]
        @test_throws ErrorException Lattice(A, positions)

        # simple cubic lattice with 3 atoms and non-default types
        A = [
            1.0 0.0 0.0;            # a1
            0.0 1.0 0.0;            # a2
            0.0 0.0 1.0             # a3
        ]
        positions = [
            0.0 0.0 0.0;             # atom 1
            -0.25 -0.25 -0.25;       # atom 2
            0.25 0.25 0.25           # atom 3
        ]
        types = [1, 2] # only two types (MUST THROW)
        @test_throws ErrorException Lattice(A, positions; types = types)
        @test isa(Lattice(A, positions; types = [1, 2, 1]), Lattice)

        # simple square lattice where atom position contain duplicates (MUST THROW)
        A = [
            1.0 0.0
            0.0 1.0
        ]
        positions = [
            0.0 0.0;                # atom 1
            0.0 0.0                 # atom 2 (duplicate of atom 1)
        ]
        @test_throws ErrorException Lattice(A, positions)

        # check constructor via `EuclideanVector` type for lattice vectors and/or for positions
        A = [
            1.0 0.0
            0.0 1.0
        ]
        vs = [EuclideanVector([1.0, 0.0]), EuclideanVector([0.0, 1.0])]
        positions = [
            0.0 0.1;
            0.2 0.3;
        ]
        ps = [EuclideanVector([0.0, 0.1]), EuclideanVector([0.2, 0.3])]
        @test (  Lattice(A, positions) == Lattice(vs, positions) == Lattice(A, ps) == Lattice(vs, ps)  )

    end




    # Test for transforming vectors between Euclidean and lattice basis
    @testset "Basis transformations" begin
        # "cubic" lattice with scaled a1 and a3 vectors
        A = [
            2.0 0.0 0.0;    # a1    
            0.0 1.0 0.0;   # a2
            0.0 0.0 5.0    # a3
        ]
        lat = Lattice(A)

        # check that the standard vectors in the lattice basis correspond to a1, a2, a3
        for i in 1:3
            v_lat_coords = zeros(3)
            v_lat_coords[i] = 1.0
            v_lat = LatticeVector(lat, v_lat_coords)
            v_euc = to_euclidean_basis(v_lat)
            @test isapprox(v_euc.coords, A[:, i]) # do we get a_i?
            @test isapprox(to_lattice_basis(lat, v_euc).coords, v_lat_coords) # do we get v_lat_coords back?
        end

        # check some other vector
        v_euc = EuclideanVector([1.0, 0.5, 2.0])
        v_lat = to_lattice_basis(lat, v_euc)
        @test isapprox(v_lat.coords, [0.5, 0.5, 0.4])
        v_euc_back = to_euclidean_basis(LatticeVector(lat, [0.5, 0.5, 0.4]))
        @test v_euc_back.coords == v_euc.coords

        # check additivity
        v0 = v_euc - v_euc_back
        @test isapprox(v0.coords, [0.0, 0.0, 0.0])
        v_lat_2 = v_lat + v_lat
        @test isapprox(v_lat_2.coords, 2.0 * v_lat.coords)
        v_euc_2 = to_euclidean_basis(v_lat_2)
        @test isapprox(v_euc_2.coords, 2.0 * v_euc.coords)
        

    end
end