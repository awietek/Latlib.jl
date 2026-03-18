# This file tests the functionalities of the `FiniteLattice` class defined inside src/lattice/.


@testset "Tests for FiniteLattice class" begin

    # Test of different constructors
    @testset "Constructors" begin
        
        # Create basic lattices for testing
        square_A = [1.0 0.0; 0.0 1.0]
        square_lattice = Lattice(square_A)
        
        cubic_A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        cubic_lattice = Lattice(cubic_A)
        
        # Test basic constructor with matrix boundary and periodicity vector
        @testset "Basic Matrix Constructor" begin
            boundary_2d = [3 0; 0 3]
            periodicity_2d = [true, true]
            
            flattice = FiniteLattice(square_lattice, boundary_2d, periodicity_2d)
            
            @test isa(flattice, FiniteLattice)
            @test flattice.lattice == square_lattice
            @test flattice.boundary == boundary_2d
            @test flattice.periodicity == periodicity_2d
            @test flattice.bravais_order == isless  # default order function
            @test flattice.atom_order === nothing  # default atom_order is nothing
            
            # Test with custom order functions
            custom_order(a, b) = a > b
            flattice_ordered = FiniteLattice(square_lattice, boundary_2d, periodicity_2d; bravais_order=custom_order)
            @test flattice_ordered.bravais_order == custom_order
            @test flattice_ordered.atom_order === nothing
            
        end
        
        
        
        # Test constructor with uniform periodicity (scalar Bool)
        @testset "Uniform Periodicity Constructor" begin
            boundary_2d = [4 0; 0 4]
            
            # Test with periodic=true (default)
            flattice_periodic = FiniteLattice(square_lattice, boundary_2d)
            @test flattice_periodic.periodicity == [true, true]
            
            # Test with periodic=false
            flattice_open = FiniteLattice(square_lattice, boundary_2d, false)
            @test flattice_open.periodicity == [false, false]
            
            # Test 3D case
            boundary_3d = [2 0 0; 0 2 0; 0 0 2]
            flattice_3d = FiniteLattice(cubic_lattice, boundary_3d, true)
            @test flattice_3d.periodicity == [true, true, true]
        end
           
        # Test constructor with LatticeVector boundaries
        @testset "LatticeVector Boundary Constructor" begin
            # Create boundary vectors as LatticeVectors
            v1 = LatticeVector(square_lattice, [3.0, 0.0])
            v2 = LatticeVector(square_lattice, [0.0, 3.0])
            boundary_vectors = [v1, v2]
            periodicity_2d = [true, false]
            
            flattice = FiniteLattice(boundary_vectors, periodicity_2d)
            @test isa(flattice, FiniteLattice)
            @test flattice.lattice == square_lattice
            @test flattice.boundary == [3 0; 0 3]
            @test flattice.periodicity == periodicity_2d
            
            # Test with uniform periodicity
            flattice_uniform = FiniteLattice(boundary_vectors, true)
            @test flattice_uniform.periodicity == [true, true]
        end
        
        # Test integer boundary matrix conversion
        @testset "Integer Matrix Boundary" begin
            boundary_int = Int64[2 0; 0 2]
            periodicity_2d = [true, true]
            
            flattice = FiniteLattice(square_lattice, boundary_int, periodicity_2d)
            @test isa(flattice, FiniteLattice)
            @test eltype(flattice.boundary) == Int64  # Should remain as Int64
        end
        
        
        
    end
    
    
    @testset "Constructor Error Cases" begin
        
        square_A = [1.0 0.0; 0.0 1.0]
        square_lattice = Lattice(square_A)
        cubic_A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        cubic_lattice = Lattice(cubic_A)
        
        # Test dimension mismatch in boundary matrix
        @testset "Boundary Matrix Dimension Errors" begin
            wrong_boundary_2x3 = [1 0 0; 0 1 0]  # 2x3 matrix for 2D lattice
            wrong_boundary_3x2 = [1 0; 0 1; 0 0]  # 3x2 matrix for 2D lattice
            periodicity_2d = [true, true]
            
            @test_throws ErrorException FiniteLattice(square_lattice, wrong_boundary_2x3, periodicity_2d)
            @test_throws ErrorException FiniteLattice(square_lattice, wrong_boundary_3x2, periodicity_2d)
            
            # Test 3D case with wrong dimensions
            wrong_boundary_3x2_for_3d = [1 0; 0 1; 0 0]
            periodicity_3d = [true, true, true]
            @test_throws ErrorException FiniteLattice(cubic_lattice, wrong_boundary_3x2_for_3d, periodicity_3d)
        end
        
        # Test periodicity vector length mismatch  
        @testset "Periodicity Vector Length Errors" begin
            boundary_2d = [3 0; 0 3]
            wrong_periodicity_3d = [true, true, true]  # 3D periodicity for 2D lattice
            wrong_periodicity_1d = [true]              # 1D periodicity for 2D lattice
            
            @test_throws ErrorException FiniteLattice(square_lattice, boundary_2d, wrong_periodicity_3d)
            @test_throws ErrorException FiniteLattice(square_lattice, boundary_2d, wrong_periodicity_1d)
        end
        
        # Test LatticeVector boundary errors
        @testset "LatticeVector Boundary Errors" begin
            # Create vectors from different lattices
            v1_square = LatticeVector(square_lattice, [2.0, 0.0])
            v2_cubic = LatticeVector(cubic_lattice, [0.0, 2.0, 0.0])
            mixed_boundary = [v1_square, v2_cubic]
            periodicity_2d = [true, true]
            
            @test_throws ErrorException FiniteLattice(mixed_boundary, periodicity_2d)
            
            # Test non-integer LatticeVector coordinates (not in lattice)
            v_non_integer = LatticeVector(square_lattice, [1.5, 0.0])
            boundary_non_integer = [v_non_integer, LatticeVector(square_lattice, [0.0, 2.0])]
            @test_throws ErrorException FiniteLattice(boundary_non_integer, periodicity_2d)
            
            # Test periodicity length mismatch with LatticeVector boundary
            v1 = LatticeVector(square_lattice, [2.0, 0.0])
            v2 = LatticeVector(square_lattice, [0.0, 2.0])
            boundary_vectors = [v1, v2]
            wrong_periodicity = [true, true, true]  # 3D periodicity for 2D vectors
            
            @test_throws ErrorException FiniteLattice(boundary_vectors, wrong_periodicity)
        end
        
    end
    
    @testset "Properties and Methods" begin
        
        square_A = [1.0 0.0; 0.0 1.0]
        square_lattice = Lattice(square_A)
        boundary_2d = [3 0; 0 3]
        periodicity_2d = [true, true]
        
        flattice = FiniteLattice(square_lattice, boundary_2d, periodicity_2d)
        
        # Test dimension function
        @test dim(flattice) == 2
        @test dim(flattice) == square_lattice.dim
        
        # Test natoms function  
        @test natoms(flattice) == natoms(square_lattice)
        @test natoms(flattice) == 1  # square lattice has 1 atom
        
    end
    
    

end