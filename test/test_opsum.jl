# This file tests the functionalities of the `Op`, `OpSum`, `neighbor_interaction`,
# and `lattice_interaction` defined in src/opsum.jl.


@testset "Op struct" begin

    @testset "Valid Constructions" begin
        # String coupling constant
        op = Op("HB", "J", [1, 2])
        @test isa(op, Op)
        @test op.type == "HB"
        @test op.cpl == "J"
        @test op.sites == [1, 2]

        # Float coupling constant
        op_f = Op("HB", 1.5, [1, 2])
        @test isa(op_f, Op)
        @test op_f.cpl == 1.5

        # Integer coupling constant
        op_i = Op("HB", 2, [3, 4])
        @test isa(op_i, Op)
        @test op_i.cpl == 2

        # Different operator types
        @test Op("Cup", "t", [5, 6]).type == "Cup"
        @test Op("Cdn", "t", [7, 8]).type == "Cdn"

        # Multi-site operator (e.g. ring exchange)
        op_ring = Op("Ring", "K", [1, 2, 3, 4])
        @test length(op_ring.sites) == 4
    end

    @testset "Invalid Constructions" begin
        # Wrong type for `type` argument
        @test_throws MethodError Op(123, "J", [1, 2])
        # Wrong type for `cpl` argument
        @test_throws MethodError Op("HB", [1, 2], [1, 2])
        # Wrong type for `sites` argument
        @test_throws MethodError Op("HB", "J", "not a vector")

        # Empty sites vector
        @test_throws ErrorException Op("HB", "J", Int64[])

        # Duplicate sites
        @test_throws ErrorException Op("HB", "J", [1, 1])
        @test_throws ErrorException Op("HB", "J", [3, 5, 3])

        # Non-positive site indices
        @test_throws ErrorException Op("HB", "J", [0, 1])
        @test_throws ErrorException Op("HB", "J", [-1, 2])
        @test_throws ErrorException Op("HB", "J", [-3])
    end

    @testset "Equality" begin
        op1 = Op("HB", "J1", [1, 2])
        op2 = Op("HB", "J1", [1, 2])
        op3 = Op("HB", "J2", [1, 2])
        op4 = Op("HB", "J1", [2, 3])
        op5 = Op("Cup", "J1", [1, 2])
        op6 = Op("HB", "J1", [2, 1])

        # Identical operators are equal
        @test op1 == op2
        @test isequal(op1, op2)

        # Operators differing in any field are not equal
        @test op1 != op3   # different coupling
        @test op1 != op4   # different sites
        @test op1 != op5   # different type
        @test op1 != op6   # different site order

        # String vs numeric coupling are different
        @test Op("HB", "1", [1, 2]) != Op("HB", 1, [1, 2])
    end

    @testset "Hashing" begin
        op1 = Op("HB", "J1", [1, 2])
        op2 = Op("HB", "J1", [1, 2])
        # Equal operators must have equal hashes
        @test hash(op1) == hash(op2)
    end
end


@testset "OpSum struct" begin

    @testset "Default (empty) Constructor" begin
        opsum = OpSum()
        @test isa(opsum, OpSum)
        @test isempty(opsum.ops)
    end

    @testset "Constructor from Vector{Op}" begin
        ops = [Op("HB", "J", [1, 2]), Op("HB", "J", [2, 3])]
        opsum = OpSum(ops)
        @test length(opsum.ops) == 2
        @test opsum.ops[1] == ops[1]
        @test opsum.ops[2] == ops[2]
    end

    @testset "Invalid Constructor Arguments" begin
        @test_throws MethodError OpSum("not a vector")
        @test_throws MethodError OpSum(42)
    end

    @testset "Op addition (OpSum + Op)" begin
        opsum = OpSum()
        opsum = opsum + Op("HB", "J", [1, 2, 3])
        @test length(opsum.ops) == 1

        opsum = opsum + Op("HB", "J", [4])
        @test length(opsum.ops) == 2

        # Adding an Op with empty sites must fail
        @test_throws ErrorException OpSum() + Op("HB", "J", Int64[])
    end

    @testset "OpSum addition (OpSum + OpSum)" begin
        os1 = OpSum() + Op("HB", "J1", [1, 2])
        os2 = OpSum() + Op("HB", "J2", [3])
        os3 = os1 + os2
        @test length(os3.ops) == 2
    end

    @testset "unique_ops!" begin
        opsum = OpSum()
        opsum = opsum + Op("HB", "J", [1, 2])
        opsum = opsum + Op("HB", "J", [1, 2])  # duplicate
        opsum = opsum + Op("HB", "J", [2, 1])  # not necessarily a duplicate
        opsum = opsum + Op("HB", "J", [3, 4])

        @test length(opsum.ops) == 4
        unique_ops!(opsum)
        @test length(opsum.ops) == 3

        # After unique_ops!, calling it again should not change anything
        unique_ops!(opsum)
        @test length(opsum.ops) == 3
    end

end




@testset "neighbor_interaction" begin

    # ================================================================
    # Setup: 2D square lattice (1 atom per cell)
    # ================================================================
    square_A = [1.0 0.0; 0.0 1.0]
    square_lattice = Lattice(square_A)

    # 3×3 periodic square lattice → 9 sites
    fl_sq = FiniteLattice(square_lattice, [3 0; 0 3], true)
    nsites_sq = 9

    # ================================================================
    # Setup: 3D cubic lattice (1 atom per cell)
    # ================================================================
    cubic_A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    cubic_lattice = Lattice(cubic_A)

    # 3×3×3 periodic cubic lattice → 27 sites
    fl_cb = FiniteLattice(cubic_lattice, [3 0 0; 0 3 0; 0 0 3], true)
    nsites_cb = 27


    @testset "2D Square – Nearest Neighbors (num_distance=1)" begin
        ops = neighbor_interaction("HB", "J", fl_sq; num_distance=1)
        @test isa(ops, OpSum)

        # Each of the 9 sites has 4 NN (PBC) → 9*4/2 = 18 unique bonds
        @test length(ops.ops) == 2 * nsites_sq

        for op in ops.ops
            @test op.type == "HB"
            @test op.cpl == "J"
            @test length(op.sites) == 2
            @test op.sites[1] < op.sites[2]  # ordered pairs
        end
    end

    @testset "2D Square – Next-Nearest Neighbors (num_distance=2)" begin
        ops = neighbor_interaction("HB", "Jd", fl_sq; num_distance=2)
        @test isa(ops, OpSum)

        # Each of the 9 sites has 4 NNN (diagonals with PBC) → 9*4/2 = 18
        @test length(ops.ops) == 2 * nsites_sq

        for op in ops.ops
            @test op.type == "HB"
            @test op.cpl == "Jd"
        end
    end

    @testset "2D Square – No Duplicate Bonds" begin
        ops = neighbor_interaction("HB", "J", fl_sq; num_distance=1)
        n_before = length(ops.ops)
        unique_ops!(ops)
        @test length(ops.ops) == n_before
    end

    @testset "3D Cubic – Nearest Neighbors (num_distance=1)" begin
        ops = neighbor_interaction("HB", "J", fl_cb; num_distance=1)
        @test isa(ops, OpSum)

        # Each of 27 sites has 6 NN (PBC) → 27*6/2 = 27*3 = 81 unique bonds
        @test length(ops.ops) == 81

        for op in ops.ops
            @test op.type == "HB"
            @test op.cpl == "J"
            @test length(op.sites) == 2
            @test op.sites[1] < op.sites[2]
        end
    end

    @testset "3D Cubic – Next-Nearest Neighbors (num_distance=2)" begin
        ops = neighbor_interaction("HB", "Jd", fl_cb; num_distance=2)
        @test isa(ops, OpSum)

        # Each of 27 sites has 12 NNN (face diagonals with PBC) → 27*12/2 = 162
        @test length(ops.ops) == 162
    end

    @testset "3D Cubic – No Duplicate Bonds" begin
        ops = neighbor_interaction("HB", "J", fl_cb; num_distance=1)
        n_before = length(ops.ops)
        unique_ops!(ops)
        @test length(ops.ops) == n_before
    end

    @testset "Invalid num_distance" begin
        @test_throws ErrorException neighbor_interaction("HB", "J", fl_sq; num_distance=0)
        @test_throws ErrorException neighbor_interaction("HB", "J", fl_sq; num_distance=-1)
    end

end





@testset "lattice_interaction" begin

    # ================================================================
    # Setup: 2D square lattice
    # ================================================================
    square_A = [1.0 0.0; 0.0 1.0]
    square_lattice = Lattice(square_A)
    fl_sq = FiniteLattice(square_lattice, [3 0; 0 3], true)
    nsites_sq = 9

    # ================================================================
    # Setup: 3D cubic lattice
    # ================================================================
    cubic_A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    cubic_lattice = Lattice(cubic_A)
    fl_cb = FiniteLattice(cubic_lattice, [3 0 0; 0 3 0; 0 0 3], true)
    nsites_cb = 27


    # ----- 2D square lattice -----
    @testset "2D Square – NN along x: cell2=[1,0]" begin
        ops = lattice_interaction("HB", "J", fl_sq, 1, 1, [1, 0])
        @test isa(ops, OpSum)
        # One bond per Bravais cell → 9 operators
        @test length(ops.ops) == nsites_sq

        for op in ops.ops
            @test op.type == "HB"
            @test op.cpl == "J"
            @test length(op.sites) == 2
            @test op.sites[1] < op.sites[2]
        end
    end

    @testset "2D Square – NN along y: cell2=[0,1]" begin
        ops = lattice_interaction("HB", "J", fl_sq, 1, 1, [0, 1])
        @test length(ops.ops) == nsites_sq
    end

    @testset "2D Square – NNN diagonal: cell2=[1,1]" begin
        ops = lattice_interaction("HB", "Jd", fl_sq, 1, 1, [1, 1])
        @test length(ops.ops) == nsites_sq
    end

    @testset "2D Square – NNN anti-diagonal: cell2=[1,-1]" begin
        ops = lattice_interaction("HB", "Jd", fl_sq, 1, 1, [1, -1])
        @test length(ops.ops) == nsites_sq
    end

    @testset "2D Square – No Duplicate Bonds" begin
        ops = lattice_interaction("HB", "J", fl_sq, 1, 1, [1, 0])
        n_before = length(ops.ops)
        unique_ops!(ops)
        @test length(ops.ops) == n_before
    end


    # ----- 3D cubic lattice -----
    @testset "3D Cubic – NN along each axis" begin
        for cell2 in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            ops = lattice_interaction("HB", "J", fl_cb, 1, 1, cell2)
            @test isa(ops, OpSum)
            @test length(ops.ops) == nsites_cb

            for op in ops.ops
                @test op.type == "HB"
                @test op.cpl == "J"
                @test length(op.sites) == 2
                @test op.sites[1] < op.sites[2]
            end
        end
    end

    @testset "3D Cubic – NNN along face diagonals" begin
        for cell2 in [[1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]]
            ops = lattice_interaction("HB", "Jd", fl_cb, 1, 1, cell2)
            @test length(ops.ops) == nsites_cb
        end
    end

    @testset "3D Cubic – No Duplicate Bonds" begin
        ops = lattice_interaction("HB", "J", fl_cb, 1, 1, [1, 0, 0])
        n_before = length(ops.ops)
        unique_ops!(ops)
        @test length(ops.ops) == n_before
    end


    # ----- LatticeVector argument for cell2 -----
    @testset "LatticeVector cell2 gives same result as Vector{Int}" begin
        cell2_lv = LatticeVector(square_lattice, [1.0, 0.0])
        ops_lv = lattice_interaction("HB", "J", fl_sq, 1, 1, cell2_lv)

        cell2_arr = [1, 0]
        ops_arr = lattice_interaction("HB", "J", fl_sq, 1, 1, cell2_arr)

        @test length(ops_lv.ops) == length(ops_arr.ops)
        @test Set(ops_lv.ops) == Set(ops_arr.ops)
    end


    # ----- Invalid arguments -----
    @testset "Invalid atom indices" begin
        # atom index 0 is out of bounds (1-indexed)
        @test_throws ErrorException lattice_interaction("HB", "J", fl_sq, 0, 1, [1, 0])
        # square lattice has only 1 atom per unit cell
        @test_throws ErrorException lattice_interaction("HB", "J", fl_sq, 1, 2, [1, 0])
        @test_throws ErrorException lattice_interaction("HB", "J", fl_sq, 2, 1, [1, 0])

        # same for cubic
        @test_throws ErrorException lattice_interaction("HB", "J", fl_cb, 0, 1, [1, 0, 0])
        @test_throws ErrorException lattice_interaction("HB", "J", fl_cb, 1, 2, [1, 0, 0])
    end

end





@testset "neighbor_interaction vs lattice_interaction" begin

    # ================================================================
    # 2D Square Lattice
    # ================================================================
    @testset "2D Square Lattice" begin
        square_A = [1.0 0.0; 0.0 1.0]
        square_lattice = Lattice(square_A)
        fl = FiniteLattice(square_lattice, [3 0; 0 3], true)

        @testset "NN bonds match" begin
            # Build NN via neighbor_interaction
            ops_nb = neighbor_interaction("HB", "J", fl; num_distance=1)

            # Build NN via lattice_interaction: x-direction + y-direction
            ops_li = lattice_interaction("HB", "J", fl, 1, 1, [1, 0]) +
                      lattice_interaction("HB", "J", fl, 1, 1, [0, 1])

            @test length(ops_nb.ops) == length(ops_li.ops)
            @test Set(ops_nb.ops) == Set(ops_li.ops)
        end

        @testset "NNN bonds match" begin
            # Build NNN via neighbor_interaction
            ops_nb = neighbor_interaction("HB", "Jd", fl; num_distance=2)

            # Build NNN via lattice_interaction: two diagonal directions
            ops_li = lattice_interaction("HB", "Jd", fl, 1, 1, [1, 1]) +
                      lattice_interaction("HB", "Jd", fl, 1, 1, [1, -1])

            @test length(ops_nb.ops) == length(ops_li.ops)
            @test Set(ops_nb.ops) == Set(ops_li.ops)
        end
    end

    # ================================================================
    # 3D Cubic Lattice
    # ================================================================
    @testset "3D Cubic Lattice" begin
        cubic_A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        cubic_lattice = Lattice(cubic_A)
        fl = FiniteLattice(cubic_lattice, [3 0 0; 0 3 0; 0 0 3], true)

        @testset "NN bonds match" begin
            ops_nb = neighbor_interaction("HB", "J", fl; num_distance=1)

            # Build NN via lattice_interaction: three axis directions
            ops_li = lattice_interaction("HB", "J", fl, 1, 1, [1, 0, 0]) +
                      lattice_interaction("HB", "J", fl, 1, 1, [0, 1, 0]) +
                      lattice_interaction("HB", "J", fl, 1, 1, [0, 0, 1])

            @test length(ops_nb.ops) == length(ops_li.ops)
            @test Set(ops_nb.ops) == Set(ops_li.ops)
        end

        @testset "NNN bonds match" begin
            ops_nb = neighbor_interaction("HB", "Jd", fl; num_distance=2)

            # Build NNN via lattice_interaction: six face-diagonal directions
            ops_li = lattice_interaction("HB", "Jd", fl, 1, 1, [1, 1, 0]) +
                      lattice_interaction("HB", "Jd", fl, 1, 1, [1, -1, 0]) +
                      lattice_interaction("HB", "Jd", fl, 1, 1, [1, 0, 1]) +
                      lattice_interaction("HB", "Jd", fl, 1, 1, [1, 0, -1]) +
                      lattice_interaction("HB", "Jd", fl, 1, 1, [0, 1, 1]) +
                      lattice_interaction("HB", "Jd", fl, 1, 1, [0, 1, -1])

            @test length(ops_nb.ops) == length(ops_li.ops)
            @test Set(ops_nb.ops) == Set(ops_li.ops)
        end
    end


end

