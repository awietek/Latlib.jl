

# This file tests the functionalities of the metric module defined in src/metric.jl.

@testset "Tests for Metric module" begin

    # ============================================================
    # Setup: 2D square lattice with non-trivial boundary [4 1; 0 3]
    # ============================================================
    square_A = [1.0 0.0; 0.0 1.0]
    square_lattice = Lattice(square_A)

    # Non-trivial (sheared) 2D boundary: b1 = 4*a1 + 1*a2, b2 = 3*a2
    boundary_2d = [4 1; 0 3]
    flattice_2d = FiniteLattice(square_lattice, boundary_2d, [true, true])

    # Simpler diagonal boundary for sanity checks
    boundary_2d_diag = [4 0; 0 4]
    flattice_2d_diag = FiniteLattice(square_lattice, boundary_2d_diag, [true, true])

    # ============================================================
    # Setup: 3D cubic lattice with non-trivial boundary [3 1 0; 0 3 1; 0 0 3]
    # ============================================================
    cubic_A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    cubic_lattice = Lattice(cubic_A)

    boundary_3d = [3 1 0; 0 3 1; 0 0 3]
    flattice_3d = FiniteLattice(cubic_lattice, boundary_3d, [true, true, true])

    boundary_3d_diag = [3 0 0; 0 3 0; 0 0 3]
    flattice_3d_diag = FiniteLattice(cubic_lattice, boundary_3d_diag, [true, true, true])

    # ============================================================
    # EuclideanMetric tests
    # ============================================================
    @testset "EuclideanMetric" begin

        @testset "2D" begin
            metric = EuclideanMetric()
            x = EuclideanVector([0.0, 0.0])
            y = EuclideanVector([3.0, 4.0])
            @test metric(x, y) ≈ 5.0  # 3-4-5 triangle

            # distance is symmetric
            @test metric(x, y) ≈ metric(y, x)

            # distance to self is zero
            @test metric(x, x) ≈ 0.0

            # unit vectors
            e1 = EuclideanVector([1.0, 0.0])
            e2 = EuclideanVector([0.0, 1.0])
            @test metric(e1, e2) ≈ sqrt(2.0)
        end

        @testset "3D" begin
            metric = EuclideanMetric()
            x = EuclideanVector([1.0, 2.0, 3.0])
            y = EuclideanVector([4.0, 6.0, 3.0])
            @test metric(x, y) ≈ 5.0  # sqrt(9 + 16 + 0)

            @test metric(x, y) ≈ metric(y, x)
            @test metric(x, x) ≈ 0.0
        end
    end

    

    # ============================================================
    # PeriodicEuclideanMetric tests
    # ============================================================
    @testset "PeriodicEuclideanMetric" begin

        @testset "2D diagonal boundary" begin
            metric = PeriodicEuclideanMetric(flattice_2d_diag)

            origin = EuclideanVector([0.0, 0.0])

            # Points separated by exactly one boundary vector should have distance 0
            b1_vec = EuclideanVector([4.0, 0.0])
            @test metric(origin, b1_vec) ≈ 0.0 atol=1e-10

            b2_vec = EuclideanVector([0.0, 4.0])
            @test metric(origin, b2_vec) ≈ 0.0 atol=1e-10

            # Nearest neighbor: distance 1 along axis
            nn = EuclideanVector([1.0, 0.0])
            @test metric(origin, nn) ≈ 1.0

            # Wrapping: point at (3.5, 0) should be closer via wrap than directly
            # Direct distance = 3.5; wrapping distance = 4.0 - 3.5 = 0.5
            x_wrap = EuclideanVector([3.5, 0.0])
            @test metric(origin, x_wrap) ≈ 0.5 atol=1e-10

            # Symmetry check
            a = EuclideanVector([1.0, 2.0])
            b = EuclideanVector([3.0, 3.5])
            @test metric(a, b) ≈ metric(b, a)
        end

        @testset "2D sheared boundary" begin
            metric = PeriodicEuclideanMetric(flattice_2d)

            origin = EuclideanVector([0.0, 0.0])

            # Boundary vectors themselves map to distance 0
            # b1 = 4*a1 + 1*a2 = (4, 1) in Euclidean coords for square lattice
            b1_euc = EuclideanVector([4.0, 1.0])
            @test metric(origin, b1_euc) ≈ 0.0 atol=1e-10

            # b2 = 3*a2 = (0, 3)
            b2_euc = EuclideanVector([0.0, 3.0])
            @test metric(origin, b2_euc) ≈ 0.0 atol=1e-10

            # Distance should be non-negative
            x = EuclideanVector([1.5, 0.7])
            y = EuclideanVector([3.2, 2.1])
            @test metric(x, y) >= 0.0

            # Periodic distance should be ≤ Euclidean distance
            euc_metric = EuclideanMetric()
            @test metric(x, y) <= euc_metric(x, y) + 1e-10

            # Self-distance is zero
            @test metric(x, x) ≈ 0.0 atol=1e-10

            # Symmetry
            @test metric(x, y) ≈ metric(y, x)
        end

        @testset "3D diagonal boundary" begin
            metric = PeriodicEuclideanMetric(flattice_3d_diag)

            origin = EuclideanVector([0.0, 0.0, 0.0])

            # Boundary vectors map to distance 0
            @test metric(origin, EuclideanVector([3.0, 0.0, 0.0])) ≈ 0.0 atol=1e-10
            @test metric(origin, EuclideanVector([0.0, 3.0, 0.0])) ≈ 0.0 atol=1e-10
            @test metric(origin, EuclideanVector([0.0, 0.0, 3.0])) ≈ 0.0 atol=1e-10

            # Wrapping in 3D: point at (2.8, 0, 0) wraps to distance 0.2
            x_wrap = EuclideanVector([2.8, 0.0, 0.0])
            @test metric(origin, x_wrap) ≈ 0.2 atol=1e-10

            # NN distance = 1.0
            @test metric(origin, EuclideanVector([1.0, 0.0, 0.0])) ≈ 1.0
        end


        # ----- ATTENTION -----
        # this test is the one that can catch the "error"
        # in the current implementation of distance_vector() inside metric.jl!
        # (explained therein)
        @testset "3D sheared boundary" begin

            metric = PeriodicEuclideanMetric(flattice_3d)

            origin = EuclideanVector([0.0, 0.0, 0.0])

            # b1 = (3, 1, 0), b2 = (0, 3, 1), b3 = (0, 0, 3) for cubic lattice
            @test metric(origin, EuclideanVector([3.0, 1.0, 0.0])) ≈ 0.0 atol=1e-10
            @test metric(origin, EuclideanVector([0.0, 3.0, 1.0])) ≈ 0.0 atol=1e-10
            @test metric(origin, EuclideanVector([0.0, 0.0, 3.0])) ≈ 0.0 atol=1e-10

            # Periodic distance ≤ Euclidean distance
            euc_metric = EuclideanMetric()
            x = EuclideanVector([0.5, 1.2, 2.1])
            y = EuclideanVector([2.3, 0.4, 0.8])
            @test metric(x, y) <= euc_metric(x, y) + 1e-10

            # Symmetry
            @test metric(x, y) ≈ metric(y, x)
        end

        @testset "LatticeVector dispatch" begin
            metric = PeriodicEuclideanMetric(flattice_2d_diag)

            # PeriodicEuclideanMetric should also accept LatticeVectors
            v1 = LatticeVector(square_lattice, [0.0, 0.0])
            v2 = LatticeVector(square_lattice, [1.0, 0.0])
            @test metric(v1, v2) ≈ 1.0

            # Wrapping: v at (3.5, 0) in lattice coords = (3.5, 0) in Euclidean for square
            v_wrap = LatticeVector(square_lattice, [3.5, 0.0])
            v_origin = LatticeVector(square_lattice, [0.0, 0.0])
            @test metric(v_origin, v_wrap) ≈ 0.5 atol=1e-10
        end
    end

    

    # ============================================================
    # distance() convenience function
    # ============================================================
    @testset "distance()" begin

        @testset "Euclidean (no flattice)" begin
            x = EuclideanVector([0.0, 0.0])
            y = EuclideanVector([3.0, 4.0])
            @test distance(x, y) ≈ 5.0
            @test distance(x, y; flattice=nothing) ≈ 5.0
        end

        @testset "Periodic (with flattice)" begin
            origin = EuclideanVector([0.0, 0.0])
            nn = EuclideanVector([1.0, 0.0])
            @test distance(origin, nn; flattice=flattice_2d_diag) ≈ 1.0

            # Wrapping
            x_wrap = EuclideanVector([3.7, 0.0])
            @test distance(origin, x_wrap; flattice=flattice_2d_diag) ≈ 0.3 atol=1e-10

            # 3D
            origin_3d = EuclideanVector([0.0, 0.0, 0.0])
            x_3d = EuclideanVector([2.8, 0.0, 0.0])
            @test distance(origin_3d, x_3d; flattice=flattice_3d_diag) ≈ 0.2 atol=1e-10
        end
    end

    

    # ============================================================
    # distance_vector() tests
    # ============================================================
    @testset "distance_vector()" begin

        using LinearAlgebra

        @testset "Euclidean" begin
            x = EuclideanVector([1.0, 2.0])
            y = EuclideanVector([4.0, 6.0])
            dv = distance_vector(x, y)
            @test dv.coords ≈ [3.0, 4.0]
        end

        @testset "Periodic wrapping 2D" begin
            origin = EuclideanVector([0.0, 0.0])
            far = EuclideanVector([3.5, 0.0])
            dv = distance_vector(origin, far; flattice=flattice_2d_diag)
            # The periodic image should give a shorter vector
            @test LinearAlgebra.norm(dv.coords) ≈ distance(origin, far; flattice=flattice_2d_diag) atol=1e-10
        end
    end

    

    # ============================================================
    # Nearest neighbors
    # ============================================================
    @testset "Nearest neighbors" begin

        # -------------------------------------------------------
        # 2D square lattice, 4×4 diagonal boundary, fully periodic
        # -------------------------------------------------------
        @testset "2D square 4×4 periodic NN" begin
            # Atom coordinates on a 4×4 periodic square lattice
            # For a 1-atom square lattice with boundary [4 0; 0 4],
            # there are 16 sites at integer coords (0..3) x (0..3).
            sites_2d = EuclideanVector[]
            for ix in 0:3, iy in 0:3
                push!(sites_2d, EuclideanVector([Float64(ix), Float64(iy)]))
            end

            # Unique distances on 4×4 periodic square lattice
            dists_2d = distances(sites_2d; flattice=flattice_2d_diag)

            # Distance 0 should be the first entry (self-distance)
            @test dists_2d[1] ≈ 0.0

            # NN distance on square lattice = 1.0
            @test dists_2d[2] ≈ 1.0

            # 2nd NN distance = sqrt(2) (diagonal)
            @test dists_2d[3] ≈ sqrt(2.0) atol=1e-10

            # 3rd NN distance = 2.0
            @test dists_2d[4] ≈ 2.0 atol=1e-10

            # --- NN bonds ---
            nn_2d = neighbors(sites_2d; num_distance=1, flattice=flattice_2d_diag)

            # 16-sites with 4 NN each => 16 * 4 / 2 unique pairs
            @test length(nn_2d) == 32

            # Every pair should have distance 1.0
            metric_2d = PeriodicEuclideanMetric(flattice_2d_diag)
            for (i, j) in nn_2d
                @test metric_2d(sites_2d[i], sites_2d[j]) ≈ 1.0
                @test i < j  # pairs are ordered
            end

            # --- 2nd NN bonds (diagonal neighbors) ---
            nnn_2d = neighbors(sites_2d; num_distance=2, flattice=flattice_2d_diag)

            # Each site has 4 diagonal neighbors on a square lattice => 16 * 4 / 2 = 32
            @test length(nnn_2d) == 32

            for (i, j) in nnn_2d
                @test metric_2d(sites_2d[i], sites_2d[j]) ≈ sqrt(2.0) atol=1e-10
            end
        end

        

        # -------------------------------------------------------
        # 2D square lattice with sheared boundary [4 1; 0 3]
        # -------------------------------------------------------
        @testset "2D square sheared boundary NN" begin
            # 12 sites in the sheared cell (det([4 1; 0 3]) = 12)
            # Generate sites via bravais_cells or manually
            sites_sheared = EuclideanVector[]
            for ix in 0:5, iy in 0:5
                trial = [Float64(ix), Float64(iy)]
                # Express in boundary basis: [4 1; 0 3]^{-T} * trial
                # Check if fractional coords are in [0, 1)
                frac = inv(Matrix{Float64}(boundary_2d)') * trial
                if all(frac .>= -1e-10) && all(frac .< 1.0 - 1e-10)
                    push!(sites_sheared, EuclideanVector(trial))
                end
            end
            @test length(sites_sheared) == 12  # det = 12

            dists_sheared = distances(sites_sheared; flattice=flattice_2d)
            @test dists_sheared[1] ≈ 0.0

            # NN distance is still 1.0 for a square lattice regardless of boundary shape
            @test dists_sheared[2] ≈ 1.0

            nn_sheared = neighbors(sites_sheared; num_distance=1, flattice=flattice_2d)

            # Each site has 4 NNs on square lattice => 12 * 4 / 2 = 24 pairs
            @test length(nn_sheared) == 24

            metric_sheared = PeriodicEuclideanMetric(flattice_2d)
            for (i, j) in nn_sheared
                @test metric_sheared(sites_sheared[i], sites_sheared[j]) ≈ 1.0
                @test i < j
            end
        end

        

        # -------------------------------------------------------
        # 3D cubic lattice, 3×3×3 diagonal boundary, fully periodic
        # -------------------------------------------------------
        @testset "3D cubic 3×3×3 periodic NN" begin
            sites_3d = EuclideanVector[]
            for ix in 0:2, iy in 0:2, iz in 0:2
                push!(sites_3d, EuclideanVector([Float64(ix), Float64(iy), Float64(iz)]))
            end
            @test length(sites_3d) == 27

            dists_3d = distances(sites_3d; flattice=flattice_3d_diag)
            @test dists_3d[1] ≈ 0.0
            @test dists_3d[2] ≈ 1.0        # NN
            @test dists_3d[3] ≈ sqrt(2.0)   # 2nd NN (face diagonal)
            @test dists_3d[4] ≈ sqrt(3.0)   # 3rd NN (body diagonal)

            nn_3d = neighbors(sites_3d; num_distance=1, flattice=flattice_3d_diag)

            # Each of 27 sites has 6 NNs => 27 * 6 / 2 = 81 pairs
            @test length(nn_3d) == 81

            metric_3d = PeriodicEuclideanMetric(flattice_3d_diag)
            for (i, j) in nn_3d
                @test metric_3d(sites_3d[i], sites_3d[j]) ≈ 1.0
                @test i < j
            end

            # 2nd NN (face-diagonal)
            nnn_3d = neighbors(sites_3d; num_distance=2, flattice=flattice_3d_diag)
            # Each site has 12 face-diagonal neighbors => 27 * 12 / 2 = 162
            @test length(nnn_3d) == 162

            for (i, j) in nnn_3d
                @test metric_3d(sites_3d[i], sites_3d[j]) ≈ sqrt(2.0) atol=1e-10
            end
        end


        # -------------------------------------------------------
        # 3D cubic lattice with sheared boundary [3 1 0; 0 3 1; 0 0 3]
        # -------------------------------------------------------
        @testset "3D cubic sheared boundary NN" begin
            # det([3 1 0; 0 3 1; 0 0 3]) = 27 sites
            sites_3d_sheared = EuclideanVector[]
            bmat = Matrix{Float64}(boundary_3d)
            for ix in -3:5, iy in -3:5, iz in -3:5
                trial = [Float64(ix), Float64(iy), Float64(iz)]
                frac = inv(bmat') * trial
                if all(frac .>= -1e-10) && all(frac .< 1.0 - 1e-10)
                    push!(sites_3d_sheared, EuclideanVector(trial))
                end
            end
            @test length(sites_3d_sheared) == 27  # det = 27

            dists_3d_sh = distances(sites_3d_sheared; flattice=flattice_3d)
            @test dists_3d_sh[1] ≈ 0.0
            @test dists_3d_sh[2] ≈ 1.0  # NN for cubic lattice

            nn_3d_sh = neighbors(sites_3d_sheared; num_distance=1, flattice=flattice_3d)

            # Still 27 sites with 6 NNs each => 81 pairs
            @test length(nn_3d_sh) == 81

            metric_3d_sh = PeriodicEuclideanMetric(flattice_3d)
            for (i, j) in nn_3d_sh
                @test metric_3d_sh(sites_3d_sheared[i], sites_3d_sheared[j]) ≈ 1.0
                @test i < j
            end
        end

        

        # -------------------------------------------------------
        # Edge cases for neighbors()
        # -------------------------------------------------------
        @testset "neighbors edge cases" begin
            sites_small = [
                EuclideanVector([0.0, 0.0]),
                EuclideanVector([1.0, 0.0]),
            ]

            # num_distance < 0 should error
            @test_throws ErrorException neighbors(sites_small; num_distance=-1, flattice=flattice_2d_diag)

            # num_distance too large should error
            @test_throws ErrorException neighbors(sites_small; num_distance=100, flattice=flattice_2d_diag)
        end
    end

    

    # ============================================================
    # distance_matrix() tests
    # ============================================================
    @testset "distance_matrix()" begin

        @testset "Euclidean 2D" begin
            pts = [
                EuclideanVector([0.0, 0.0]),
                EuclideanVector([1.0, 0.0]),
                EuclideanVector([0.0, 1.0]),
            ]
            dm = distance_matrix(pts)

            # Diagonal should be zero
            for i in 1:3
                @test dm[i, i] ≈ 0.0
            end

            # Symmetry
            for i in 1:3, j in 1:3
                @test dm[i, j] ≈ dm[j, i]
            end

            @test dm[1, 2] ≈ 1.0
            @test dm[1, 3] ≈ 1.0
            @test dm[2, 3] ≈ sqrt(2.0)
        end

        @testset "Periodic 2D" begin
            pts = [
                EuclideanVector([0.0, 0.0]),
                EuclideanVector([3.5, 0.0]),
                EuclideanVector([0.0, 3.5]),
            ]
            dm_per = distance_matrix(pts; flattice=flattice_2d_diag)

            # Wrapping distances
            @test dm_per[1, 2] ≈ 0.5 atol=1e-10
            @test dm_per[1, 3] ≈ 0.5 atol=1e-10

            # Symmetry
            for i in 1:3, j in 1:3
                @test dm_per[i, j] ≈ dm_per[j, i]
            end
        end

        @testset "Periodic 3D" begin
            pts_3d = [
                EuclideanVector([0.0, 0.0, 0.0]),
                EuclideanVector([1.0, 0.0, 0.0]),
                EuclideanVector([2.8, 0.0, 0.0]),
            ]
            dm_3d = distance_matrix(pts_3d; flattice=flattice_3d_diag)
            @test dm_3d[1, 1] ≈ 0.0
            @test dm_3d[1, 2] ≈ 1.0
            @test dm_3d[1, 3] ≈ 0.2 atol=1e-10
        end
    end

    # ============================================================
    # distances() uniqueness tests
    # ============================================================
    @testset "distances()" begin
        # On a 4×4 periodic square lattice
        sites_2d = EuclideanVector[]
        for ix in 0:3, iy in 0:3
            push!(sites_2d, EuclideanVector([Float64(ix), Float64(iy)]))
        end

        dists = distances(sites_2d; flattice=flattice_2d_diag)

        # Should be sorted
        @test issorted(dists)

        # No duplicates (unique)
        for i in 2:length(dists)
            @test dists[i] - dists[i-1] > 1e-11
        end

        # First entry is self-distance 0
        @test dists[1] ≈ 0.0

        # Euclidean distances (no periodicity) — should include larger values
        dists_euc = distances(sites_2d)
        @test maximum(dists_euc) >= maximum(dists)  # periodic wrapping reduces max distance
    end
    

end
