# ----- 2D Lattices -----

# square lattice
A_square = [1.0 0.0;
            0.0 1.0]
positions_square = [0.0 0.0]
const square = Lattice(A_square, positions_square)

# triangular lattice
theta = pi/3
a1 = [cos(theta) +sin(theta)]
a2 = [cos(theta) -sin(theta)]
A_tri = Matrix(vcat(a1, a2))
println("A_tri = ", A_tri)
positions_tri = [0.0 0.0]
const triangular = Lattice(A_tri, positions_tri)

# kagome lattice
A_kagome = [
        1.0 0.0;
        0.5 sqrt(3)/2;
    ]
positions_kagome = [
    0.0 0.0;
    0.5 0.0;
    0.0 0.5;
    ]
const kagome = Lattice(A_kagome, positions_kagome)

