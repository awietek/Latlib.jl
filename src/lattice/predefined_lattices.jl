# ----------------------------------------------------
#                     2D LATTICES
# ----------------------------------------------------


# ----- square lattice -----
A_square = [1.0 0.0;
            0.0 1.0]
pos_square = [0.0 0.0]
const square = Lattice(A_square, pos_square)

# ----- triangular lattice -----
theta = pi/3
a1 = [cos(theta) +sin(theta)]
a2 = [cos(theta) -sin(theta)]
A_tri = Matrix(vcat(a1, a2))
pos_tri = [0.0 0.0]
const triangular = Lattice(A_tri, pos_tri)

# ----- Shastry-Sutherland lattice -----
# (square Bravais lattice with 4 atoms per unit cell)
A_ss = [
    1.0 0.0;
    0.0 1.0
]
pos_ss = [0.0 0.0;
                0.0 0.5;
                0.5 0.0;
                0.5 0.5]
const shastry_sutherland = Lattice(A_ss, pos_ss)

# ----- kagome lattice -----
A_kagome = [
        1.0 0.0;
        0.5 sqrt(3)/2;
    ]
pos_kagome = [
    0.0 0.0;
    0.5 0.0;
    0.0 0.5;
    ]
const kagome = Lattice(A_kagome, pos_kagome)


# ----------------------------------------------------
#                     3D LATTICES
# ----------------------------------------------------

# ----- hyperhoneycomb lattice -----
A_hyp = [
        2.0 4.0 0.0;  # a1
        3.0 3.0 2.0;  # a2
        -1.0 1.0 2.0; # a3
    ]
pos_hyp_eucl_coords = [
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 2.0, 1.0],
        [2.0, 3.0, 1.0],
    ]
# by using the EuclideanVector constructor, we can give the atom coordiates
# in standard basis and dont have to convert to lattice basis
const hyperhoneycomb = Lattice(A_hyp, EuclideanVector.(pos_hyp_eucl_coords))

