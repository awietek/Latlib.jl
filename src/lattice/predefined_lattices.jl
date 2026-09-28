# ----------------------------------------------------
#                     2D LATTICES
# ----------------------------------------------------


# ----- square lattice -----
A_square = [1.0 0.0;
            0.0 1.0]
pos_square = [0.0 0.0]
"""
    square

Square lattice with lattice vectors ``\\mathbf{a}_1 = (1, 0)``, ``\\mathbf{a}_2 = (0, 1)``
and one atom per unit cell.
"""
const square = Lattice(A_square, pos_square)

# ----- triangular lattice -----
theta = pi/3
a1 = [cos(theta) +sin(theta)]
a2 = [cos(theta) -sin(theta)]
A_tri = Matrix(vcat(a1, a2))
pos_tri = [0.0 0.0]
"""
    triangular

Triangular lattice with lattice vectors ``\\mathbf{a}_{1,2} = (\\cos\\frac{\\pi}{3}, \\pm\\sin\\frac{\\pi}{3})``
and one atom per unit cell.
"""
const triangular = Lattice(A_tri, pos_tri)

# ----- Honeycomb lattice -----
theta = pi/6
a1 = [cos(theta) +sin(theta)]
a2 = [cos(theta) -sin(theta)]
A_honeycomb = Matrix(vcat(a1, a2))
pos_honeycomb = [0.0 0.0;
                 1/3 1/3]
"""
    honeycomb

Honeycomb lattice with lattice vectors ``\\mathbf{a}_{1,2} = (\\cos\\frac{\\pi}{6}, \\pm\\sin\\frac{\\pi}{6})``
and two atoms per unit cell at ``(0, 0)`` and ``(1/3, 1/3)`` (lattice basis).
"""
const honeycomb = Lattice(A_honeycomb, pos_honeycomb)

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
"""
    shastry_sutherland

Shastry-Sutherland lattice: square Bravais lattice with four atoms per unit cell at
``(0, 0)``, ``(0, 1/2)``, ``(1/2, 0)``, ``(1/2, 1/2)`` (lattice basis). The dimer bonds
have to be added explicitly with [`lattice_interaction`](@ref), see `examples/Shastry_Sutherland`.
"""
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
"""
    kagome

Kagome lattice with lattice vectors ``\\mathbf{a}_1 = (1, 0)``, ``\\mathbf{a}_2 = (1/2, \\sqrt{3}/2)``
and three atoms per unit cell at ``(0, 0)``, ``(1/2, 0)``, ``(0, 1/2)`` (lattice basis).
"""
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
"""
    hyperhoneycomb

Three-dimensional hyperhoneycomb lattice with lattice vectors
``\\mathbf{a}_1 = (2, 4, 0)``, ``\\mathbf{a}_2 = (3, 3, 2)``, ``\\mathbf{a}_3 = (-1, 1, 2)``
and four atoms per unit cell at ``(0,0,0)``, ``(1,1,0)``, ``(1,2,1)``, ``(2,3,1)`` (Cartesian coordinates).
"""
const hyperhoneycomb = Lattice(A_hyp, EuclideanVector.(pos_hyp_eucl_coords))

