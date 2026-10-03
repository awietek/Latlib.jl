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
# a1 = [cos(theta) +sin(theta)]
# a2 = [cos(theta) -sin(theta)]
a1 = [1 0]
a2 = [cos(theta) sin(theta)]
A_tri = Matrix(vcat(a1, a2))
pos_tri = [0.0 0.0]
"""
    triangular

Triangular lattice with lattice vectors ``\\mathbf{a}_1 = (1, 0)``, ``\\mathbf{a}_2 = (1/2, \\sqrt{3}/2)``
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

# ----- Shastry-Sutherland lattice with its non-symmorphic symmetry p4gm -----
# orthogonal dimers centered at (0, 0) along [11] and at (1/2, 1/2) along [1-1]
# (lattice basis). The dimer half-length δ = (√3 - 1)/4 and the lattice constant
# √(2 + √3) make dimer and inter-dimer bonds equally long (length 1).
a_ss_ns = sqrt(2 + sqrt(3))
A_ss_ns = [a_ss_ns 0.0;
           0.0     a_ss_ns]
δ_ss_ns = (sqrt(3) - 1) / 4
pos_ss_ns = [δ_ss_ns        δ_ss_ns;
             -δ_ss_ns       -δ_ss_ns;
             0.5 + δ_ss_ns  0.5 - δ_ss_ns;
             0.5 - δ_ss_ns  0.5 + δ_ss_ns]
"""
    shastry_sutherland_non_symmorphic

Shastry-Sutherland lattice in the geometry of orthogonal dimers, which carries the
non-symmorphic plane group p4gm of the Shastry-Sutherland model. Square Bravais lattice
with lattice constant ``\\sqrt{2 + \\sqrt{3}}`` and four atoms per unit cell at
``\\pm\\delta\\,(1, 1)`` and ``(1/2, 1/2) \\pm \\delta\\,(1, -1)`` (lattice basis) with
``\\delta = (\\sqrt{3} - 1)/4``. Dimer bonds and inter-dimer bonds all have length 1, so every
site has five nearest neighbors (snub-square tiling).

Unlike [`shastry_sutherland`](@ref), whose sites lie on a square grid (plane group p4mm, with
the dimers to be added as interactions), the sites alone already have the symmetry of the model.
"""
const shastry_sutherland_non_symmorphic = Lattice(A_ss_ns, pos_ss_ns)

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

# ----- Lieb lattice -----
A_lieb = [1.0 0.0;
          0.0 1.0]
pos_lieb = [0.0 0.0;
            0.5 0.0;
            0.0 0.5]
"""
    lieb

Lieb lattice: square Bravais lattice with three atoms per unit cell, one at the corner
``(0, 0)`` and two at the edge centers ``(1/2, 0)`` and ``(0, 1/2)`` (lattice basis).
Plane group p4mm.
"""
const lieb = Lattice(A_lieb, pos_lieb)

# ----- trellis lattice -----
# two-leg ladders along x (legs and rungs of length 1); neighboring ladders are shifted
# by half a leg along x and coupled by zigzag bonds, here also of length 1
A_trellis = [1.0 0.0;
             0.5 1.0 + sqrt(3)/2]
pos_trellis_eucl_coords = [[0.0, 0.0],
                           [0.0, 1.0]]
"""
    trellis

Trellis lattice: two-leg ladders along ``x`` with legs and rungs of length 1. Neighboring
ladders are shifted by half a leg and coupled through zigzag bonds, which also have
length 1. Lattice vectors ``\\mathbf{a}_1 = (1, 0)``, ``\\mathbf{a}_2 = (1/2, 1 + \\sqrt{3}/2)``
and two atoms per unit cell at ``(0, 0)`` and ``(0, 1)`` (Cartesian coordinates).
Plane group c2mm.
"""
const trellis = Lattice(A_trellis, EuclideanVector.(pos_trellis_eucl_coords))

# ----- maple-leaf lattice -----
# 1/7-depleted triangular lattice (nearest-neighbor distance 1)
A_maple_leaf = sqrt(7) * [0.5 sqrt(3)/2;
                          1.0 0.0]
pos_maple_leaf = [0.0   0.0;
                  -1/7  3/7;
                  -2/7  6/7;
                  1/7   4/7;
                  4/7   2/7;
                  2/7   1/7]
"""
    maple_leaf

Maple-leaf lattice: triangular lattice (nearest-neighbor distance 1) with one out of
seven sites removed. Lattice vectors ``\\mathbf{a}_1 = \\sqrt{7}\\,(1/2, \\sqrt{3}/2)``,
``\\mathbf{a}_2 = \\sqrt{7}\\,(1, 0)`` and six atoms per unit cell. Plane group p6 (chiral).
"""
const maple_leaf = Lattice(A_maple_leaf, pos_maple_leaf)


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

# ----- cubic lattices (conventional cubic cell of side length 1) -----
A_simple_cubic = [1.0 0.0 0.0;
                  0.0 1.0 0.0;
                  0.0 0.0 1.0]
"""
    simple_cubic

Simple cubic lattice with lattice vectors ``\\mathbf{a}_1 = (1, 0, 0)``, ``\\mathbf{a}_2 = (0, 1, 0)``,
``\\mathbf{a}_3 = (0, 0, 1)`` and one atom per unit cell. Space group Pm-3m.
"""
const simple_cubic = Lattice(A_simple_cubic)

A_bcc = [-0.5  0.5  0.5;
          0.5 -0.5  0.5;
          0.5  0.5 -0.5]
"""
    bcc

Body-centered cubic lattice (conventional cubic cell of side 1) with primitive lattice vectors
``\\mathbf{a}_1 = (-1, 1, 1)/2``, ``\\mathbf{a}_2 = (1, -1, 1)/2``, ``\\mathbf{a}_3 = (1, 1, -1)/2``
and one atom per unit cell. Space group Im-3m.
"""
const bcc = Lattice(A_bcc)

A_fcc = [0.0 0.5 0.5;
         0.5 0.0 0.5;
         0.5 0.5 0.0]
"""
    fcc

Face-centered cubic lattice (conventional cubic cell of side 1) with primitive lattice vectors
``\\mathbf{a}_1 = (0, 1, 1)/2``, ``\\mathbf{a}_2 = (1, 0, 1)/2``, ``\\mathbf{a}_3 = (1, 1, 0)/2``
and one atom per unit cell. Space group Fm-3m.
"""
const fcc = Lattice(A_fcc)

"""
    diamond

Diamond lattice: face-centered cubic Bravais lattice with the lattice vectors of [`fcc`](@ref)
and two atoms per unit cell at ``(0, 0, 0)`` and ``(1/4, 1/4, 1/4)`` (lattice basis).
Space group Fd-3m (non-symmorphic).
"""
const diamond = Lattice(A_fcc, [0.0 0.0 0.0;
                                0.25 0.25 0.25])

"""
    pyrochlore

Pyrochlore lattice of corner-sharing tetrahedra: face-centered cubic Bravais lattice with
the lattice vectors of [`fcc`](@ref) and four atoms per unit cell at ``(0, 0, 0)``,
``(1/2, 0, 0)``, ``(0, 1/2, 0)`` and ``(0, 0, 1/2)`` (lattice basis).
Space group Fd-3m (non-symmorphic).
"""
const pyrochlore = Lattice(A_fcc, [0.0 0.0 0.0;
                                   0.5 0.0 0.0;
                                   0.0 0.5 0.0;
                                   0.0 0.0 0.5])

# ----- hexagonal lattices -----
A_simple_hexagonal = [1.0  0.0         0.0;
                      -0.5 sqrt(3)/2   0.0;
                      0.0  0.0         1.0]
"""
    simple_hexagonal

Simple hexagonal lattice (stacked triangular layers) with lattice vectors
``\\mathbf{a}_1 = (1, 0, 0)``, ``\\mathbf{a}_2 = (-1/2, \\sqrt{3}/2, 0)``, ``\\mathbf{a}_3 = (0, 0, 1)``
and one atom per unit cell. Space group P6/mmm.
"""
const simple_hexagonal = Lattice(A_simple_hexagonal)

A_hcp = [1.0  0.0        0.0;
         -0.5 sqrt(3)/2  0.0;
         0.0  0.0        sqrt(8/3)]
"""
    hcp

Hexagonal close-packed lattice with lattice vectors ``\\mathbf{a}_1 = (1, 0, 0)``,
``\\mathbf{a}_2 = (-1/2, \\sqrt{3}/2, 0)``, ``\\mathbf{a}_3 = (0, 0, \\sqrt{8/3})`` (ideal ratio
``c/a``, all twelve nearest neighbors at distance 1) and two atoms per unit cell at
``(1/3, 2/3, 1/4)`` and ``(2/3, 1/3, 3/4)`` (lattice basis). Space group P6_3/mmc (non-symmorphic).
"""
const hcp = Lattice(A_hcp, [1/3 2/3 1/4;
                            2/3 1/3 3/4])

