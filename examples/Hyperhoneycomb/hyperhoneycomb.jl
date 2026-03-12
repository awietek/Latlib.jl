using Revise
using Latlib

# Hyperhoneycomb: 3D lattice with 4 atoms per unit cell
fl = FiniteLattice(hyperhoneycomb, [2 0 0; 0 2 0; 0 0 2])

plot_3d(fl; show_neighbors=true, show_unit_cell=true)