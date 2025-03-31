using Latlib
let 
    ### Lattice Dimensions    
    L = 6
    W = 4
    
    #### Matrix with unit vectors in each collum;
    BravaisVec = [1 0; 0 1]
    SimTorus =   [L÷2 0; 0 W÷2]

    Basis = [0 0; 0.0 0.5; 0.5 0.0; 0.5 0.5]' # 1 Atom per unit cell!

    ## Infite Bravais Lattice:
    InFlattice = Latlib.Lattice(BravaisVec, Basis)
    
    ## Finite Bravais Lattice inside of Simulation Torus:
    FiniteLat = Latlib.FiniteLattice(InFlattice,SimTorus,[false, true])

    H = OpSum()
    
    H += neighbor_bonds("HB", "J", FiniteLat; num_distance = 1)
    H += lattice_bonds("HB", "Jd", FiniteLat, [1, 0, 0], [4, 0, 0])
    H += lattice_bonds("HB", "Jd", FiniteLat, [3, 0, 1], [2, 1, 0])
    
    Latlib.GLMakie.activate!(inline=true)
    ## Print Lattice

    plot_opsum(H,FiniteLat)

    write_opsum_to_toml!(H,"shastry_sutherland.toml",index_zero=false)
end