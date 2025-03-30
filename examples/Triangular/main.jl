using Latlib
using Printf
using GLMakie
using TOML


let 
    ### Lattice Dimensions    
    L = 12
    W = 4
    
    #### Matrix with unit vectors in each collum;
    BravaisVec = [1 0.5; 0 sqrt(3)/2]
    SimTorus = [L 0; 0 W]
    Basis = [0 0]' # 1 Atom per unit cell!
    
    ## Infite Bravais Lattice:
    InFlattice = Latlib.Lattice(BravaisVec, Basis)
    

    ## Finite Bravais Lattice inside of Simulation Torus:
    FiniteLat = Latlib.FiniteLattice(InFlattice,SimTorus,[false, true])


    H = OpSum()

    H += neighbor_bonds("HB", "J1", FiniteLat; num_distance = 1)

    GLMakie.activate!(inline=true)
    ## Print Lattice
    Latlib.plotOpSum(H,FiniteLat)

    write_OpSum_to_toml!(H,"triangular.toml",index_zero=false)

end