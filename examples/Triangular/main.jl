using Latlib

function order_a2first_triangular(A::Vector{Float64},B::Vector{Float64})

    BravaisVec = [1 0.5; 0.0 sqrt(3)/2]
    P1 = round.(BravaisVec \ A)
    P2 = round.(BravaisVec \ B)
    return (P1[1] < P2[1]) || ((P1[1] == P2[1]) && (P1[2] < P2[2]))
end


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
    FiniteLat = Latlib.FiniteLattice(InFlattice,SimTorus,[false, true],order_yfirst)


    H = OpSum()
    
    H += neighbor_bonds("HB", "J1", FiniteLat; num_distance = 1)

    Latlib.GLMakie.activate!(inline=true)
    ## Print Lattice
    plot_opsum(H,FiniteLat)

    write_opsum_to_toml!(H,"triangular_J1.toml",index_zero=false)

end