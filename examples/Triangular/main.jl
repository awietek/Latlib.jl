using Latlib
using Printf
using GLMakie

function order_xfirst(A::Vector{Float64},B::Vector{Float64})
    return (A[2] < B[2]) || ((A[2] == B[2]) && (A[1] < B[1]))
end


let 
    ### Lattice Dimensions    
    L = 12
    W = 4
    
    #### Matrix with unit vectors in each collum;
    BravaisVec = [1 0.5; 0 sqrt(3)/2]
    SimTorus =   [L 0; 0 W]
    Basis = [0 0]' # 1 Atom per unit cell!
    
    ## Infite Bravais Lattice:
    InFlattice = Latlib.Lattice(BravaisVec, Basis)
    

    ## Finite Bravais Lattice inside of Simulation Torus:
   

    ## Get coordinates
    #op1 = Latlib.Op("HB","Jd",[1 0])

    ## Print Lattice
 
    
    
end


# Assuming OpSum is a struct with a field `ops` that holds the vector of elements


### Lattice Dimensions    
L = 8
W = 4

#### Matrix with unit vectors in each collum;
BravaisVec = [1 0.5; 0 sqrt(3)/2]
SimTorus =   [L 0; 0 W]
Basis = [0 0]' # 1 Atom per unit cell!

## Infite Bravais Lattice:
InFlattice = Latlib.Lattice(BravaisVec, Basis)


FiniteLat = Latlib.FiniteLattice(InFlattice,SimTorus,[false, true])

j1_bonds1 = neighbor_bonds("HB", "Jd", FiniteLat; num_distance = 2)

# Correct usage of neighbors with num_distance keyword
neighbor_list = neighbors(FiniteLat; num_distance = 1)

GLMakie.activate!(inline=true)
Latlib.plot(FiniteLat)