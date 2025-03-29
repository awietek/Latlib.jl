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
    FiniteLattice = Latlib.FiniteLattice(InFlattice,SimTorus,[false, true])

    ## Get coordinates
    #op1 = Latlib.Op("HB","Jd",[1 0])

    ## Print Lattice
    GLMakie.activate!(inline=true)
    Latlib.plot(FiniteLattice)
    
    
end

A = OpSum()

@show A

# Create Op objects
op1 = Op("HB", "Jd", [1, 0])
op2 = Op("HB", "Jd", [2, 1])
op3 = Op("HB", "Jd", [3, 2])

# Add Op objects to OpSum using +=

A += op2
A += op3


    N = size(coords)[2]
    
    
    jd_bonds1 = lattice_bonds("HB", "Jd", lat, [1, 0, 0], [4, 0, 0]; periodic_dims=pd)
    
    jd_bonds2 = lattice_bonds("HB", "Jd", lat, [3, 0, 1], [2, 1, 0]; periodic_dims=pd)
    
    jd_bonds = vcat([jd_bonds1, jd_bonds2]...)
    
    jd_sites = sort(unique(vcat([b.sites for b in jd_bonds]...)))
    