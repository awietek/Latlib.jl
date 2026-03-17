using LinearAlgebra

function order_xfirst(A::Vector{Float64},B::Vector{Float64})
    return (A[1] < B[1]) || ((A[1] == B[1]) && (A[2] < B[2]))
end

# sort 
function order_yfirst(A::Vector{Float64}, B::Vector{Float64})
    return (A[2] < B[2]) || ((A[2] == B[2]) && (A[1] < B[1]))
end


