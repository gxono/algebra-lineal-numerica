using LinearAlgebra, Symbolics


# ACTIVIDAD 1.a

M = diagm(0 => [0.8, -0.8, 0.8], 1 => [1.0, 1.0])

@info opnorm(M, 1)
@info opnorm(M, 2)
@info opnorm(M, Inf)




# ACTIVIDAD 1.b 
x0 = [1.0,1.0,1.0]
ρₘ = maximum(abs.(eigvals(M)))




# ACTIVIDAD 2 
@variables ε
function A(v)
  A = diagm(4, 3, -1 => fill(v, 3))
  A[1,:] .= 1.0
  return A
end

M = A(ε)
b = [1.0,0.0,0.0,0.0]


cond(A(10^(-4)))