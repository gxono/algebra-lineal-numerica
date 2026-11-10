---
tags: [capitulo4, algoritmo, teorema]
---

# Método SOR (Successive Over-Relaxation)

> [!idea]
> [[Metodo de Gauss-Seidel|Gauss-Seidel]], pero "exagerando" el paso por un factor $ω$ (relajación). Con $ω=1$ es exactamente Gauss-Seidel; $ω>1$ sobre-relaja (converge más rápido si se elige bien); $ω\in(0,2)$ es la única región donde puede converger. Para $A$ sdp, converge para **todo** $ω\in(0,2)$ (teorema de Ostrowski-Reich, ver [[4.13 - Convergencia global de SOR para A sdp|Problema 4.13]]).

## Construcción

$$
x_i^{k+1} = (1-\omega)x_i^k + \frac{\omega}{a_{ii}}\Big(b_i-\sum_{j<i}a_{ij}x_j^{k+1}-\sum_{j>i}a_{ij}x_j^k\Big),
$$
es decir, $\Big(\frac{D}{\omega}+L\Big)x_{k+1} = \Big(\frac{1-\omega}{\omega}D-U\Big)x_k+b$, de donde
$$
M_{SOR} = \Big(\frac D\omega+L\Big)^{-1}\Big(\frac{1-\omega}\omega D-U\Big).
$$

## Condición necesaria: $ω∈(0,2)$

> [!idea]
> Para **cualquier** $A$, si SOR converge entonces $0<\omega<2$. Intuición: el determinante de $M_{SOR}$ es $\prod\lambda_i(M_{SOR}) = \prod(1-\omega)$ (uno por cada autovalor, viene de la estructura triangular del armado), así que $\rho(M_{SOR})\geq|1-\omega|$; si $\omega\notin(0,2)$, $|1-\omega|\geq1$ y no puede converger.

## Condición suficiente: $A$ sdp y $ω∈(0,2)$ (Ostrowski-Reich)

Demostración completa en [[4.13 - Convergencia global de SOR para A sdp|Problema 4.13]] (es, de hecho, exactamente ese problema).

## Algoritmo

```julia
"""
    SOR(A, b, ω; tol=1e-6, maxiter=10_000)
"""
function SOR(A, b, ω; tol=1e-6, maxiter=10_000)
    n = length(b)
    D = Diagonal(A)
    L = tril(A, -1)
    U = triu(A, 1)
    N = D/ω + L
    Msor = N \ ((1-ω)/ω * D - U)
    x = zeros(n)
    err = Float64[]
    for k in 1:maxiter
        xn = Msor * x + N \ b
        push!(err, norm(xn - x, 1))
        x = xn
        err[end] < tol && break
    end
    return x, err
end
```

## Elegir $ω$

No hay una fórmula universal simple (salvo para familias especiales de matrices, fuera del alcance del curso); en la práctica se prueban varios valores y se mide $\rho(M_{SOR})$ empíricamente (el cociente $\|x_{k+1}-x_k\|/\|x_k-x_{k-1}\|$ se estabiliza en $\approx\rho(M_{SOR})$) - ver [[4.15 - Jacobi Gauss-Seidel SOR en sistemas concretos|Problema 4.15]].

## Ver también
- [[Metodo de Gauss-Seidel]]
- [[Comparacion metodos iterativos estacionarios]]
