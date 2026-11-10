---
tags: [capitulo4, algoritmo, teorema]
---

# Método de Gauss-Seidel

> [!idea]
> Igual que [[Metodo de Jacobi|Jacobi]], pero usando los valores **ya actualizados** en esta misma iteración apenas están disponibles (en vez de esperar a la iteración siguiente). Converge más rápido que Jacobi, a costa de perder el paralelismo (cada $x_i^{k+1}$ depende de $x_1^{k+1},\dots,x_{i-1}^{k+1}$).

## Construcción

$$
x_i^{k+1} = \frac1{a_{ii}}\Big(b_i-\sum_{j<i}a_{ij}x_j^{k+1}-\sum_{j>i}a_{ij}x_j^k\Big) \iff (L+D)x_{k+1} = -Ux_k+b \iff x_{k+1}=\underbrace{-(L+D)^{-1}U}_{M_{GS}}x_k+\underbrace{(L+D)^{-1}b}_{c}.
$$

## Convergencia: e.d.d. es suficiente

Si $A$ es e.d.d., $ρ(M_{GS})<1$ y por lo tanto converge. Demostración completa: [[Teorema de convergencia de Gauss-Seidel con e.d.d.]]. Si además la razón de dominancia $\alpha$ es uniforme entre filas, se puede probar algo más fuerte y cuantitativo: $\|M_{GS}\|_\infty\leq\alpha$ (la misma cota que Jacobi) - demostración completa: [[Cota de norma infinito de Gauss-Seidel para matrices e.d.d.]].

> [!pendiente]
> La demostración **clásica** (la del primer link arriba, $\rho(M_{GS})<1$) es más sutil que la de Jacobi: no sale de acotar una norma de $M_{GS}$ directamente, hace falta trabajar con el autovector y la fila $k$ donde se alcanza el máximo. Pero **sí existe** una cota directa de $\|M_{GS}\|_\infty$ (segundo link arriba, con un argumento distinto - inducción fila por fila) que da la misma cota $\alpha$ que Jacobi. Vale la pena repasar ambas antes del parcial.

## Algoritmo

```julia
"""
    gauss_seidel(A, b; tol=1e-6, maxiter=10_000)

Método de Gauss-Seidel (forward sweep).
"""
function gauss_seidel(A, b; tol=1e-6, maxiter=10_000)
    n = length(b)
    L = tril(A, -1)
    D = Diagonal(A)
    U = triu(A, 1)
    LDinv = inv(L + D)      # en la práctica: resolver el triangular, no invertir
    Mgs = -LDinv * U
    x = zeros(n)
    err = Float64[]
    for k in 1:maxiter
        xn = Mgs * x + LDinv * b
        push!(err, norm(xn - x, 1))
        x = xn
        err[end] < tol && break
    end
    return x, err
end
```

> [!pendiente]
> **Dos errores típicos al implementar esto** (de la experiencia real resolviendo el capítulo):
> 1. `tril(A)`/`triu(A)` en Julia **incluyen** la diagonal - hace falta `tril(A,-1)` y `triu(A,1)` para la parte estrictamente triangular.
> 2. Si se arma $Mgs = -LDinv*U$ (ya con el signo puesto), la iteración es $xn = Mgs*x + LDinv*b$, **sin** volver a negar $Mgs*x$ dentro del lazo - negar dos veces invierte el signo de la convergencia.

## Ver también
- [[Metodo de Jacobi]]
- [[Metodo SOR]] (generaliza Gauss-Seidel con un parámetro $ω$)
- [[4.11 - Convergencia Jacobi y Gauss-Seidel caso particular|Problema 4.11]]
