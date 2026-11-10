---
tags: [capitulo4, algoritmo, teorema]
---

# Método de Jacobi

> [!idea]
> La versión más simple: despejar $x_i$ de la ecuación $i$ usando **solo** los valores viejos ($x_k$) de las demás componentes. $M_J = -D^{-1}(L+U)$.

## Construcción

Con $A=L+D+U$, la ecuación $i$ de $Ax=b$ es $\sum_j a_{ij}x_j=b_i$; despejando $x_i$:
$$
x_i^{k+1} = \frac{1}{a_{ii}}\Big(b_i-\sum_{j\neq i}a_{ij}x_j^k\Big), \qquad\text{es decir}\qquad x_{k+1}=\underbrace{-D^{-1}(L+U)}_{M_J}x_k+\underbrace{D^{-1}b}_{c}.
$$

## Convergencia: e.d.d. es suficiente

Si $A$ es e.d.d., $‖M_J‖_∞<1$ y por lo tanto converge para todo $x_0$ (ver [[Espectro y radio espectral]]). Demostración completa: [[Teorema de convergencia de Jacobi con e.d.d.]].

## Algoritmo

```julia
"""
    jacobi(A, b; tol=1e-6, maxiter=10_000)

Método de Jacobi. Devuelve la aproximación final y el historial de error
‖x_{k+1}-x_k‖ (norma 1).
"""
function jacobi(A, b; tol=1e-6, maxiter=10_000)
    n = length(b)
    Dinv = 1 ./ diag(A)
    x = zeros(n)
    err = Float64[]
    for k in 1:maxiter
        xn = Dinv .* (b - A*x + diag(A).*x)   # = D⁻¹(b - (L+U)x)
        push!(err, norm(xn - x, 1))
        x = xn
        err[end] < tol && break
    end
    return x, err
end
```

> [!pendiente]
> Cuidado al implementar: $A*x$ ya incluye el término diagonal $a_{ii}x_i$, que **no** debe estar en la suma $\sum_{j\neq i}a_{ij}x_j$. Un bug típico (y real: me pasó resolviendo este capítulo) es olvidar sacar la diagonal, lo que da una iteración distinta a la de Jacobi.

## Ver también
- [[Metodos iterativos estacionarios - marco general]]
- [[Metodo de Gauss-Seidel]] (usa los valores *nuevos* ya calculados en la misma iteración - converge más rápido)
- [[2.11 - Matrices e.d.d. son no singulares|Problema 2.11]] (ya probaba que e.d.d. ⟹ no singular; acá además da convergencia)
