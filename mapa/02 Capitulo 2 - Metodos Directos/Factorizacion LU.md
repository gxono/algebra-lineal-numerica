---
tags: [capitulo2, teorema, algoritmo]
---

# Factorización LU

> [!idea]
> La [[Metodo de Gauss - eliminacion|eliminación de Gauss]], vista con cuidado, produce $A=LU$ con $L$ triangular inferior unitaria (los multiplicadores) y $U$ triangular superior (lo que queda). Si hay que resolver $Ax=b$ para **muchos** $b$ con la misma $A$, factorizar una vez ($≅2n³/3$) y resolver dos triangulares por cada $b$ ($O(n²)$) es mucho más barato que repetir Gauss entero cada vez.

## De dónde sale

El paso $j$ de Gauss es multiplicar a izquierda por una **matriz de Gauss** $Gʲ = I - gⱼeⱼᵀ$, con $gⱼ=(0,…,0,m_{j+1,j},…,m_{nj})ᵀ$. Como $(Gʲ)⁻¹ = I + gⱼeⱼᵀ$ (ver [[2.6 - Matrices de Gauss|Problema 2.6]]), y los $gⱼ$ no interfieren entre sí,
$$
(G^1)^{-1}(G^2)^{-1}\cdots(G^{n-1})^{-1} = I+\sum_{j=1}^{n-1}g_je_j^T = L,
$$
la matriz triangular inferior unitaria con los multiplicadores $mᵢⱼ$ debajo de la diagonal. Y $U = G^{n-1}\cdots G^1 A$, de donde $A = LU$.

## Teorema (existencia y unicidad)

Si $\det(A_{1:k,1:k})\neq0$ para todo $k=1,\dots,n-1$, existe $L$ triangular inferior unitaria y $U$ triangular superior con $A=LU$; si además $A$ es invertible, la factorización es única. Demostración completa: [[Teorema de existencia y unicidad de LU]].

> [!pendiente]
> **No siempre existe factorización LU.** Contraejemplo clásico: $A=\begin{bmatrix}0&1\\1&0\end{bmatrix}$ no tiene factorización LU (ver [[2.7 - Matriz sin factorizacion LU|Problema 2.7]]). La hipótesis de los menores principales no nulos es necesaria, y por eso hace falta el [[Pivoteo y matrices de permutacion|pivoteo]].

## Algoritmo

```julia
"""
    eleu(M)

Factorización LU sin pivoteo. Devuelve L (triangular inferior unitaria)
y U (triangular superior) tales que M = L*U.
"""
function eleu(M::AbstractMatrix)
    n = size(M, 1)
    U = float(copy(M))
    L = Matrix{Float64}(I, n, n)
    for c in 1:n-1
        @views for r in c+1:n
            m = U[r, c] / U[c, c]
            U[r, :] .-= m .* U[c, :]
            L[r, c] = m
        end
    end
    return L, U
end
```

Resolver $Ax=b$ con $A=LU$: $Ly=b$ (sustitución hacia adelante, $n²-n$ op.), luego $Ux=y$ (hacia atrás, $n²$ op.) - ver [[Matrices triangulares - sustitucion]].

## Por qué conviene factorizar una vez

|                                              | costo                             |
| -------------------------------------------- | --------------------------------- |
| $k$ sistemas, Gauss completo cada vez        | $k·(2n³/3)$                       |
| factorizar 1 vez + $k$ pares de triangulares | $2n³/3 + k·(2n²-n) ≅ 2n³/3 + kn²$ |

Para $k$ grande, la segunda opción gana claramente: el costo por sistema adicional baja de $O(n³)$ a $O(n²)$.

## Ver también
- [[Metodo de Gauss - eliminacion]]
- [[Pivoteo y matrices de permutacion]]
- [[Factorizacion de Cholesky]] (el caso sdp, con $L=Uᵀ$)
- [[Costo computacional metodos directos]]
