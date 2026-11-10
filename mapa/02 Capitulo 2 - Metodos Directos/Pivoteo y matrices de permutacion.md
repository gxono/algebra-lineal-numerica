---
tags: [capitulo2, teorema, algoritmo]
---

# Pivoteo y matrices de permutación

> [!idea]
> Si un pivote $a_{jj}^{(j-1)}$ da cero (o chico), Gauss falla o amplifica errores de redondeo. El pivoteo intercambia filas para usar siempre el elemento de mayor módulo disponible como pivote. El resultado es $PA=LU$ con $P$ una matriz de permutación.

## Pivoteo parcial vs. total

- **Pivoteo parcial** (el que se usa en la práctica): en el paso $j$, intercambiar la fila $j$ con la fila $k≥j$ que tiene el elemento de mayor $|·|$ en la columna $j$, *dentro* de las filas $j..n$.
- **Pivoteo total**: busca el máximo en toda la submatriz $A^{(j-1)}_{j:n,j:n}$ e intercambia también columnas. Más caro, y en la práctica no mejora significativamente sobre el parcial.

Dividir por un pivote chico amplifica errores de redondeo; el pivoteo evita eso asegurando que siempre se divide por el elemento más grande posible.

## Teorema: $PA=LU$

Con pivoteo, $PA=LU$, con $L$ triangular inferior unitaria y $|ℓᵢⱼ|≤1$ (porque el pivote elegido es siempre el de mayor módulo). Demostración completa: [[Teorema de PA=LU con pivoteo]].

## Algoritmo

```julia
"""
    eleupiv(A)

Factorización LU con pivoteo parcial. Devuelve L (unitaria inferior),
U (superior), p (vector de permutación) tales que L*U = A[p,:].
"""
function eleupiv(A::AbstractMatrix)
    n = size(A, 1)
    L, U, p = zeros(n, n), float(copy(A)), collect(1:n)
    for c in 1:n-1
        idx = argmax(abs.(U[c:n, c])) + c - 1
        p[[idx, c]] = p[[c, idx]]
        L[[idx, c], 1:c-1] = L[[c, idx], 1:c-1]
        U[[idx, c], :] = U[[c, idx], :]
        L[c+1:n, c] = U[c+1:n, c] / U[c, c]
        U[c+1:n, c] .= 0
        U[c+1:n, c+1:n] .-= L[c+1:n, c] * U[c, c+1:n]'
    end
    return LowerTriangular(L + I), UpperTriangular(U), p
end
```

Resolver $Ax=b$: $L*(U*x) = A[p,:]*x = b[p]$, así que `y = triinf(L, b[p])`, `x = trisup(U, y)`.

## Propiedades de las matrices de permutación

Sea $P$ la matriz que se obtiene intercambiando las filas $i,j$ de la identidad. $P$ es simétrica, $PP=I$, y $PA$/$AP$ intercambian filas/columnas de $A$. Demostraciones: [[Matrices de permutacion - P es simetrica|$P$ es simétrica]], [[Matrices de permutacion - PP=I|$PP=I$]], [[Matrices de permutacion - PA y AP intercambian filas y columnas|$PA$/$AP$ intercambian filas/columnas]].

## Ver también
- [[Factorizacion LU]]
- [[Metodo de Gauss - eliminacion]]
- [[2.8 - Inestabilidad sin pivoteo|Problema 2.8]] / [[2.9 - Estabilidad con pivoteo|Problema 2.9]] (el ejemplo numérico que muestra *por qué* hace falta pivotear)
