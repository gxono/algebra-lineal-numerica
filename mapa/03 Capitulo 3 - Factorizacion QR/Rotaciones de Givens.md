---
tags: [capitulo3, definicion, algoritmo]
---

# Rotaciones de Givens

> [!idea]
> Para matrices **Hessenberg** (casi triangulares, con una sola subdiagonal no nula), no hace falta [[Matrices de Householder|Householder]]: una rotación de $2×2$ por columna alcanza para anular esa única entrada sobrante. Mucho más barato, y es la herramienta clave para la implementación eficiente de [[GMRes - implementacion y rotaciones de Givens|GMRes]] en el Capítulo 5.

## Definición 3.15: matriz Hessenberg superior

$h_{ij}=0$ para todo $i>j+1$ (nula por debajo de la primera subdiagonal).

## Definición 3.16: rotación de Givens

$$
G = \begin{bmatrix}c&-s\\s&c\end{bmatrix}, \qquad c^2+s^2=1,
$$
ortogonal, y rota $(c,-s)$ hacia $(1,0)$. La versión $n×n$, $G_j(c,s)$, es la identidad con ese bloque $2×2$ en las filas/columnas $j,j+1$.

## La idea: anular la subdiagonal, una entrada por vez

Para anular $h_{21}$ de $H$, se rota el vector $(h_{11},h_{21})$:
$$
c_1 = \frac{h_{11}}{\sqrt{h_{11}^2+h_{21}^2}}, \qquad s_1 = \frac{-h_{21}}{\sqrt{h_{11}^2+h_{21}^2}}.
$$
$G_1H$ deja un único no-cero en la posición $(1,1)$ de la primera columna, y **solo modifica las dos primeras filas**. Repitiendo con $G_2$ sobre $(h_{22},h_{32})$ de la matriz actualizada, y así sucesivamente, tras $p$ pasos ($p=m-1$ si $n=m$, $p=m$ si $n>m$):
$$
\underbrace{G_p\cdots G_1}_{Q^T}H = R \ \text{triangular superior}.
$$

Cada $G_k$ solo afecta dos filas: el costo por rotación es $O(1)$ por columna que toca, mucho más barato que una reflexión de Householder completa.

## Algoritmo

```julia
"""
    rotacion_de_givens(M)

Reduce M (upper Hessenberg) a triangular superior, aplicando una rotación
de Givens por columna. Devuelve la matriz triangularizada (no arma Q).
"""
function rotacion_de_givens(M::AbstractMatrix)
    A = copy(M)
    nrow, ncol = size(A)
    @views for j in 1:ncol-1
        a, b = A[j, j], A[j+1, j]
        nrm = hypot(a, b)
        c, s = a / nrm, -b / nrm
        for k in j:ncol
            aj, aj1 = A[j, k], A[j+1, k]
            A[j, k]   = c * aj - s * aj1
            A[j+1, k] = s * aj + c * aj1
        end
    end
    return A
end
```

## Ver también
- [[Matrices de Householder]] (la alternativa general, más cara para Hessenberg)
- [[GMRes - implementacion y rotaciones de Givens]] (Capítulo 5): la aplicación estrella de Givens - factorizar $QR$ de $H_k$ de forma incremental, actualizando $R_k$ a $R_{k+1}$ con **una sola** rotación nueva por iteración.
- [[3.7 - Conteo Givens en Hessenberg|Problema 3.7]]
