---
tags: [capitulo2, definicion, teorema, algoritmo]
---

# Factorización de Cholesky. Matrices sdp

> [!idea]
> Si $A$ es sdp, existe $R$ triangular superior con $A=RᵀR$ (en vez de $A=LU$ con dos factores distintos). Esto cuesta la **mitad** de Gauss/[[Factorizacion LU|LU]] ($≅n³/3$ en vez de $2n³/3$) porque se explota la simetría, y es automáticamente estable (no hace falta pivotear).

## Definición: sdp

$A∈ℝⁿˣⁿ$ es **simétrica y definida positiva (sdp)** si $Aᵀ=A$ y $xᵀAx>0$ para todo $x≠0$. Es **semidefinida positiva** si alcanza con $xᵀAx≥0$.

## Propiedades de las matrices sdp (muy preguntables)

- $A$ sdp ⟹ $A$ invertible ($Ax=0 ⟹ xᵀAx=0 ⟹ x=0$).
- $A$ sdp ⟹ diagonalizable ortogonalmente con todos los autovalores reales **positivos**: $A=UΛUᵀ$.
- $RᵀR$ es siempre simétrica semidefinida positiva, para cualquier $R$; es sdp si las columnas de $R$ son linealmente independientes (ver [[2.13 - Propiedades de matrices sdp|Problema 2.13]]).
- Si $A$ es sdp y se escribe por bloques $A=\begin{bmatrix}α&aᵀ\\a&A_*\end{bmatrix}$, entonces $α>0$ y $A_*$ es sdp (hereditario).
- Todos los elementos de la diagonal de una matriz sdp son positivos: $a_{ii}=e_i^TAe_i>0$.
- Toda matriz se escribe como parte simétrica + antisimétrica; si $xᵀAx>0$ para todo $x≠0$, entonces la **parte simétrica** $(A+Aᵀ)/2$ es sdp (la parte antisimétrica no aporta a $xᵀAx$, porque $xᵀAx$ es un escalar y su "transpuesta" es él mismo).

## Teorema de existencia

Si $A$ es sdp, existe $R$ triangular superior con $A=RᵀR$. Demostración completa (constructiva, por inducción en $n$): [[Teorema de existencia de Cholesky]].

**Observación:** la factorización con diagonal positiva es **única**. Además, $A=RᵀR = (RᵀD⁻¹)(DR) = LU$ con $D=\text{diag}(R)$: Cholesky es un caso particular de LU (con $L=Rᵀ$, salvo el reescalado por $D$).

## Algoritmo

```julia
"""
    mi_cholesky!(A)

Calcula R triangular superior tal que A = RᵀR, sobreescribiendo A.
Explota la simetría: solo actualiza el triángulo superior.
"""
function mi_cholesky!(A::AbstractMatrix)
    n = size(A, 1)
    for c in 1:n-1
        ρ = sqrt(A[c, c])
        A[c, c] = ρ
        A[c, c+1:n] ./= ρ                      # fila c de R
        A[c+1:n, c] .= 0
        @views A[c+1:n, c+1:n] .-= A[c, c+1:n] * A[c, c+1:n]'   # actualización rango 1
    end
    A[n, n] = sqrt(A[n, n])
    return UpperTriangular(A)
end
```

En Julia, `cholesky(A)` (de $LinearAlgebra$) ya está optimizada y es la que conviene usar en la práctica.

## Costo

Explotando la simetría, Cholesky cuesta $≅n³/3$ operaciones: **la mitad** de Gauss/LU ($2n³/3$), y sin necesidad de pivotear (siempre es estable para matrices sdp). Ver [[Costo computacional metodos directos]].

## Ver también
- [[Factorizacion LU]]
- [[2.13 - Propiedades de matrices sdp|Problema 2.13]]
- [[2.14 - Cholesky recursiva|Problema 2.14]]
- En el Capítulo 3: las ecuaciones normales $AᵀAx=Aᵀb$ se resuelven con Cholesky (ver [[Minimos cuadrados y ecuaciones normales]]).
- En el Capítulo 5: [[Precondicionadores genericos|Cholesky incompleta]] como precondicionador.
