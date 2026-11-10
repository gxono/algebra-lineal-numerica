---
tags: [capitulo3, definicion, teorema, algoritmo]
---

# Mínimos cuadrados y ecuaciones normales

> [!idea]
> Ajustar un modelo lineal a más puntos de los que hacen falta ($n>m$ incógnitas) no tiene solución exacta; se busca minimizar $‖Ax-b‖₂$. Eso equivale a pedir que el residuo $Ax-b$ sea ortogonal a todas las columnas de $A$, lo cual da el sistema cuadrado $AᵀAx=Aᵀb$ (ecuaciones normales).

## El problema

Ajustar $(X_i,Y_i)$, $i=1,\dots,n$, con una función (ej. cuadrática $y=\alpha X^2+\beta X+\gamma$) minimizando
$$
J(\alpha,\beta,\gamma) = \sum_{i=1}^n[(\alpha X_i^2+\beta X_i+\gamma)-Y_i]^2 = \|Ax-b\|_2^2,
$$
con $A$ la matriz de diseño (columnas $X_i^2, X_i, 1$), $x=(\alpha,\beta,\gamma)$, $b=(Y_i)$.

## De la minimización a un sistema lineal

$Ax$ es combinación lineal de las columnas de $A$: minimizar $\|Ax-b\|_2$ equivale a hallar el punto del espacio columna de $A$ más cercano a $b$. Con la norma euclídea eso significa que $Ax-b$ debe ser ortogonal a cada columna de $A$:
$$
\text{col}_j(A)^T(Ax-b)=0,\ j=1,\dots,m \iff A^T(Ax-b)=0 \iff \boxed{A^TAx=A^Tb}.
$$

Si las columnas de $A$ son linealmente independientes, $AᵀA$ es sdp (ver [[2.13 - Propiedades de matrices sdp|Problema 2.13]]) y el sistema tiene solución única $x=(AᵀA)^{-1}Aᵀb$; $(AᵀA)^{-1}Aᵀ$ se llama **pseudoinversa** de $A$ (`pinv` en Julia).

## Algoritmo (vía [[Factorizacion de Cholesky|Cholesky]])

```julia
"""
    mincuad(A, b)

Resuelve min‖Ax-b‖₂ vía las ecuaciones normales AᵀAx=Aᵀb, factorizando
AᵀA con Cholesky.
"""
function mincuad(A::AbstractMatrix, b::AbstractVector)
    R = mi_cholesky!(A' * A)        # A'A = R'R
    y = triinf(R', A' * b)
    return trisup(R, y)
end
```

> [!pendiente]
> Las ecuaciones normales son la opción **más barata** ($≅nm²$, la mitad que QR) pero la **menos estable**: $κ₂(AᵀA)=κ₂(A)²$, así que si $A$ está mal condicionada, se pierde el doble de cifras significativas. Ver [[Estabilidad y numero de condicion]] para la comparación completa y cuándo conviene cada vía.

## Ver también
- [[Matrices ortogonales]] (la vía alternativa: factorizar $A=QR$)
- [[Factorizacion de Cholesky]]
- [[3.2 - Funcion mincuad|Problema 3.2]], [[3.3 - Ajuste lineal precio vivienda|Problema 3.3]]
