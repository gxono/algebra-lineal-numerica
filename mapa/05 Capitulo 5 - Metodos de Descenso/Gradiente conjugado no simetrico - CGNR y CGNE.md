---
tags: [capitulo5, algoritmo]
---

# Gradiente conjugado para problemas no simétricos: CGNR y CGNE

> [!idea]
> Truco simple para usar GC (que necesita $A$ sdp) en matrices no simétricas: aplicarlo a las *ecuaciones normales* ($A^TA$ o $AA^T$, siempre sdp si $A$ es invertible). El costo: $\kappa_2(A^TA)=\kappa_2(A)^2$ (igual que en [[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]], Capítulo 3) - puede converger mucho más lento. **[[GMRes - definicion y propiedad de minimizacion|GMRes]] suele ser preferible** si $A$ está mal condicionada.

## CGNR (minimiza el Residuo)

$Ax=b \iff A^TAx=A^Tb$ ($A^TA$ sdp). GC sobre este sistema minimiza $\|x-x_*\|_{A^TA}^2=\|A(x-x_*)\|_2^2=\|Ax-b\|_2^2$ - minimiza el **residuo** $\|Ax-b\|_2$. No hace falta formar $A^TA$: alcanza con la rutina $ξ\to A^T(A\xi)$.

## CGNE (minimiza el Error)

$Ax=b$, con $x=A^Ty$: $AA^Ty=b$ ($AA^T$ sdp). GC sobre esto minimiza $\|y-y_*\|_{AA^T}^2=\|A^T(y-y_*)\|_2^2=\|x-x_*\|_2^2$ - minimiza el **error**.

## Algoritmo

Ambos reutilizan `gradconj` ([[Implementacion del gradiente conjugado]]) sobre la acción matrix-free correspondiente - nunca se forma $A^TA$ ni $AA^T$ explícitamente.

```julia
function cgnr(A, b, x0; tol=1e-8, max_iter=1000)
    AtA(v) = A'*(A*v)              # acción de A^TA, sin formarla
    return gradconj(AtA, A'b, x0; tol=tol, max_iter=max_iter)
end

function cgne(A, b, y0; tol=1e-8, max_iter=1000)
    AAt(v) = A*(A'v)               # acción de AA^T, sin formarla
    y, iter = gradconj(AAt, b, y0; tol=tol, max_iter=max_iter)
    return A'y, iter                # x = A^T y
end
```

## El problema

$$
\kappa_2(A^TA)=\kappa_2(A)^2.
$$
Si $κ₂(A)$ ya es grande, elevarlo al cuadrado puede hacer que CGNR/CGNE necesiten muchas más iteraciones que GMRes sobre el sistema original.

## Ver también
- [[Minimos cuadrados y ecuaciones normales]] (Capítulo 3, el mismo fenómeno)
- [[GMRes - definicion y propiedad de minimizacion]]
- [[5.18 - Problema con conveccion|Problema 5.18]]
