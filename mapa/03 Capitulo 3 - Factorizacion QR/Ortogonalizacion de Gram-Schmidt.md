---
tags: [capitulo3, teorema, algoritmo]
---

# Ortogonalización de Gram-Schmidt

> [!idea]
> Construcción directa de $A=Q₁R₁$: ir ortonormalizando las columnas de $A$ una por una, restando a cada nueva columna sus componentes en la dirección de las ya procesadas.

## El procedimiento

Dados $v_1,\dots,v_m$ li:
$$
u_1=\frac{v_1}{\|v_1\|}, \qquad w_k = v_k-\sum_{i=1}^{k-1}(v_k\cdot u_i)u_i, \qquad u_k=\frac{w_k}{\|w_k\|}, \quad k=2,\dots,m.
$$
Cada $u_k$ es combinación lineal de $v_1,\dots,v_k$ con coeficiente de $v_k$ no nulo, así que
$$
\begin{bmatrix}v_1&\cdots&v_m\end{bmatrix} = \begin{bmatrix}u_1&\cdots&u_m\end{bmatrix}R_1, \qquad R_1 \text{ triangular superior}.
$$

## Teorema 3.5

Si las columnas de $A\in\mathbb{R}^{n\times m}$ son li, existe $Q_1\in\mathbb{R}^{n\times m}$ ortogonal y $R_1\in\mathbb{R}^{m\times m}$ triangular superior con $A=Q_1R_1$. *(Completando $Q_1$ a una base ortonormal de $\mathbb{R}^n$ se obtiene la versión cuadrada $A=QR$ de la Observación 3.6.)* Demostración completa (constructiva, por inducción): [[Teorema 3.5]].

## Algoritmo

```julia
"""
    gram_schmidt!(M::AbstractMatrix)

Ortonormaliza las columnas de M in-place (versión modificada: resta
proyecciones ya actualizadas, más estable que la versión clásica).
"""
function gram_schmidt!(M::AbstractMatrix)
    n, m = size(M)
    M[:, 1] ./= norm(@view M[:, 1])
    @views for j in 2:m
        for k in 1:j-1
            M[:, j] .-= dot(M[:, j], M[:, k]) .* M[:, k]
        end
        normalize!(M[:, j])
    end
    return M
end
```

> [!pendiente]
> Esta versión ("modificada") resta las proyecciones de a una, recalculando con la columna ya parcialmente corregida - eso la hace más estable numéricamente que el Gram-Schmidt "clásico" (que proyecta usando la columna original sin corregir). Comparar con el fenómeno análogo en [[Rotaciones de Givens]] / el remedio de Gram-Schmidt modificado para [[GMRes - definicion y propiedad de minimizacion|GMRes]] (Capítulo 5).

## Clásico vs. modificado: evidencia numérica

La diferencia entre ambas versiones **no se nota** con columnas bien condicionadas - hace falta un caso extremo para que se vea. La [[Parcial Modelo - Problema 2|matriz de Läuchli]] del parcial modelo es ideal para esto: Gram-Schmidt clásico (proyectando siempre contra la columna *original*, como en la implementación de abajo) se rompe mucho antes que las ecuaciones normales, y lo hace **sin avisar** - no lanza una excepción como Cholesky, directamente devuelve una respuesta incorrecta.

```julia
# clásico: todas las proyecciones usan la columna original A[:,j] (no la ya corregida)
function gram_schmidt_clasico(A)
    n, m = size(A)
    Q = zeros(n, m)
    R = zeros(m, m)
    for j in 1:m
        v = A[:, j]
        for i in 1:j-1
            R[i, j] = Q[:, i]' * A[:, j]
            v -= R[i, j] * Q[:, i]
        end
        R[j, j] = norm(v)
        Q[:, j] = v / R[j, j]
    end
    return Q, R
end
```

Para la matriz de Läuchli ($\varepsilon$ de $10^{-2}$ a $10^{-9}$): el error relativo del clásico ya es notable en $\varepsilon=10^{-6}$ ($7\times10^{-5}$, mientras Householder sigue en $10^{-16}$), y para $\varepsilon\leq10^{-8}$ da una respuesta con error relativo $\approx\sqrt2$ - completamente incorrecta. Tabla completa y comparación con los otros dos métodos: [[Parcial Modelo - Problema 2]].

## Costo

Ortonormalizar las $m$ columnas cuesta $≅2nm²$ operaciones (ver [[Estabilidad y numero de condicion]] para la tabla comparativa con [[Matrices de Householder|Householder]] y ecuaciones normales).

## Ver también
- [[Matrices ortogonales]]
- [[Matrices de Householder]] (alternativa más estable para obtener QR)
- [[3.1 - Gram-Schmidt propiedades|Problema 3.1]]
- [[Parcial Modelo - Problema 2]] (clásico vs. modificado, caso extremo)
