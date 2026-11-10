---
tags: [capitulo2, definicion, algoritmo]
---

# Matrices triangulares: sustitución hacia adelante y hacia atrás

> [!idea]
> Si $A$ ya es triangular, resolver $Ax=b$ cuesta solo $O(n²)$ (despejar una incógnita por vez), contra $O(n³)$ de una matriz llena. Por eso **todo** el capítulo apunta a reducir $A$ a forma triangular una sola vez.

## Definiciones

- **Triangular inferior**: $L ∈ ℝⁿˣⁿ$ con $ℓᵢⱼ=0$ si $i<j$.
- **Triangular inferior unitaria**: además $ℓᵢᵢ=1$ para todo $i$.
- **Triangular superior**: $U ∈ ℝⁿˣⁿ$ con $uᵢⱼ=0$ si $i>j$.

$\det(L)=ℓ₁₁ℓ₂₂⋯ℓₙₙ$, así que $Lx=b$ tiene solución única sii todos los $ℓᵢᵢ≠0$.

## Sustitución hacia adelante (forward substitution)

Para $Lx=b$ con $L$ triangular inferior:
$$
x_1 = b_1/\ell_{11}, \qquad x_i = \Big(b_i-\sum_{j=1}^{i-1}\ell_{ij}x_j\Big)\Big/\ell_{ii}, \quad i=2,\dots,n.
$$

```julia
"""
    triinf(L, b)

Resuelve Lx = b suponiendo L triangular inferior cuadrada,
con todos los elementos de la diagonal no nulos.
Un solo lazo for (acceso por filas, vectorizado con `dot`).
"""
function triinf(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = length(b)
    x = copy(b)
    @inbounds for i in 1:n
        x[i] -= dot(@view(L[i, 1:i-1]), @view(x[1:i-1]))
        x[i] /= L[i, i]
    end
    return x
end
```

## Sustitución hacia atrás (backward substitution)

Para $Ux=b$ con $U$ triangular superior, análogo pero empezando por $xₙ$:
$$
x_n = b_n/u_{nn}, \qquad x_i = \Big(b_i-\sum_{j=i+1}^{n}u_{ij}x_j\Big)\Big/u_{ii}, \quad i=n-1,\dots,1.
$$

```julia
"""
    trisup(U, b)

Resuelve Ux = b suponiendo U triangular superior cuadrada,
con todos los elementos de la diagonal no nulos.
"""
function trisup(U::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = length(b)
    x = copy(b)
    @inbounds for i in n:-1:1
        x[i] -= dot(@view(U[i, i+1:n]), @view(x[i+1:n]))
        x[i] /= U[i, i]
    end
    return x
end
```

## Conteo de operaciones

Contando por columnas (equivalente a por filas, solo cambia el patrón de acceso a memoria):

| Momento | mult. | sumas | div. |
|---|---|---|---|
| iteración $k$ (quedan $n-k$ filas por actualizar) | $n-k$ | $n-k$ | 1 |
| **total** | $(n-1)n/2$ | $(n-1)n/2$ | $n$ |

Total: $(n-1)n + n = n²$ operaciones. **Todas** las variantes (acceso por filas, por columnas, con $@simd$/$@views$) hacen exactamente el mismo número de operaciones: lo que cambia es *cómo* se accede a la memoria (ver [[Costo computacional metodos directos]]).

> [!pendiente]
> En la práctica, ninguna de estas implementaciones manuales le gana a `L \ b` de Julia: `\` despacha a rutinas de LAPACK/BLAS (multi-hilo, bloqueo por caché). Implementarlas a mano es para *entender* el algoritmo, no para usarlas en producción.

## Ver también
- [[Metodo de Gauss - eliminacion]]
- [[Factorizacion LU]]
- [[Resolver sistemas lineales en Julia]]
