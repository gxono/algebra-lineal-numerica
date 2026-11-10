---
tags: [capitulo5, algoritmo, definicion]
---

# Precondicionadores genéricos

> [!idea]
> Se construyen mirando solo las *entradas* de $A$ (sin usar de dónde viene la matriz). Mejoran la **constante** $c_2/c_1$ pero rara vez el **orden**: si $\kappa_2(A)=O(h^{-2})$, en general sigue valiendo $c_2/c_1=O(h^{-2})$ (con constante más chica). Para lograr $O(1)$ hace falta un precondicionador que use la estructura - ver [[Multigrid en 1D - transferencia e implementacion|multigrid]].

## Definición

Una matriz $M\approx A^{-1}$, barata de aplicar, que se usa para resolver $MAx=Mb$ en vez de $Ax=b$ directamente, con el objetivo de que $\kappa_2(MA)\ll\kappa_2(A)$ y el método iterativo converja mucho más rápido.

## La receta: reciclar un método estacionario

Si $M$ invertible, $Ax=b \iff x=(I-MA)x+Mb$, que converge sii $\rho(I-MA)<1$ (Capítulo 4). La idea: elegir $M$ tal que $I-MA=B$, la matriz de iteración de un método estacionario conocido.

## Diagonal ([[Metodo de Jacobi|Jacobi]])

$M_J=D^{-1}$. Acción trivial ($n$ divisiones). Si $A$ ya tiene diagonal constante (ej. el laplaciano discreto), $M_J$ no cambia nada.

## [[Metodo de Gauss-Seidel|Gauss-Seidel]] simétrico

Combina un barrido de Gauss-Seidel hacia adelante y uno hacia atrás:
$$
M_{SGS} = (L^T+D)^{-1}D(L+D)^{-1}.
$$
Acción: resolver un triangular inferior, multiplicar por diagonal, resolver un triangular superior.

Esta fórmula no es arbitraria: $I-M_{SGS}A$ coincide exactamente con la matriz de iteración de aplicar un barrido de Gauss-Seidel hacia adelante seguido de uno hacia atrás. Demostración completa: [[Identidad del precondicionador Gauss-Seidel simetrico]].

```julia
diagonal_precon(A) = v -> v ./ diag(A)              # acción de M_J = D⁻¹

function sgs_precon(A)         # A simétrica (L+D = LowerTriangular(A), L^T+D = UpperTriangular(A))
    D = Diagonal(A)
    return v -> UpperTriangular(A) \ (D * (LowerTriangular(A) \ v))   # acción de M_SGS
end
```

## Factorización incompleta ([[Factorizacion de Cholesky|Cholesky]] incompleta)

Si $A=R^TR+E$ con $E$ "chica" (fill-in descartado): $M=(R^TR)^{-1}$. Acción: dos sustituciones triangulares (`icholesky`, [[5.11 - GCP y Cholesky incompleta|Problema 5.11]]). Costo de preparación: factorizar **una vez**, antes de iterar.

> [!idea]
> Para matrices **tridiagonales**, Cholesky no tiene fill-in: la "incompleta" coincide con la completa, y el precondicionador es $M=A^{-1}$ **exacto** - un caso degenerado (converge en 1 sola iteración, demostración completa en [[5.11 - GCP y Cholesky incompleta|Problema 5.11]](b)), no ilustrativo de lo que pasa con precondicionadores genéricos en general. El caso interesante es **2D** (el laplaciano de 5 puntos), donde sí hay fill-in real que se descarta.

## Tabla comparativa

> [!comparacion]
> | Precondicionador | $M$ | Acción $ξ→Mξ$ | Preparación |
> |---|---|---|---|
> | Ninguno | `I` | gratis | ninguna |
> | Diagonal (Jacobi) | $D^{-1}$ | $n$ divisiones | ninguna |
> | Gauss-Seidel simétrico | $(L^T+D)^{-1}D(L+D)^{-1}$ | 2 triangulares + 1 diagonal | ninguna |
> | Cholesky incompleta | $(R^TR)^{-1}$ | 2 triangulares | factorizar 1 vez |
> | Basado en estructura ([[Multigrid en 1D - transferencia e implementacion|multigrid]]) | - | $O(n)$ | depende del problema |

## Experimento numérico

Laplaciano 2D ($N=20$, $n=400$, [[5.7 - Comparacion metodos sistema 2D|Problema 5.7]]), GC con cada precondicionador hasta $\|r_k\|_2<10^{-10}$:

![[cap5-precondicionadores.png]]

**La curva "Ninguno" queda exactamente tapada por "Diagonal (Jacobi)"**: este `laplaciano2d` tiene diagonal constante (todas las entradas $=4$), así que $M_J=D^{-1}$ es un múltiplo escalar de $I$ y no cambia el condicionamiento relativo de nada - el mismo fenómeno, ya visto, de [[5.14 - Comparacion de los cuatro metodos|Problema 5.14]](b) con el problema 1D ("GCP-Jacobi = GC, exactamente"). Gauss-Seidel simétrico y Cholesky incompleta sí ayudan, y Cholesky incompleta es el mejor de los genéricos (consistente con que explota más información de $A$ que un simple barrido).

```julia
# requiere icholesky, laplaciano2d, diagonal_precon, sgs_precon -
# ver [[5.11 - GCP y Cholesky incompleta|Problema 5.11]], [[5.7 - Comparacion metodos sistema 2D|Problema 5.7]]

# versión instrumentada de gradconjp: igual al algoritmo original, guardando ‖r_k‖₂ en cada paso
function gradconjp_hist(A, b, precon, x0; tol=1e-10, max_iter=60)
    Mv(v) = precon isa Function ? precon(v) : precon*v
    x = copy(x0)
    r = b - A*x
    z = Mv(r)
    ρ = r'z
    resid = Float64[norm(r)]
    d = copy(z)
    ρold = 0.0
    iter = 0
    while norm(r) > tol && iter < max_iter
        if iter > 0
            β = ρ/ρold
            d = z + β*d
        end
        w = A*d
        α = ρ/(d'w)
        x += α*d
        r -= α*w
        z = Mv(r)
        ρold = ρ
        ρ = r'z
        push!(resid, norm(r))
        iter += 1
    end
    return x, resid
end

N = 20
A = laplaciano2d(N)
b = A * ones(N^2)
L = icholesky(A)

_, resid_none = gradconjp_hist(A, b, v -> v, zeros(N^2))
_, resid_jac  = gradconjp_hist(A, b, diagonal_precon(A), zeros(N^2))
_, resid_sgs  = gradconjp_hist(A, b, sgs_precon(A), zeros(N^2))
_, resid_ic   = gradconjp_hist(A, b, v -> L' \ (L \ v), zeros(N^2))

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10, xlabel="iteración", ylabel="‖r_k‖₂",
          title="Precondicionadores genéricos, laplaciano 2D (N=$N)")
lines!(ax, 0:length(resid_none)-1, resid_none, label="Ninguno")
lines!(ax, 0:length(resid_jac)-1, resid_jac, label="Diagonal (Jacobi)")
lines!(ax, 0:length(resid_sgs)-1, resid_sgs, label="Gauss-Seidel simétrico")
lines!(ax, 0:length(resid_ic)-1, resid_ic, label="Cholesky incompleta")
axislegend(ax)
fig
```

## Ver también
- [[Precondicionamiento y GCP]]
- [[Multigrid en 1D - transferencia e implementacion]]
- [[Factorizacion de Cholesky]] (Capítulo 2)
- [[Tablas Comparativas - Precondicionadores]] - estos tres genéricos comparados contra multigrid, con ventajas y desventajas
