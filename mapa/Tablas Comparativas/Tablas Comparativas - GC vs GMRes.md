---
tags: [capitulo5, comparacion]
---

# GC vs. [[GMRes - definicion y propiedad de minimizacion|GMRes]] (y costo por iteración)

> [!comparacion]
> | | GC | GMRes |
> |---|---|---|
> | Hipótesis | $A$ sdp | $A$ no singular |
> | Minimiza | $\|x-x_*\|_A$ (error, norma $A$) | $\|b-Ax\|_2$ (residuo, norma 2) |
> | Polinomios | $\|x_k-x_*\|_A=\min_p\|p(A)(x_0-x_*)\|_A$ | $\|r_k\|_2=\min_p\|p(A)r_0\|_2$ |
> | Terminación | $≤n$ pasos | $≤n$ pasos |
> | Recurrencia | dos términos | no - guarda toda la base |
> | Memoria/iter. | $O(n)$, fija | $O(nk)$, crece con $k$ |

## Costo por iteración, en detalle

> [!comparacion]
> | | GC ($A$ sdp) | GMRes (iteración $k$) |
> |---|---|---|
> | Productos mat-vec | 1 | 1 |
> | Productos escalares | 2 | $k+1$ |
> | Operaciones axpy | 3 | $k$ |
> | Rotaciones de [[Rotaciones de Givens|Givens]] | - | $k$, $O(k)$ |
> | Memoria | 4 vectores | $k+1$ vectores + $O(k^2)$ para $H$ |
> | Costo total hasta iter. $k$ | $O(nk)$ | $O(nk^2)$ |

La diferencia la da la recurrencia de dos términos, que GC obtiene gratis por ser $A$ sdp. Sin ella, cada iteración de GMRes es más cara que la anterior - por eso en GMRes es **todavía más importante** precondicionar bien: con pocas iteraciones es la única manera de mantener acotados memoria y costo.

## Experimento numérico

Problema modelo 1D ($A_h$ sdp, $N=100$) - el caso en el que **ambos** métodos son aplicables:

![[cap5-gc-vs-gmres.png]]

Las curvas son prácticamente indistinguibles hasta que GC "colapsa" por terminación finita un par de iteraciones antes que GMRes - confirma lo que dice la tabla de arriba: para $A$ sdp, **la cantidad de iteraciones es comparable**; la diferencia real está en el costo por iteración (ver [[5.17 - Comparar GC y GMRes|Problema 5.17]]).

```julia
# requiere Ah - ver [[Multigrid en 1D - problema modelo y suavizado]]
# versiones instrumentadas de gradconj y gmres: igual a los algoritmos originales
# ([[Implementacion del gradiente conjugado]], [[GMRes - implementacion y rotaciones de Givens]]),
# guardando ‖r_k‖₂ en cada paso

function gradconj_resid_hist(A, b, x0; tol=1e-10, max_iter=150)
    x = copy(x0)
    r = b - A*x
    ρ = r'r
    resid = Float64[norm(r)]
    d = copy(r)
    ρold = 0.0
    iter = 0
    while norm(r) > tol && iter < max_iter
        if iter > 0
            β = ρ/ρold
            d = r + β*d
        end
        w = A*d
        α = ρ/(d'w)
        x += α*d
        r -= α*w
        ρold = ρ
        ρ = r'r
        push!(resid, norm(r))
        iter += 1
    end
    return x, resid
end

function gmres_resid_hist(A, b, x0; tol=1e-10, kmax=length(b))
    r0 = b - A*x0
    β = ρ = norm(r0)
    V = zeros(length(b), kmax+1)
    V[:,1] = r0/ρ
    H = zeros(kmax+1, kmax)
    g = zeros(kmax+1)
    g[1] = β
    c = zeros(kmax)
    s = zeros(kmax)
    resid = Float64[ρ]
    k = 0
    while ρ > tol*norm(b) && k < kmax
        k += 1
        w = A*V[:,k]
        for j in 1:k
            H[j,k] = w'V[:,j]
            w -= H[j,k]*V[:,j]
        end
        H[k+1,k] = norm(w)
        V[:,k+1] = w/H[k+1,k]
        for j in 1:k-1
            H[j:j+1,k] = [c[j] -s[j]; s[j] c[j]] * H[j:j+1,k]
        end
        ν = hypot(H[k,k], H[k+1,k])
        c[k] = H[k,k]/ν
        s[k] = -H[k+1,k]/ν
        H[k,k] = ν
        H[k+1,k] = 0
        g[k:k+1] = [c[k] -s[k]; s[k] c[k]] * g[k:k+1]
        ρ = abs(g[k+1])
        push!(resid, ρ)
    end
    y = UpperTriangular(H[1:k,1:k]) \ g[1:k]
    return x0 + V[:,1:k]*y, resid
end

N = 100
A = Ah(N)
b = A * ones(N)

_, resid_gc    = gradconj_resid_hist(A, b, zeros(N))
_, resid_gmres = gmres_resid_hist(A, b, zeros(N); kmax=150)

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10, xlabel="iteración", ylabel="‖r_k‖₂ / ‖b‖₂",
          title="GC vs. GMRes, A sdp (problema modelo 1D, N=$N)")
lines!(ax, 0:length(resid_gc)-1, resid_gc ./ norm(b), label="GC")
lines!(ax, 0:length(resid_gmres)-1, resid_gmres ./ norm(b), label="GMRes")
axislegend(ax)
fig
```

## Ver también
- [[Implementacion del gradiente conjugado]] · [[GMRes - implementacion y rotaciones de Givens]]
