---
tags: [capitulo5, comparacion]
---

# GC vs. GCP ([[Gradiente conjugado - definicion y propiedades|gradiente conjugado]] precondicionado)

> [!comparacion]
> | | GC | GCP |
> |---|---|---|
> | Inicialización | $r=b-Ax$, $ρ=r^Tr$ | $r=b-Ax$, $z=Mr$, $ρ=r^Tz$ |
> | Dirección | $d=r+\beta d$ | $d=z+\beta d$ |
> | Actualizar $r$ | $r=r-\alpha w$ | $r=r-\alpha w$, luego $z=Mr$ |
> | Escalar | $\rho=r^Tr=\|r\|_2^2$ | $\rho=r^Tz=\|r\|_M^2$ |
> | Parada | $\|r\|_2<\text{tol}\,\|b\|_2$ | $\|r\|_M<\text{tol}\,\|b\|_M$ |
> | Memoria | 4 vectores | 5 vectores |
> | Costo/iteración | 1 mat-vec | 1 mat-vec + 1 acción $Mξ$ |
> | Tasa | $\left(\dfrac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\right)^k$ | $\left(\dfrac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\right)^k$ |

GC es el caso particular $M=I$ de GCP (ahí $z=r$ y todo colapsa).

## Experimento numérico

Laplaciano 2D ($N=20$, $n=400$ - ver [[5.7 - Comparacion metodos sistema 2D|Problema 5.7]]), precondicionador: [[Precondicionadores genericos|Cholesky incompleta]] (`icholesky`, [[5.11 - GCP y Cholesky incompleta|Problema 5.11]]):

![[cap5-gc-vs-gcp.png]]

GC necesita 42 iteraciones para $\|r_k\|_2<10^{-10}$; GCP, 25 - consistente con la tabla real de [[5.11 - GCP y Cholesky incompleta|Problema 5.11]] (que usa $N$ más chicos, pero muestra la misma tendencia).

```julia
# requiere gradconjp (gradiente conjugado precondicionado), icholesky, laplaciano2d -
# ver [[Precondicionamiento y GCP]], [[5.11 - GCP y Cholesky incompleta|Problema 5.11]], [[5.7 - Comparacion metodos sistema 2D|Problema 5.7]]

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

_, resid_gc  = gradconjp_hist(A, b, v -> v, zeros(N^2))           # M=I: GC liso
L = icholesky(A)
_, resid_gcp = gradconjp_hist(A, b, v -> L' \ (L \ v), zeros(N^2))  # M=(R'R)⁻¹

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10, xlabel="iteración", ylabel="‖r_k‖₂",
          title="GC vs. GCP con Cholesky incompleta, laplaciano 2D (N=$N)")
lines!(ax, 0:length(resid_gc)-1, resid_gc, label="GC (sin precondicionar)")
lines!(ax, 0:length(resid_gcp)-1, resid_gcp, label="GCP (Cholesky incompleta)")
axislegend(ax)
fig
```

## Ver también
- [[Implementacion del gradiente conjugado]] · [[Precondicionamiento y GCP]]
