---
tags: [capitulo5, comparacion]
---

# GC vs. [[Descenso mas pronunciado|descenso más pronunciado]] (SD)

> [!comparacion]
> | Método | Cota del factor de reducción por iteración | $κ=100$ |
> |---|---|---|
> | [[Metodo de Richardson|Richardson]] óptimo / SD | $\left(\dfrac{\kappa-1}{\kappa+1}\right)^k$ | $(0.980)^k$ |
> | GC (Chebyshev) | $2\left(\dfrac{\sqrt\kappa-1}{\sqrt\kappa+1}\right)^k$ | $2\,(0.818)^k$ |

Para reducir el error en un factor $10^{-2}$ con $κ=100$: SD necesita $\approx230$ iteraciones, GC necesita $\approx30$. La diferencia entre $κ$ y $\sqrt κ$ en el exponente crece con el mal condicionamiento.

**GC es, en cada paso, al menos tan bueno como SD:** el afín sobre el que minimiza SD ($x_k+\alpha r_k$) está contenido en el afín sobre el que minimiza GC ($x_0+\mathcal K_{k+1}$) - ver [[Gradiente conjugado - definicion y propiedades]].

## Experimento numérico

Problema modelo 1D ($A_h$, $N=100$, $\kappa_2(A_h)\approx4134$), $x_0=0$, $x_*=\mathbb{1}$:

![[cap5-gc-vs-sd.png]]

SD apenas se mueve en 300 iteraciones (factor $\approx0.9995$ por paso, consistente con $\frac{\kappa-1}{\kappa+1}$); GC decae lento al principio y después cae en picada a precisión de máquina alrededor de la iteración 50 - la combinación de la cota de Chebyshev (decrece $\approx0.97$ por paso) con la **terminación finita** ([[GC como metodo directo]]): con $N=100$ y buena parte del espectro "cubierto" por el espacio de Krylov, el método colapsa de golpe en vez de decaer suavemente hasta el final.

```julia
# requiere descmaspron, gradconj (gradiente conjugado), Ah -
# ver [[Descenso mas pronunciado]], [[Implementacion del gradiente conjugado]], [[Multigrid en 1D - problema modelo y suavizado]]

# versiones instrumentadas: igual al algoritmo original, guardando ‖x_k - x_*‖_A en cada paso
function descmaspron_hist(A, b, x0, xstar; tol=1e-12, max_iter=300)
    x = copy(x0)
    r = b - A*x
    errA = Float64[]
    iter = 0
    while norm(r) > tol && iter < max_iter
        push!(errA, sqrt((x - xstar)' * A * (x - xstar)))
        Ar = A*r
        α = (r'r) / (r'Ar)
        x += α*r
        r -= α*Ar
        iter += 1
    end
    return x, errA
end

function gradconj_hist(A, b, x0, xstar; tol=1e-12, max_iter=300)
    x = copy(x0)
    r = b - A*x
    ρ = r'r
    errA = Float64[]
    d = copy(r)
    ρold = 0.0
    iter = 0
    while sqrt(ρ) > tol && iter < max_iter
        push!(errA, sqrt((x - xstar)' * A * (x - xstar)))
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
        iter += 1
    end
    return x, errA
end

N = 100
A = Ah(N)
xstar = ones(N)
b = A * xstar

_, errA_sd = descmaspron_hist(A, b, zeros(N), xstar)
_, errA_gc = gradconj_hist(A, b, zeros(N), xstar)

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10, xlabel="iteración", ylabel="‖x_k - x*‖_A",
          title="Descenso más pronunciado vs. gradiente conjugado (N=$N)")
lines!(ax, errA_sd, label="Descenso más pronunciado")
lines!(ax, errA_gc, label="Gradiente conjugado")
axislegend(ax)
fig
```

## Ver también
- [[Descenso mas pronunciado]] · [[Cotas de error del gradiente conjugado]]
