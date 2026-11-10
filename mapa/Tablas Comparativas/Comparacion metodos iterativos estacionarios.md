---
tags: [capitulo4, comparacion]
---

# Comparación: métodos iterativos estacionarios (Capítulo 4)

> [!comparacion]
> | Método | $M$ (matriz de iteración) | Usa valores nuevos en la misma iteración | Condición suficiente de convergencia |
> |---|---|---|---|
> | [[Metodo de Jacobi\|Jacobi]] | $-D^{-1}(L+U)$ | No | $A$ e.d.d. ($\|M_J\|_\infty<1$) |
> | [[Metodo de Gauss-Seidel\|Gauss-Seidel]] | $-(L+D)^{-1}U$ | Sí | $A$ e.d.d., **o** $A$ sdp (converge siempre) |
> | [[Metodo SOR\|SOR]] | $(D/\omega+L)^{-1}\big(\frac{1-\omega}\omega D-U\big)$ | Sí, + relajación $ω$ | $A$ sdp y $\omega\in(0,2)$ (Ostrowski-Reich); **necesario** $\omega\in(0,2)$ siempre |
> | [[Metodo de Richardson\|Richardson]] | $I-\alpha A$ | - | $A$ sdp y $0<\alpha<2/\lambda_{\max}$; óptimo $\alpha^*=2/(\lambda_{\min}+\lambda_{\max})$ |

## Velocidad relativa (regla empírica, no una garantía general)

En los ejemplos trabajados ([[4.15 - Jacobi Gauss-Seidel SOR en sistemas concretos|Problema 4.15]]), con $ω$ bien elegido:
$$
\text{SOR (}\omega\text{ óptimo)} \ \text{más rápido que}\ \text{Gauss-Seidel} \ \text{más rápido que}\ \text{Jacobi}.
$$
Esto **no** es un teorema general (depende de la matriz), pero es la tendencia típica: cada método "ve" información más actualizada que el anterior.

## Lo que todos tienen en común: $ρ(M)$ manda

La velocidad de convergencia asintótica de **cualquiera** de estos métodos es $\rho(M)^k$: cuantas más cifras decimales se quieran ganar, más iteraciones hacen falta, en una relación logarítmica (análoga a la Observación 5.20 del Capítulo 5 para el [[Gradiente conjugado - definicion y propiedades|gradiente conjugado]]). Lo que diferencia a los métodos es **qué tan chico** logran hacer $ρ(M)$ para una $A$ dada.

## Experimento numérico

Problema modelo 1D ($A_h$, $N=100$), $\alpha$ de [[Metodo de Richardson|Richardson]] óptimo ($\lambda_{\min}$, $\lambda_{\max}$ exactos vía `eigmin`/`eigmax`):

![[cap4-comparacion-iterativos.png]]

**[[Metodo de Jacobi|Jacobi]] y Richardson óptimo son prácticamente indistinguibles** en este problema: para $A_h$ (diagonal constante $2/h^2$), el paso implícito de Jacobi ($1/a_{ii}=h^2/2$) y el $\alpha$ óptimo de Richardson ($2/(\lambda_{\min}+\lambda_{\max})\approx h^2/2$ cuando $\lambda_{\max}\gg\lambda_{\min}$) coinciden casi exactamente - Jacobi *es*, en este caso particular, prácticamente Richardson con el $\alpha$ óptimo. [[Metodo de Gauss-Seidel|Gauss-Seidel]] mejora un poco sobre ambos. [[Metodo SOR|SOR]] con $\omega=1.9$ es dramáticamente mejor que los tres - el comportamiento típico de SOR bien afinado sobre un laplaciano discreto.

```julia
# requiere jacobi, gauss_seidel, SOR, richardson, Ah -
# ver [[Metodo de Jacobi]], [[Metodo de Gauss-Seidel]], [[Metodo SOR]], [[Metodo de Richardson]], [[Multigrid en 1D - problema modelo y suavizado]]

N = 100
A = Ah(N)
b = A * ones(N)
λmin, λmax = eigmin(A), eigmax(A)
αopt = 2/(λmin + λmax)

_, errJ  = jacobi(A, b; tol=1e-8, maxiter=3000)
_, errGS = gauss_seidel(A, b; tol=1e-8, maxiter=3000)
_, errSOR = SOR(A, b, 1.9; tol=1e-8, maxiter=3000)
_, errR  = richardson(A, b, αopt; tol=1e-8, maxiter=3000)

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10, xlabel="iteración", ylabel="‖x_{k+1}-x_k‖₁",
          title="Métodos iterativos estacionarios, problema modelo 1D (N=$N)")
lines!(ax, errJ, label="Jacobi")
lines!(ax, errGS, label="Gauss-Seidel")
lines!(ax, errSOR, label="SOR (ω=1.9)")
lines!(ax, errR, label="Richardson (α óptimo)")
axislegend(ax)
fig
```

## El límite de todos los métodos estacionarios

Ninguno de estos métodos cambia $M$ entre iteraciones. El [[Capitulo 5 - MOC|Capítulo 5]] completo nace de preguntarse qué se gana si se permite que la "dirección de búsqueda" cambie en cada paso (steepest descent, gradiente conjugado): la respuesta es una mejora drástica en la tasa de convergencia, especialmente para matrices mal condicionadas.

## Ver también
- [[Capitulo 4 - MOC]]
- [[Tablas Comparativas - Indice]]
