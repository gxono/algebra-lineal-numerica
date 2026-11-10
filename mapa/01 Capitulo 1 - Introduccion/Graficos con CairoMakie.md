---
tags: [capitulo1]
---

# Gráficos con CairoMakie

> [!idea]
> `lines(x, y)` / `scatter(x, y)` crean una figura nueva y devuelven `(fig, ax, plt)`; `lines!(ax, x, y)` / `scatter!(ax, x, y)` agregan sobre un `Axis` existente (el `!` es la convención de Julia para "función que modifica/muta algo", acá modifica la figura). Se usa todo el curso para graficar error vs. iteraciones, curvas de nivel, etc.

## Curvas básicas

```julia
using CairoMakie

x = 0:0.01:2π
lines(x, sin.(x))                       # función 1D

t = 0:0.01:2π
lines(cos.(t), sin.(t))                 # curva paramétrica (círculo)

t = 0:0.01:6π
lines(cos.(t), sin.(t), t)              # curva 3D (hélice) - Makie detecta Axis3 automáticamente
```

## Varias curvas, títulos, leyenda

```julia
x = 0:0.1:3
fig = Figure()
ax = Axis(fig[1, 1], title="Recta y parábola", xlabel="x", ylabel="y")
lines!(ax, x, x .+ 1, label="y = x + 1")
lines!(ax, x, x .^ 2, label="y = x²")
axislegend(ax)
fig
```

## Escala logarítmica (muy usado en capítulos 2, 4 y 5 para graficar error vs. iteraciones)

```julia
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10, xlabel="iteración", ylabel="error")
scatterlines!(ax, 1:length(err), err)   # "err", no "error": ese nombre ya es una función de Base
fig
```
> [!pendiente]
> No llamar `error` a la variable que guarda el historial de error - `error` ya es la función de Base que lanza excepciones (`error("mensaje")`). Si se la sobreescribe sin querer, cualquier código que dependa de tirar una excepción deja de funcionar, y el propio `length(error)` tira `MethodError` porque `error` sigue siendo una función, no un vector. Convención usada en todo el resto de la vault: `err`.

## Curvas de nivel (`contour`)

```julia
fig = Figure()
ax = Axis(fig[1, 1])
contour!(ax, xs, ys, φ; levels=25)      # φ(x,y), xs/ys rangos
scatter!(ax, [x_star[1]], [x_star[2]])  # marcar un punto, ej. la solución
fig
```

## Rectas de referencia (`hlines!`/`vlines!`)

```julia
hlines!(ax, [valor], color=:red, linestyle=:dash, label="...")   # recta horizontal
vlines!(ax, [valor], color=:gray, linestyle=:dot)                 # recta vertical
```
Útil para marcar cotas teóricas sobre un gráfico (ej. $1/\lambda_{\min}$, $1/\lambda_{\max}$ en [[Parcial Modelo - Problema 4|Parcial Modelo - Problema 4(f)]]).

> [!pendiente]
> `lines!` **no** acepta `marker` como atributo (da error) - para puntos marcados sobre la curva, usar `scatterlines!` (línea + marcadores) o `scatter!` (solo marcadores) en su lugar.

## Ver también
- [[Comparacion metodos iterativos estacionarios]] (gráficos de error de [[Metodo de Jacobi|Jacobi]]/[[Metodo de Gauss-Seidel|Gauss-Seidel]]/[[Metodo SOR|SOR]], Capítulo 4)
- [[Tablas Comparativas - GC vs SD]] (Capítulo 5)
- [[Referencia de Julia]]
