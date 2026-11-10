---
tags: [capitulo3, comparacion, definicion]
---

# Estabilidad y número de condición ([[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]])

> [!idea]
> Número de condición $κ₂(A)=σ₁/σₘ$ (cociente de valores singulares extremos). Las ecuaciones normales elevan al cuadrado el condicionamiento ($κ₂(AᵀA)=κ₂(A)²$); QR por [[Matrices de Householder|Householder]] no amplifica nada, porque usa solo isometrías. Por eso QR es más estable, aunque más caro.

## Número de condición de $A$ rectangular

$$
\kappa_2(A) = \sqrt{\kappa_2(A^TA)} = \sqrt{\frac{\lambda_{\max}(A^TA)}{\lambda_{\min}(A^TA)}} = \frac{\sigma_1}{\sigma_m},
$$
con $\sigma_i$ los valores singulares de $A$ (raíces de los autovalores de $AᵀA$). Así, $\kappa_2(A^TA)=\kappa_2(A)^2$: **formar $AᵀA$ siempre empeora el condicionamiento**, elevándolo al cuadrado. Demostración completa: [[Lema 4.17]] (Capítulo 4).

## Estabilidad hacia atrás (backward stability)

Un algoritmo es *estable hacia atrás* si la solución calculada $\hat x$ (en punto flotante) es la solución **exacta** de un problema levemente perturbado ($\|\Delta A\|/\|A\|,\|\Delta b\|/\|b\|\cong\varepsilon_{maq}$). La relación general entre error hacia adelante ($\|\hat x-x\|$) y hacia atrás es
$$
\text{error hacia adelante} \cong \kappa_2(A)\times(\text{error hacia atrás}).
$$

## La comparación

> [!comparacion]
> | Método | Costo ($n≫m$) | Error hacia adelante típico |
> |---|---|---|
> | Ecuaciones normales ([[Factorizacion de Cholesky|Cholesky]] de $AᵀA$) | $\cong nm^2$ | $\cong\kappa_2(A)^2\,\varepsilon_{maq}$ |
> | QR por [[Ortogonalizacion de Gram-Schmidt|Gram-Schmidt]] | $\cong2nm^2$ | intermedio (peor que Householder si se usa la versión clásica) |
> | QR por Householder | $\cong2nm^2-\frac23m^3\cong2nm^2$ | $\cong\kappa_2(A)\,\varepsilon_{maq}$ (óptimo) |
>
> **Por qué Householder es estable:** $Q$ es producto de reflexiones, que son isometrías ($\|Qx\|_2=\|x\|_2$): no amplifican el error de redondeo. Las ecuaciones normales, en cambio, forman $AᵀA$ explícitamente, cuyo condicionamiento ya es $\kappa_2(A)^2$ **antes** de factorizar nada: aunque Cholesky se resuelva perfectamente, el error hacia adelante hereda ese $\kappa_2(A)^2$.

## Experimento numérico

Vandermonde con $a=\text{range}(1,w,\text{length}=200)$, $w$ creciente (mismo setup que [[3.6 - Estabilidad experimental minimos cuadrados|Problema 3.6]], variando $w$ para barrer $\kappa_2(A)$):

![[cap3-estabilidad-minimos-cuadrados.png]]

`mincuad` (ecuaciones normales) es consistentemente 2-3 órdenes de magnitud menos preciso que `mincuadqr` (Householder) en **todo** el rango de $\kappa_2(A)$ probado - la brecha absoluta confirma [[3.6 - Estabilidad experimental minimos cuadrados|Problema 3.6]], ahora con $\kappa_2(A)$ variando explícitamente en el eje $x$ en vez de fijo.

```julia
# requiere mincuadqr (de qrhouse/qrsolve) -
# ver [[Matrices de Householder]], [[3.4 - mincuadqr y ajuste polinomial|Problema 3.4]]

mincuad_normales(A, b) = cholesky(Symmetric(A'A)) \ (A'b)

n = 200
anchos = [2.0, 3.0, 4.5, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0]
kappas, err_normales, err_qr = Float64[], Float64[], Float64[]
for w in anchos
    a = range(1.0, w, length=n)
    A = Matrix([a.^3 a.^2 a ones(n)])
    xtrue = ones(4)
    b = A * xtrue
    push!(kappas, cond(A))
    push!(err_normales, norm(mincuad_normales(A, b) - xtrue))
    push!(err_qr, norm(mincuadqr(A, b)[1] - xtrue))
end

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="κ₂(A)", ylabel="‖x_computado - x_exacto‖₂",
          title="Estabilidad: ecuaciones normales vs. QR (Householder)")
scatterlines!(ax, kappas, err_normales, label="mincuad (ecuaciones normales)")
scatterlines!(ax, kappas, err_qr, label="mincuadqr (Householder)")
axislegend(ax, position=:lt)
fig
```

Un ejemplo todavía más extremo (hasta el punto de que las ecuaciones normales **no logran dar una respuesta**, no solo pierden precisión): la matriz de Läuchli en [[Parcial Modelo - Problema 2]], donde para $\varepsilon\leq10^{-8}$ la factorización de Cholesky de $A^TA$ falla directamente porque $A^TA$ calculada en punto flotante resulta singular.

## Regla práctica

- $A$ bien condicionada y se necesita velocidad: ecuaciones normales (la mitad de operaciones).
- $A$ mal condicionada, o se necesita precisión: QR por Householder.

## Ver también
- [[Minimos cuadrados y ecuaciones normales]]
- [[Matrices de Householder]]
- En el Capítulo 4: [[Numero de condicion y error relativo]] generaliza esta idea a sistemas cuadrados.
