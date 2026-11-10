---
tags: [capitulo2, comparacion]
---

# Costo computacional: métodos directos (Capítulo 2)

> [!comparacion]
> | Método | Costo (orden dominante) | Notas |
> |---|---|---|
> | Sustitución (triangular ya factorizada) | $n²$ | [[Matrices triangulares - sustitucion]] |
> | Eliminación de Gauss (resolver 1 sistema) | $2n³/3$ | [[Metodo de Gauss - eliminacion]] |
> | Factorización [[Factorizacion LU|LU]] (una vez) | $2n³/3$ | luego cada sistema adicional cuesta solo $2n²$ |
> | Factorización de [[Factorizacion de Cholesky|Cholesky]] ($A$ sdp) | $n³/3$ | mitad que LU, explota la simetría |
> | $k$ sistemas con Gauss repetido | $k·(2n³/3)$ | sin reusar factorización |
> | $k$ sistemas, factorizando 1 vez (LU) | $2n³/3 + k·n²$ | gana para $k$ grande |
> | Sistema tridiagonal (LU especializado) | $O(n)$ | ver `tridilu`, [[2.5 - Factorizacion LU tridiagonal|Problema 2.5]] |

## La pregunta que conecta todo el capítulo

**¿Por qué no usar siempre Cholesky?** Porque exige que $A$ sea sdp. Si $A$ es general, hace falta LU con [[Pivoteo y matrices de permutacion|pivoteo]]. Si $A$ no es cuadrada (sistema rectangular, [[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]]), hace falta QR (Capítulo 3) en vez de Cholesky directamente sobre $A$, por motivos de **estabilidad** - ver [[Estabilidad y numero de condicion]].

## Experimento numérico: tiempo real, no solo flops

Benchmark real (`@belapsed`) de `metodo_de_gauss` (denso), `mi_cholesky!` (sdp) y `tridilu` (tridiagonal) para $n$ creciente:

![[cap2-costo-metodos-directos.png]]

Gauss y Cholesky siguen de cerca la recta de referencia $n^3$ (log-log); `tridilu` ni se nota en la misma escala - para $n=2000$ sigue corriendo en microsegundos, mientras que a Gauss/Cholesky ya les cuesta llegar limpiamente a $n=300$. Confirma empíricamente la tabla de arriba, no solo en el conteo teórico de operaciones sino en tiempo de reloj real.

```julia
# requiere metodo_de_gauss, mi_cholesky!, tridilu -
# ver [[Metodo de Gauss - eliminacion]], [[Factorizacion de Cholesky]], [[2.5 - Factorizacion LU tridiagonal|Problema 2.5]]
using BenchmarkTools

ns_densas = [30, 50, 80, 120, 170, 230, 300]
ns_tridi  = [30, 50, 80, 120, 170, 230, 300, 500, 800, 1300, 2000]

function tiempo_chol(n)
    A = rand(n,n)
    A = A*A' + n*I
    return @belapsed mi_cholesky!(copy($A))
end

function tiempo_tridi(n)
    A = Tridiagonal(fill(-1.0,n-1), fill(2.0,n), fill(-1.0,n-1))
    return @belapsed tridilu($A)
end

t_gauss = [(@belapsed metodo_de_gauss($(rand(n,n)+n*I), $(rand(n)))) for n in ns_densas]
t_chol  = [tiempo_chol(n) for n in ns_densas]
t_tridi = [tiempo_tridi(n) for n in ns_tridi]

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="n", ylabel="tiempo (s)",
          title="Costo real: Gauss vs. Cholesky vs. LU tridiagonal")
scatterlines!(ax, ns_densas, t_gauss, label="metodo_de_gauss (denso, O(2n³/3))")
scatterlines!(ax, ns_densas, t_chol, label="mi_cholesky! (sdp, O(n³/3))")
scatterlines!(ax, ns_tridi, t_tridi, label="tridilu (tridiagonal, O(n))")
axislegend(ax, position=:lt)
fig
```

## Ver también
- [[Capitulo 2 - MOC]]
- [[Tablas Comparativas - Indice]]
