---
tags: [capitulo1]
---

# Resolver sistemas lineales en Julia: el operador `\`

> [!idea]
> `x = A \ b` resuelve $Ax=b$ eligiendo automáticamente el mejor algoritmo según el *tipo* de $A$ (triangular, tridiagonal, simétrica, sdp, densa genérica, rectangular → [[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]]). Es el punto de comparación de todo el curso: los capítulos 2 a 5 son, en el fondo, "qué hace `\` por dentro, y cuándo conviene escribir algo distinto a mano".

## Ejemplo básico

```julia
using LinearAlgebra

A = [3 5 1; 1 -6 1; 1 0 1]
b = [8, 4, 2]
x = A \ b
```

## `A \ b` vs. `inv(A) * b`

```julia
using BenchmarkTools

A = rand(10,10)
b = rand(10)
@benchmark A \ b          # resuelve el sistema directamente (rápido, estable)
@benchmark inv(A) * b     # calcula la inversa completa y multiplica (más lento, menos estable)
```

> [!pendiente]
> **Nunca calcular `inv(A)` para resolver un sistema.** `A \ b` es más rápido (no hace falta la inversa completa, basta con una factorización) y más estable numéricamente (ver [[Estabilidad y numero de condicion|estabilidad hacia atrás]] en el Capítulo 3). Esto se repite como principio en todo el curso: la inversa explícita casi nunca es necesaria ni conveniente.

## Ejemplo con matriz tridiagonal (preludio del Capítulo 2 y 4)

```julia
using CairoMakie

N = 100
A = zeros(N,N)
A[1,1] = 1
A[N,N] = 1
for i in 2:N-1
    A[i, i-1:i+1] = [-1, 2, -1]
end
b = vcat(0, fill(1/N^2, N-2), 0)

x = A \ b
lines(1:N, x)
```

Esta matriz tridiagonal `tridiag(-1,2,-1)` reaparece constantemente: es la discretización de $-u''=f$ (ver [[Multigrid en 1D - problema modelo y suavizado]] en el Capítulo 5) y el ejemplo canónico de matriz sparse, sdp y mal condicionada.

## Ver también
- [[Factorizacion LU]] (qué hace `\` internamente para matrices generales)
- [[Factorizacion de Cholesky]] (qué hace `\` para matrices sdp)
- [[Referencia de Julia]]
