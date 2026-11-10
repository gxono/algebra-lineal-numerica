---
tags: [capitulo1, definicion]
---

# Sintaxis básica de Julia

> [!idea]
> Vectores son columna por convención en todo el curso. $A[i,:]$ es fila $i$, $A[:,j]$ es columna $j$. `tril`/`triu`/`diag`/`Diagonal` se usan todo el tiempo desde el Capítulo 2.

## Matrices y vectores

```julia
A = [1 2 3; 4 5 6; 7 8 9; 10 11 12]   # matriz 4x3, separador de filas ";"
v = -8:5                       # range (no vector denso)
w = range(-8, 5, 10)            # 10 puntos equiespaciados entre -8 y 5
collect(v)                      # materializa el range en un Vector

size(A)                         # (4,3)
A[2,3]                           # elemento fila 2, columna 3
A[3,:]                           # fila 3 (como vector)
A[:,3]                           # columna 3
A[[2,4], :]                      # filas 2 y 4
A[2:4, :]                        # filas 2 a 4
```

## Concatenar y operar

```julia
[A b]           # concatenación horizontal (hcat) - requiere b columna del mismo n° de filas
[A; b']         # concatenación vertical (vcat)
A .+ 1          # suma elemento a elemento (broadcasting)
A .^ 2          # potencia elemento a elemento (≠ A^2, que es producto matricial)
```

> [!pendiente]
> Ojo con la diferencia $A^2$ (producto matricial $A*A$) vs $A.^2$ (cada entrada al cuadrado). Es un error típico de parcial.

## Matrices especiales

```julia
I(4)               # identidad 4x4 (tipo UniformScaling/Diagonal)
Diagonal(v)        # matriz diagonal eficiente a partir de un vector
diagm(v)           # matriz diagonal densa
diagm(1 => v)      # v en la sobre-diagonal
triu(A), tril(A)   # parte triangular superior/inferior (incluye diagonal)
triu(A,1)          # estrictamente superior (sin diagonal)
tril(A,-1)         # estrictamente inferior (sin diagonal)
```

> [!pendiente]
> **Diferencia clave con Octave:** en Octave, `diag` tiene doble propósito (extrae la diagonal de una matriz, o construye una matriz diagonal a partir de un vector). En Julia, `diag` **solo extrae** la diagonal; para construir hay que usar `Diagonal(v)` o `diagm(v)`.

## Funciones y despacho

```julia
factorial(n)                 # ya viene en Base
comb(n, k) = factorial(n) / (factorial(k) * factorial(n - k))   # función de una línea
```

Julia usa *despacho múltiple*: se pueden definir varios métodos para la misma función según el tipo de los argumentos (por ejemplo, `eleu(A::Tridiagonal)` vs `eleu(A::Bidiagonal)` en el [[Capitulo 2 - MOC|Capítulo 2]]), y Julia elige automáticamente la versión más específica y más eficiente según el tipo de la matriz.

## Ver también
- [[Resolver sistemas lineales en Julia]]
- [[Graficos con CairoMakie]]
- [[Referencia de Julia]] (qué función usar para cada cuenta, built-in vs. implementación propia, Symbolics.jl)
