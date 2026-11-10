---
tags: [capitulo1]
---

# Sentencias de control en Julia

> [!idea]
> `for`/`while` normales, `n < 0 && return 0` como atajo de `if` de una línea, y todas las funciones de los capítulos 2-5 están escritas con este estilo: nombre, docstring con la fórmula/algoritmo, cuerpo con bucles explícitos tipados.

## Funciones, recursión y atajos

```julia
function fact1(n)
    n < 0 && return 0        # atajo: "si se cumple la condición, evaluar lo de la derecha"
    iszero(n) && return 1

    f = 1
    for i in 2:n
        f *= i
    end
    return f
end

fact2(n) = n < 0 ? 0 : reduce(*, 1:n)   # equivalente, con reduce
fact3 = factorial                        # ya viene en Base
```

## `while` con criterio de parada (patrón que se repite en todo el curso: Newton, punto fijo, [[Metodo de Jacobi|Jacobi]], GC, [[GMRes - definicion y propiedad de minimizacion|GMRes]]...)

```julia
function puntofijo_sqrt2(tol=1e-8)
    iter(n) = 0.5 * (n + 2/n)          # iteración de punto fijo para sqrt(2)
    old, new = 1, iter(1)
    n_iter = 0
    while abs(old - new) > tol
        old, new = new, iter(new)
        n_iter += 1
    end
    return new, n_iter
end
```
*(Envuelto en una función a propósito: un `while` que reasigna variables en el **nivel superior** de un script - no dentro de una función - choca con una regla de scoping de Julia ["soft vs. hard global scope"](https://docs.julialang.org/en/v1/manual/variables-and-scoping/#On-soft-scope) y tira `UndefVarError` al correrlo como archivo, aunque funcione bien pegado directamente en el REPL. Dentro de una función no hay ambigüedad.)*

Este patrón - `while` que compara una aproximación nueva contra la vieja, con una tolerancia - es exactamente la estructura de **todos** los métodos iterativos del curso: [[Metodo de Jacobi]], [[Metodo de Gauss-Seidel]], [[Implementacion del gradiente conjugado|gradiente conjugado]], [[GMRes - implementacion y rotaciones de Givens|GMRes]], etc.

## Ver también
- [[Sintaxis basica de Julia]]
