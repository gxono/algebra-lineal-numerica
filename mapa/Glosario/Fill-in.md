---
tags: [glosario, capitulo2, capitulo5]
---

# Fill-in

## Definición

Entradas que eran cero en $A$ y dejan de serlo durante una factorización (Gauss/LU/Cholesky). En una matriz [[Sparse|sparse]] sin estructura especial, el fill-in puede ser tan severo que la matriz factorizada termine casi densa - por eso conviene explotar la estructura (tridiagonal, banda) o usar métodos iterativos en vez de factorizar. [[Precondicionadores genericos|Cholesky incompleta]] es, justamente, Cholesky que descarta el fill-in a propósito.

## Ver también
- [[2.5 - Factorizacion LU tridiagonal|Problema 2.5]]
- [[Precondicionadores genericos]]
- [[Glosario - Indice|Glosario]]
