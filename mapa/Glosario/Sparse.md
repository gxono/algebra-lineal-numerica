---
tags: [glosario, capitulo2, capitulo5]
---

# Sparse (matriz)

## Definición

Matriz donde la gran mayoría de las entradas son cero (típicamente $O(n)$ no-nulos en vez de los $O(n^2)$ de una matriz densa) - por ejemplo la tridiagonal `tridiag(-1,2,-1)` que aparece en varios capítulos. Importa porque factorizar una matriz sparse puede generar [[Fill-in|fill-in]] y destruir esa estructura; por eso para sistemas grandes y sparse suele convenir un método iterativo (GC, GMRes) o uno que respete la estructura (tridiagonal, multigrid) antes que factorizar directo.

## Ver también
- [[Fill-in]]
- [[Resolver sistemas lineales en Julia]]
- [[Glosario - Indice|Glosario]]
