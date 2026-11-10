---
tags: [glosario, capitulo5]
aliases: [matrix-free, function handle]
---

# Function handle / matrix-free

## Definición

Pasarle a un algoritmo (GC, GMRes) la *acción* $\xi\mapsto A\xi$ como función de Julia en vez de la matriz $A$ en memoria. Sirve cuando $A$ ni siquiera conviene construirla explícitamente (por ejemplo, el laplaciano discreto: conocer la fórmula de la acción alcanza, no hace falta la matriz completa).

## Ver también
- [[5.6 - Funcion matvec tridiagonal|Problema 5.6]]
- [[Implementacion del gradiente conjugado]]
- [[Glosario - Indice|Glosario]]
