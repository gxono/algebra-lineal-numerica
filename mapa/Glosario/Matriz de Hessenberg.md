---
tags: [glosario, capitulo3, capitulo5, capitulo7]
aliases: [Hessenberg]
---

# Matriz de Hessenberg

## Definición

Matriz (casi) triangular superior: todas las entradas por debajo de la subdiagonal son cero ($h_{ij}=0$ si $i>j+1$). Aparece en GMRes porque el proceso de Arnoldi produce exactamente una matriz así, lo que permite triangularla con una sola rotación de [[Rotaciones de Givens|Givens]] por iteración en vez de recalcular todo. En el Capítulo 7 es la forma a la que se reduce $A$ (por semejanza, con Householder, $O(n^3)$) antes del algoritmo QR: el paso QR la **conserva** y cuesta $O(n^2)$. Hessenberg + hermítica = tridiagonal.

**Ojo:** Hessenberg × triangular es Hessenberg, pero Hessenberg × Hessenberg **no** ([[7.13 - Producto de matrices Hessenberg|Problema 7.13]]).

## Ver también
- [[GMRes - implementacion y rotaciones de Givens]]
- [[3.7 - Conteo Givens en Hessenberg|Problema 3.7]]
- [[Algoritmo QR - Hessenberg, shifts y deflacion]]
- [[7.14 - QR preserva Hessenberg|Problema 7.14]]
- [[Glosario - Indice|Glosario]]
