---
tags: [capitulo5, moc]
---

# Capítulo 5 - Índice de problemas

| # | Tema |
|---|---|
| [[5.1 - Raiz cuadrada de matriz sdp\|5.1]] | Existencia y unicidad de $A^{1/2}$ |
| [[5.2 - Identidad del funcional de energia\|5.2]] | $\phi(y)=\phi(x_*)+\frac12\|y-x_*\|_A^2$ |
| [[5.3 - Programar descenso y GC\|5.3]] | `descmaspron`, `gradconj` |
| [[5.4 - Curvas de nivel ejemplo 2x2\|5.4]] | Ejemplo $2×2$: SD vs. GC |
| [[5.5 - Comparacion metodos sistema 1D\|5.5]] | [[Metodo de Jacobi\|Jacobi]]/GS/[[Metodo SOR\|SOR]]/[[Metodo de Richardson\|Richardson]]/SD/GC en el problema modelo 1D |
| [[5.6 - Funcion matvec tridiagonal\|5.6]] | `matvec` sin construir $A$ |
| [[5.7 - Comparacion metodos sistema 2D\|5.7]] | Lo mismo que 5.5, en 2D (laplaciano de 5 puntos) |
| [[5.8 - Cota pesimista autovalor aislado\|5.8]] | Cota general vs. cota específica (autovalor aislado) |
| [[5.9 - Espectro de S2A y SAS\|5.9]] | $\sigma(S^2A)=\sigma(SAS)$ |
| [[5.10 - Cota de kappa2 SAS\|5.10]] | $\kappa_2(SAS)\leq c_2/c_1$ |
| [[5.11 - GCP y Cholesky incompleta\|5.11]] | `gradconjp`, comparación con [[Precondicionadores genericos\|Cholesky incompleta]] |
| [[5.12 - Multigrid como metodo iterativo\|5.12]] | Tests de consistencia, relación variacional, $O(n)$ |
| [[5.13 - Multigrid como precondicionador\|5.13]] | $M=$ ciclo en V; $κ₂(MA_h)$ acotado |
| [[5.14 - Comparacion de los cuatro metodos\|5.14]] | GC, GCP-Jacobi, [[Multigrid en 1D - transferencia e implementacion\|multigrid]], GCP-MG - tabla real del profesor |
| [[5.15 - Programar GMRes con Givens\|5.15]] | Implementación eficiente |
| [[5.16 - Programar GMRes con restart\|5.16]] | `GMRes(m)` - tabla real del profesor |
| [[5.17 - Comparar GC y GMRes\|5.17]] | En el sistema del Problema 5.5 |
| [[5.18 - Problema con conveccion\|5.18]] | [[GMRes - definicion y propiedad de minimizacion\|GMRes]], CGNE, CGNR - tabla real del profesor |
| [[5.19 - Demostrar Corolario 5.40\|5.19]] | $k$ autovalores distintos ⟹ GMRes termina en $≤k$ |
| [[5.20 - Demostrar Corolario 5.41\|5.20]] | $b$ combinación de $k$ autovectores ⟹ termina en $≤k$ |

Volver a [[Capitulo 5 - MOC]].
