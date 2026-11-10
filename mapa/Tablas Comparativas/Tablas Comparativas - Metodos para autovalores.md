---
tags: [comparacion, capitulo7]
---

# Métodos para autovalores: potencias vs. potencias inversas vs. Rayleigh vs. QR

> [!comparacion]
> | | Potencias | Potencias inversas (shift fijo $\mu$) | Cociente de Rayleigh (shift actualizado) | QR (Hessenberg + shift + deflación) |
> |---|---|---|---|---|
> | Qué calcula | autovalor **dominante** + autovector | el autovalor **más cercano a $\mu$** + autovector | **un** autovalor + autovector (no se elige cuál) | **todos** los autovalores (forma de Schur; autovectores en $P_k$) |
> | Hipótesis | $\vert \lambda_1\vert >\vert \lambda_2\vert $, $y_1^Hu_0\neq0$ | $\mu$ más cerca de $\lambda_*$ que del resto | $u_0$ razonable | Teo. 7.11 (sin shift); en la práctica, genérico |
> | Convergencia | lineal, $\vert \lambda_2/\lambda_1\vert $ | lineal, $\rho_\mu\propto\vert \lambda_*-\mu\vert $ | **cúbica** (hermítica), cuadrática (general) | lineal sin shift; ~cúbica con shift |
> | Costo por paso | 1 mat-vec | 2 triangulares ($O(n^2)$), tras **una** LU | 1 factorización nueva ($O(n^3)$ denso) | $O(n^2)$ (Hessenberg) |
> | Costo total | $k\cdot\text{nnz}(A)$, $k$ puede ser enorme | $O(n^3)+k\,O(n^2)$ | pocos pasos $\times O(n^3)$ | $O(n^3)$ todos |
> | Sirve para matrices enormes y ralas | **sí** | sí, con un solver iterativo | sí, con un solver iterativo | no (llena la matriz) |
> | Falla si | no hay dominante (par complejo conjugado, $\pm\lambda$) o $\vert \lambda_2\vert \approx\vert \lambda_1\vert $ | $\mu$ equidistante de dos autovalores | se va a un autovalor no deseado | shift de Rayleigh con espectro simétrico respecto de $a_{nn}$ ⟹ usar Wilkinson |

## Números de los problemas ($A_N=\text{tridiag}(-1,2,-1)$, tolerancia $10^{-10}$–$10^{-12}$)

| $N$ | Potencias | Pot. inversas $\mu=0$ | Pot. inversas, shift actualizado | QR sin shift | QR + Rayleigh | QR + Wilkinson |
|---|---|---|---|---|---|---|
| 5 | 96 | 17 | 4 | 64 | 45 | 8 |
| 15 | 641 | 17 | 7 (a $\lambda_3$) | 412 | se estanca | 33 |
| 50 | 5756 | 16 | 6 (a $\lambda_4$) | 3401 | se estanca | 102 |

(Potencias e inversas dan **un** autovalor; QR da los $N$.) Fuentes: [[7.2 - Potencias en la matriz A_N|7.2]], [[7.4 - Potencias inversas en la matriz A_N|7.4]], [[7.6 - Potencias inversas con y sin shift en A_N|7.6]], [[7.9 - QR con shift y deflacion|7.9]].

## Hermítica vs. general

| | $A$ hermítica | $A$ general |
|---|---|---|
| Autovalores | reales | complejos (pares conjugados si $A$ es real) |
| Autovectores | base ortonormal | pueden no formar base (defectuosa) |
| Cociente de Rayleigh | error **cuadrático** en el autovector | error lineal |
| Potencias: tasa del autovalor | $\vert \lambda_2/\lambda_1\vert ^2$ | $\vert \lambda_2/\lambda_1\vert $ |
| Rayleigh / QR con shift | cúbica | cuadrática |
| Límite de QR | **diagonal**; $P_k\to$ autovectores | triangular (Schur) |
| Hessenberg + QR | tridiagonal, se mantiene | Hessenberg llena |

## Ver también
- [[Capitulo 7 - MOC]]
- [[Tablas Comparativas - Indice]]
