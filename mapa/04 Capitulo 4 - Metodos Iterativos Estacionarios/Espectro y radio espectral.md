---
tags: [capitulo4, definicion, teorema]
---

# Espectro y radio espectral

> [!idea]
> $ρ(A)=\max\{|\lambda|:\lambda\in\sigma(A)\}$. Es la cantidad que **realmente** determina si un método iterativo estacionario converge (no cualquier norma de $M$). Siempre $ρ(A)≤‖A‖$ para cualquier norma inducida, pero puede ser estrictamente menor - por eso $ρ$ es la caracterización exacta, y una norma chica es solo una condición *suficiente*.

## Definición

$\rho(A) := \max\{|\lambda| : \lambda \text{ autovalor de } A\}$.

## $ρ(A)≤‖A‖$ para toda norma inducida

Demostración completa: [[Radio espectral acotado por norma inducida]]. De hecho vale algo más fuerte - $\rho(A)$ es el **ínfimo** de $\|A\|$ sobre todas las normas inducidas (puede acercarse tanto como se quiera, aunque en general no se alcanza): [[Teorema 4.20]].

> [!pendiente]
> La desigualdad puede ser **estricta**: $ρ(A)$ puede valer $0$ aunque $A≠0$ (ejemplo: $A$ triangular con diagonal nula y algún elemento fuera de la diagonal no nulo - sus autovalores son los de la diagonal, todos cero, pero $‖A‖≠0$). Ver [[4.7 - Radio espectral y norma|Problema 4.7]].

## Por qué importa: caracterización exacta de convergencia

Para la iteración $x_{k+1}=Mx_k+c$, el error $e_k=x_k-x_*$ satisface $e_{k+1}=Me_k$, así que $e_k=M^ke_0$. La iteración converge **para todo $x_0$** si y solo si $M^k\to0$, lo cual ocurre **si y solo si $ρ(M)<1$** (no hace falta que exista una norma con $‖M‖<1$, aunque si existe, alcanza) - las cuatro formas de decir esto son equivalentes, demostración completa: [[Teorema 4.23]]. Más fuerte todavía, y con el mismo sabor, la **fórmula de Gelfand**:
$$
\rho(A) = \lim_{k\to\infty}\|A^k\|^{1/k}
$$
para cualquier norma, lo cual muestra que el comportamiento asintótico de $M^k$ está gobernado exactamente por $ρ(M)$, no por $‖M‖$. Demostración completa: [[Teorema 4.24]].

## Autovalores de una matriz triangular

Los autovalores de una matriz triangular son exactamente los elementos de la diagonal ($\det(A-\lambda I)$ es producto de $(a_{ii}-\lambda)$, ver [[4.6 - Autovalores de matriz triangular|Problema 4.6]]).

## Ver también
- [[Metodos iterativos estacionarios - marco general]]
- [[Normas de vectores y matrices]]
