---
tags: [capitulo5, teorema]
---

# Teorema 5.42 (si $A$ está cerca de la identidad)

Si $\|I-A\|_2\leq\rho<1$: $\|r_k\|_2\leq\rho^k\|r_0\|_2$. Combinado con el [[Numero de condicion y error relativo|Lema 4.10]]: $\|e_k\|_2\leq\kappa_2(A)\rho^k\|e_0\|_2$.

> [!demostracion]-
> Sea $p_k(z):=(1-z)^k$. Entonces $p_k\in\overline{\mathcal P}_k$ (grado $k$, $p_k(0)=1$) y $p_k(A)=(I-A)^k$. Por el Corolario 5.37 (del [[Teorema 5.36]]):
> $$\|r_k\|_2 \leq \|p_k(A)\|_2\|r_0\|_2 = \|(I-A)^k\|_2\|r_0\|_2 \leq \|I-A\|_2^k\|r_0\|_2 \leq \rho^k\|r_0\|_2,$$
> usando submultiplicatividad de $\|\cdot\|_2$ en el último paso. La segunda desigualdad se obtiene aplicando el Lema 4.10 ($\|e_k\|_2/\|x_*\|_2\leq\kappa_2(A)\|r_k\|_2/\|b\|_2$, reescrito con $\|e_0\|_2$ en vez de $\|x_*\|_2$) a la cota recién obtenida.

## Ver también
- [[GMRes - precondicionamiento]] (nota de tema completa)
- [[Teorema 5.36]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
