---
tags: [capitulo4, teorema]
---

# Teorema 4.24 (fórmula de Gelfand)

Sea $A\in\mathbb R^{n\times n}$. Entonces, cualquiera sea la norma inducida $\|\cdot\|$,
$$
\rho(A) = \lim_{k\to\infty}\|A^k\|^{1/k}.
$$

> [!demostracion]-
> **Paso 1: $\rho(A)\leq\|A^k\|^{1/k}$ para todo $k$.** Sea $\lambda\in\sigma(A)$ con autovector $x\neq0$. Entonces $Ax=\lambda x\implies A^kx=\lambda^kx$, así que
> $$\|A^kx\|=|\lambda|^k\|x\| \quad\Longrightarrow\quad |\lambda|^k = \frac{\|A^kx\|}{\|x\|} \leq \|A^k\|$$
> (usando $\|A^kx\|\leq\|A^k\|\|x\|$ por consistencia). Así $|\lambda|\leq\|A^k\|^{1/k}$, y como $\lambda$ era cualquier autovalor, $\rho(A)\leq\|A^k\|^{1/k}$ para todo $k$ - en particular, $\rho(A)\leq\liminf_k\|A^k\|^{1/k}$.
>
> **Paso 2: el límite de $\|A^k\|^{1/k}$ es exactamente $\rho(A)$.** Sea $\varepsilon>0$ y definamos
> $$\tilde A := \frac{1}{\rho(A)+\varepsilon}A \quad\Longrightarrow\quad \rho(\tilde A)=\frac{\rho(A)}{\rho(A)+\varepsilon}<1$$
> (los autovalores de $\tilde A$ son los de $A$ escalados por $\frac{1}{\rho(A)+\varepsilon}$). Por el [[Teorema 4.23]] (equivalencia (1)$\Rightarrow$(3)), $\lim_{k\to\infty}\|\tilde A^k\|=0$. Por lo tanto, existe $k_0\in\mathbb N$ tal que
> $$\|\tilde A^k\| = \frac{\|A^k\|}{(\rho(A)+\varepsilon)^k} < 1, \qquad \forall k\geq k_0,$$
> es decir, $\|A^k\|^{1/k} < \rho(A)+\varepsilon$ para todo $k\geq k_0$. Combinando con el Paso 1 ($\rho(A)\leq\|A^k\|^{1/k}$ para todo $k$):
> $$\rho(A) \leq \|A^k\|^{1/k} < \rho(A)+\varepsilon, \qquad \forall k\geq k_0.$$
> Como $\varepsilon>0$ era arbitrario, esto es exactamente $\lim_{k\to\infty}\|A^k\|^{1/k}=\rho(A)$.

> [!idea]
> Esta fórmula es la justificación rigurosa de una afirmación que se usa informalmente en todo el capítulo: "el comportamiento asintótico de $A^k$ está gobernado por $\rho(A)$, no por una norma particular de $A$". Es más fuerte que el [[Teorema 4.23]]: no solo dice *si* $A^k\to0$, sino a qué *tasa* (asintóticamente, $\|A^k\|\approx\rho(A)^k$) - es la base de por qué, en la práctica, se usa $\rho(M)$ (y no $\|M\|$) como la medida exacta de qué tan rápido converge un método iterativo estacionario.

## Ver también
- [[Teorema 4.23]] (usado en el Paso 2)
- [[Espectro y radio espectral]] (nota de tema completa, donde se cita esta fórmula)
- [[Demostraciones - Capitulo 4|Demostraciones]]
