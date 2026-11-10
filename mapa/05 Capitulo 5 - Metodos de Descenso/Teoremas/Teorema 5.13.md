---
tags: [capitulo5, teorema]
---

# Teorema 5.13

$$
\|x_k-x_*\|_A = \min_{p\in\overline{\mathcal P}_k}\|p(A)(x_0-x_*)\|_A, \qquad \overline{\mathcal P}_k:=\{p\in\mathcal P_k:p(0)=1\}.
$$

> [!demostracion]-
> Por el [[Funcional de energia y caracterizacion|Lema 5.3]], $x_k$ minimiza $\phi$ sobre $x_0+\mathcal K_k$ $\iff$ $x_k$ minimiza $\|\cdot-x_*\|_A$ sobre $x_0+\mathcal K_k$. Basta entonces reescribir ese conjunto. Si $w\in x_0+\mathcal K_k$, $w=x_0+\sum_{i=0}^{k-1}\gamma_iA^ir_0$ para ciertos $\gamma_0,\dots,\gamma_{k-1}\in\mathbb{R}$. Usando $r_0=b-Ax_0=A(x_*-x_0)$:
> $$w-x_* = (x_0-x_*)+\sum_{i=0}^{k-1}\gamma_iA^i\big(A(x_*-x_0)\big) = \Big(I-\sum_{i=0}^{k-1}\gamma_iA^{i+1}\Big)(x_0-x_*) = p(A)(x_0-x_*),$$
> con $p(z):=1-\sum_{i=0}^{k-1}\gamma_iz^{i+1}\in\overline{\mathcal P}_k$ (grado $\leq k$, $p(0)=1$). Recíprocamente, todo $p\in\overline{\mathcal P}_k$ se escribe así para algún $w\in x_0+\mathcal K_k$. Luego minimizar $\|w-x_*\|_A$ sobre $w\in x_0+\mathcal K_k$ es lo mismo que minimizar $\|p(A)(x_0-x_*)\|_A$ sobre $p\in\overline{\mathcal P}_k$, y como $x_k$ es ese mínimo, se obtiene la igualdad.

## Ver también
- [[Cotas de error del gradiente conjugado]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
