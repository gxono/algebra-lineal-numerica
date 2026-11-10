---
tags: [capitulo5, teorema]
---

# Teorema 5.1

Sea $A$ sdp, $b\in\mathbb{R}^n$, $\phi(x)=\frac12x^TAx-b^Tx$. El mínimo de $\phi$ en $\mathbb{R}^n$ se alcanza **solo** en la solución $x_*$ de $Ax=b$, y vale $-\frac12b^Tx_*$.

> [!demostracion]-
> Con $x_*=A^{-1}b$ y $z=y-x_*$:
> $$\phi(y)=\phi(x_*+z)=\tfrac12(x_*+z)^TA(x_*+z)-b^T(x_*+z) = \phi(x_*)+\tfrac12z^TAz,$$
> (los términos lineales se cancelan usando $Ax_*=b$). Como $A$ sdp, $\tfrac12z^TAz>0$ para $z\neq0$, así que $\phi(y)>\phi(x_*)$ para $y\neq x_*$.

**Identidad clave (eq. 5.1):** para todo $y$,
$$
\phi(y) = \phi(x_*)+\tfrac12\|y-x_*\|_A^2, \qquad \|z\|_A:=\sqrt{z^TAz}.
$$

## Ver también
- [[Funcional de energia y caracterizacion]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
