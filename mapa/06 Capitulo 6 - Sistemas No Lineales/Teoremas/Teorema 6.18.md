---
tags: [capitulo6, teorema]
---

# Teorema 6.18 (desigualdad del valor medio)

Sea $D\subset\mathbb{R}^n$ convexo y $G:D\to\mathbb{R}^n$ continuamente diferenciable en $D$. Entonces
$$
\|G(y)-G(x)\|\leq\max_{t\in[0,1]}\|G'(x+t(y-x))\|\,\|y-x\|,\qquad\forall x,y\in D,
$$
donde $\|\cdot\|$ es una norma de $\mathbb{R}^n$ y también la norma matricial inducida por ella.

> [!demostracion]-
> Por el [[Teorema 6.16]], $G(y)-G(x)=\int_0^1G'(x+t(y-x))(y-x)\,dt$. El integrando es continuo en $t$, así que por el [[Teorema 6.17]] y la consistencia de la norma inducida ($\|Bv\|\leq\|B\|\|v\|$):
> $$\begin{aligned} \|G(y)-G(x)\|&\leq\int_0^1\|G'(x+t(y-x))(y-x)\|\,dt\leq\int_0^1\|G'(x+t(y-x))\|\,\|y-x\|\,dt\\ &\leq\max_{t\in[0,1]}\|G'(x+t(y-x))\|\,\|y-x\|. \end{aligned}$$
> (El máximo existe porque $t\mapsto\|G'(x+t(y-x))\|$ es continua en el compacto $[0,1]$.)

## Ver también
- [[Corolario 6.19]] (la consecuencia práctica)
- [[Criterio de contraccion via la derivada]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
