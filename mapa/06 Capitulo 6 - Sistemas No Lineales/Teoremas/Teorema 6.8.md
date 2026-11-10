---
tags: [capitulo6, teorema]
---

# Teorema 6.8 (Ostrowski)

Sea $G:\mathbb{R}^n\to\mathbb{R}^n$ diferenciable en $x_*$, con $x_*=G(x_*)$. Si $\rho(G'(x_*))<1$, entonces $x_*$ es un punto de atracción de $G$.

> [!demostracion]-
> Sea $\sigma:=\rho(G'(x_*))<1$ y fijemos $\varepsilon>0$ con $\sigma+2\varepsilon<1$. Por el [[Teorema 4.20]] existe una norma inducida con
> $$\|G'(x_*)\|<\sigma+\varepsilon.$$
> Como $G$ es diferenciable en $x_*$, existe $\delta>0$ tal que
> $$\|G(x)-G(x_*)-G'(x_*)(x-x_*)\|<\varepsilon\|x-x_*\|\qquad\text{si }x\in B(x_*,\delta).$$
> Sea $x_0\in B(x_*,\delta)$ y $x_{k+1}=G(x_k)$. Si $x_k\in B(x_*,\delta)$:
> $$\begin{aligned} \|x_{k+1}-x_*\|=\|G(x_k)-G(x_*)\|&\leq\|G(x_k)-G(x_*)-G'(x_*)(x_k-x_*)\|+\|G'(x_*)(x_k-x_*)\|\\ &<\varepsilon\|x_k-x_*\|+(\sigma+\varepsilon)\|x_k-x_*\|=(\sigma+2\varepsilon)\|x_k-x_*\|<\|x_k-x_*\|, \end{aligned}$$
> así que $x_{k+1}\in B(x_*,\delta)$. Por inducción todos los iterados quedan en la bola, y además
> $$\|x_{k+1}-x_*\|\leq(\sigma+2\varepsilon)^{k+1}\|x_0-x_*\|\to0,$$
> porque $\sigma+2\varepsilon<1$. Luego $x_k\to x_*$.

> [!idea]
> La norma "buena" la da el Teorema 4.20, y puede ser una norma rara (no $\|\cdot\|_1,\|\cdot\|_2,\|\cdot\|_\infty$). Pero como en $\mathbb{R}^n$ todas las normas son equivalentes, la **convergencia** (no la monotonía del error) vale en cualquier norma.

## Ver también
- [[Iteracion de punto fijo y teorema de Ostrowski]]
- [[Teorema 4.20]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
