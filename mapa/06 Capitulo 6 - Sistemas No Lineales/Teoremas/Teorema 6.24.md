---
tags: [capitulo6, teorema]
---

# Teorema 6.24 (convergencia cuadrática de Newton)

Sea $D\subset\mathbb{R}^n$ abierto, $F:D\to\mathbb{R}^n$ diferenciable en $D$, $F(x_*)=0$, y supongamos que existen $\delta,L>0$ con
$$
\|F'(x)-F'(x_*)\|\leq L\|x-x_*\|,\qquad\forall x\in B(x_*,\delta). \tag{6.5}
$$
Si $F'(x_*)$ es invertible, existe $c>0$ tal que si $x_0$ está suficientemente cerca de $x_*$,
$$
\|x_{k+1}-x_*\|\leq c\|x_k-x_*\|^2,\qquad k=0,1,2,\dots
$$

> [!demostracion]-
> **Paso 1: (6.5) implica**
> $$\|F(x)-F(x_*)-F'(x_*)(x-x_*)\|\leq\frac{L}{2}\|x-x_*\|^2,\qquad x\in B(x_*,\delta). \tag{6.6}$$
> La bola es convexa, así que por el [[Teorema 6.16]] $F(x)-F(x_*)=\int_0^1F'(x_*+t(x-x_*))(x-x_*)\,dt$. Restando $F'(x_*)(x-x_*)=\int_0^1F'(x_*)(x-x_*)\,dt$ y usando el [[Teorema 6.17]] y (6.5) (el punto $x_*+t(x-x_*)$ está a distancia $t\|x-x_*\|$ de $x_*$):
> $$\begin{aligned} \|F(x)-F(x_*)-F'(x_*)(x-x_*)\|&\leq\int_0^1\|F'(x_*+t(x-x_*))-F'(x_*)\|\,\|x-x_*\|\,dt\\ &\leq\int_0^1Lt\|x-x_*\|^2\,dt=\frac{L}{2}\|x-x_*\|^2. \end{aligned}$$
>
> **Paso 2: cota de $G(x)-x_*$.** Con $G(x)=x-F'(x)^{-1}F(x)$ y $F(x_*)=0$:
> $$\begin{aligned} \|G(x)-G(x_*)\|&=\big\|-F'(x)^{-1}\big[F(x)-F(x_*)-F'(x)(x-x_*)\big]\big\|\\ &=\big\|-F'(x)^{-1}\big[F(x)-F(x_*)-F'(x_*)(x-x_*)-(F'(x)-F'(x_*))(x-x_*)\big]\big\|\\ &\leq\|F'(x)^{-1}\|\Big[\tfrac{L}{2}\|x-x_*\|^2+L\|x-x_*\|^2\Big]=\tfrac{3L}{2}\|F'(x)^{-1}\|\,\|x-x_*\|^2. \end{aligned}$$
>
> **Paso 3: controlar $\|F'(x)^{-1}\|$.** $F'$ es continua en $x_*$ (por (6.5)) y $F'(x_*)$ invertible; por el [[Corolario 6.22]] existe $\delta'>0$ con $\|F'(x)^{-1}\|\leq2\|F'(x_*)^{-1}\|$ en $B(x_*,\delta')$. Entonces
> $$\|G(x)-x_*\|\leq\underbrace{3L\|F'(x_*)^{-1}\|}_{=:c}\|x-x_*\|^2,\qquad x\in B(x_*,\delta)\cap B(x_*,\delta').$$
>
> **Paso 4.** Por el [[Teorema 6.23]], $x_*$ es punto de atracción: existe $\delta''\leq\min\{\delta,\delta'\}$ tal que si $x_0\in B(x_*,\delta'')$ todos los $x_k$ quedan ahí. Luego $\|x_{k+1}-x_*\|=\|G(x_k)-x_*\|\leq c\|x_k-x_*\|^2$.

## Ver también
- [[Metodo de Newton - convergencia cuadratica]] (Observación 6.25, ejemplos donde falla (6.5))
- [[Demostraciones - Capitulo 6|Demostraciones]]
