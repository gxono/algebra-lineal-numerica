---
tags: [capitulo6, corolario]
---

# Corolario 6.19

Sea $D\subset\mathbb{R}^n$ convexo y $G:D\to\mathbb{R}^n$ continuamente diferenciable en $D$. Si existe $\alpha>0$ tal que $\|G'(x)\|\leq\alpha$ para todo $x\in D$, entonces
$$
\|G(x)-G(y)\|\leq\alpha\|x-y\|,\qquad\forall x,y\in D.
$$
En particular, si $\alpha<1$, $G$ es una contracción en $D$.

> [!demostracion]-
> Dados $x,y\in D$, todos los puntos $x+t(y-x)$, $t\in[0,1]$, están en $D$ (convexidad), así que $\|G'(x+t(y-x))\|\leq\alpha$. Por el [[Teorema 6.18]],
> $$\|G(y)-G(x)\|\leq\max_{t\in[0,1]}\|G'(x+t(y-x))\|\,\|y-x\|\leq\alpha\|y-x\|.$$

## Ver también
- [[Criterio de contraccion via la derivada#Receta para los ejercicios|Receta para usarlo en los ejercicios]]
- [[6.4 - Region para el teorema de contraccion]], [[6.5 - Contraccion para h chico]], [[6.2 - Contraccion local si rho menor que 1]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
