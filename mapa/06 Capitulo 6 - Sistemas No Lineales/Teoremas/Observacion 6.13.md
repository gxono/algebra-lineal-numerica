---
tags: [capitulo6, observacion]
---

# Observación 6.13 (contracción lineal $\iff\rho(M)<1$)

Si $D=\mathbb{R}^n$ y $G(x)=Mx+c$, entonces $G$ es una contracción en alguna norma $\iff\rho(M)<1$.

> [!demostracion]-
> **($\Rightarrow$)** Si $\|G(x)-G(y)\|\leq\gamma\|x-y\|$ con $\gamma<1$ para todo $x,y$, entonces, usando la norma inducida y el [[Radio espectral acotado por norma inducida|radio espectral acotado por la norma]]:
> $$\rho(M)\leq\|M\|=\sup_{z\neq0}\frac{\|Mz\|}{\|z\|}=\sup_{z\neq0}\frac{\|G(z)-G(0)\|}{\|z\|}\leq\gamma<1.$$
> **($\Leftarrow$)** Si $\rho(M)<1$, por el [[Teorema 4.20]] existe una norma inducida con $\|M\|<1$, y
> $$\|G(x)-G(y)\|=\|M(x-y)\|\leq\|M\|\,\|x-y\|,$$
> así que $G$ es contracción en esa norma con $\gamma=\|M\|$.

Por eso el [[Teorema 6.12]] generaliza al caso no lineal el [[Teorema de convergencia de metodos iterativos estacionarios|teorema de convergencia del Capítulo 4]].

## Ver también
- [[Contracciones y teorema de punto fijo]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
