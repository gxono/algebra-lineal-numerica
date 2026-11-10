---
tags: [capitulo6, moc]
---

# Formulario - Capítulo 6 (Sistemas no lineales)

## Jacobiano

$$
F'(x)=\Big(\frac{\partial f_i}{\partial x_j}(x)\Big)_{ij},\qquad\lim_{x\to x_0}\frac{\|F(x)-F(x_0)-F'(x_0)(x-x_0)\|}{\|x-x_0\|}=0.
$$

## [[Iteracion de punto fijo y teorema de Ostrowski|Punto fijo]]

$$
x_{k+1}=G(x_k).\qquad\textbf{Ostrowski: }\rho(G'(x_*))<1\implies x_*\text{ atrae (local, solo suficiente).}
$$

## [[Contracciones y teorema de punto fijo|Contracción]] (Teorema 6.12)

$D$ cerrado, $G(D)\subset D$, $\|G(x)-G(y)\|\leq\gamma\|x-y\|$, $\gamma<1$ ⟹ único $x_*\in D$ y
$$
\|x_{k+1}-x_*\|\leq\gamma\|x_k-x_*\|,\qquad\|x_k-x_*\|\leq\frac{\gamma^k}{1-\gamma}\|x_1-x_0\|,\qquad k\geq\frac{\log\big(\text{tol}(1-\gamma)/\|x_1-x_0\|\big)}{\log\gamma}.
$$

## [[Criterio de contraccion via la derivada|Criterio vía la derivada]]

$D$ convexo:
$$
G(y)-G(x)=\int_0^1G'(x+t(y-x))(y-x)\,dt,\qquad\|G(y)-G(x)\|\leq\max_{t\in[0,1]}\|G'(x+t(y-x))\|\,\|y-x\|.
$$
$\|G'(x)\|\leq\alpha<1$ en $D$ ⟹ contracción. En la práctica: $\|G'\|_\infty=\max_i\sum_j|\partial_jg_i|$.

Máximos de una variable que aparecen seguido: $\max\frac{|u|}{(1+u^2)^2}=\frac{9}{16\sqrt3}$ (en $u=\frac1{\sqrt3}$); $\max|xe^{-x^2}|=\frac{1}{\sqrt{2e}}$; $\max\frac{|x|}{1+x^4}=\frac{3^{3/4}}{4}$; $\max\frac{|x|}{1+x^2}=\frac12$.

## Lemas de perturbación

$$
\|M\|<1\implies\|(I-M)^{-1}\|\leq\frac{1}{1-\|M\|};\qquad\|E\|<\frac{1}{\|A^{-1}\|}\implies\|(A+E)^{-1}\|\leq\frac{\|A^{-1}\|}{1-\|A^{-1}\|\|E\|}.
$$
Cerca de $x_*$ (con $F'$ continua, $F'(x_*)$ invertible): $\|F'(x)^{-1}\|\leq2\|F'(x_*)^{-1}\|$.

## [[Metodo de Newton - definicion y buena definicion|Newton]]

$$
F'(x_k)\,h=F(x_k),\qquad x_{k+1}=x_k-h;\qquad G(x)=x-F'(x)^{-1}F(x),\ G'(x_*)=0.
$$
- $F'$ continua en $x_*$, $F'(x_*)$ invertible ⟹ local, **superlineal**.
- $+\ \|F'(x)-F'(x_*)\|\leq L\|x-x_*\|$ ⟹ $\|x_{k+1}-x_*\|\leq c\|x_k-x_*\|^2$, con $c=3L\|F'(x_*)^{-1}\|$.
- Cota clave: $\|F(x)-F(x_*)-F'(x_*)(x-x_*)\|\leq\frac L2\|x-x_*\|^2$.
- **Simplificado** $x_{k+1}=x_k-F'(x_0)^{-1}F(x_k)$: lineal, factor $\approx\rho\big(I-F'(x_0)^{-1}F'(x_*)\big)$.
- 1D, raíz de multiplicidad $m$: lineal con tasa $1-\frac1m$.

## [[Convergencia global de Newton - Newton-Baluev|Newton-Baluev]]

$F$ convexa, $\mathcal C^1$, $F'(x)^{-1}\geq0$ para todo $x$, existe un cero ⟹ único, convergencia **global** y $x_*\leq\cdots\leq x_{k+1}\leq x_k\leq\cdots\leq x_1$. En 1D: $F$ convexa con $F'>0$.

Convexa: $f(y)-f(x)\geq\nabla f(x)^T(y-x)$.

## Orden de convergencia

$$
p\approx\frac{\log(e_{k+1}/e_k)}{\log(e_k/e_{k-1})}.
$$
