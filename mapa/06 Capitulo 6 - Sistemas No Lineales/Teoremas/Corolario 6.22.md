---
tags: [capitulo6, corolario]
---

# Corolario 6.22

Sea $C:D\subset\mathbb{R}^n\to\mathbb{R}^{n\times n}$ continua en $x_*\in D$ con $C(x_*)$ invertible. Entonces existe $\delta>0$ tal que para todo $x\in B(x_*,\delta)$, $C(x)$ es invertible y
$$
\|C(x)^{-1}\|\leq2\|C(x_*)^{-1}\|.
$$

> [!demostracion]-
> Escribimos $C(x)=\underbrace{C(x_*)}_{=:A}+\underbrace{(C(x)-C(x_*))}_{=:E}$. Por continuidad de $C$ en $x_*$, existe $\delta>0$ con
> $$\|C(x)-C(x_*)\|<\frac{1}{2\|C(x_*)^{-1}\|}\qquad\text{si }x\in B(x_*,\delta).$$
> En particular $\|E\|<1/\|A^{-1}\|$, y por el [[Lema 6.21]] $C(x)$ es invertible con
> $$\|C(x)^{-1}\|\leq\frac{\|C(x_*)^{-1}\|}{1-\|C(x_*)^{-1}\|\,\|C(x)-C(x_*)\|}<\frac{\|C(x_*)^{-1}\|}{1-\|C(x_*)^{-1}\|\frac{1}{2\|C(x_*)^{-1}\|}}=2\|C(x_*)^{-1}\|.$$

Se usa con $C=F'$: el jacobiano sigue siendo invertible (y con inversa controlada) cerca de una raíz donde lo es.

## Ver también
- [[Teorema 6.23]], [[Teorema 6.24]], [[Derivada de la iteracion de Newton en la raiz]], [[Convergencia local de Newton simplificado]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
