---
tags: [glosario, capitulo6]
aliases: [contracción, Lipschitz]
---

# Contracción

## Definición

$G:D\subset\mathbb{R}^n\to\mathbb{R}^n$ es una **contracción** en $D$ (con la norma $\|\cdot\|$) si existe $\gamma\in(0,1)$ con $\|G(x)-G(y)\|\leq\gamma\|x-y\|$ para todo $x,y\in D$ (Lipschitz con constante menor que 1). Si además $D$ es cerrado y $G(D)\subset D$, $G$ tiene un único punto fijo en $D$ y la iteración converge a él desde cualquier $x_0\in D$ (Teorema 6.12). Se verifica en la práctica acotando $\|G'\|<1$ en un convexo (Corolario 6.19).

## Ver también
- [[Contracciones y teorema de punto fijo]]
- [[Criterio de contraccion via la derivada]]
- [[Glosario - Indice|Glosario]]
