---
tags: [capitulo5, teorema, definicion]
---

# El funcional de energía y la caracterización de $Ax=b$

> [!idea]
> Para $A$ sdp, resolver $Ax=b$ es **exactamente lo mismo** que minimizar $\phi(x)=\frac12x^TAx-b^Tx$ en $\mathbb{R}^n$. Esto convierte todo el capítulo en un problema de optimización, y $\|y-x_*\|_A^2$ (la "distancia" en la norma de energía) mide directamente qué tan lejos está $y$ del óptimo.

## Teorema 5.1

Sea $A$ sdp, $b\in\mathbb{R}^n$, $\phi(x)=\frac12x^TAx-b^Tx$. El mínimo de $\phi$ en $\mathbb{R}^n$ se alcanza **solo** en la solución $x_*$ de $Ax=b$, y vale $-\frac12b^Tx_*$. Demostración completa: [[Teorema 5.1]].

**Identidad clave (eq. 5.1):** para todo $y$,
$$
\phi(y) = \phi(x_*)+\tfrac12\|y-x_*\|_A^2, \qquad \|z\|_A:=\sqrt{z^TAz}.
$$
($\|\cdot\|_A$ es una norma de verdad, no solo una expresión parecida a una - demostración: [[Producto escalar A y norma A]].)

## Lema 5.3 (lema de caracterización)

Sea $S\subset\mathbb{R}^n$ cualquier conjunto (en particular, un afín $x_0+\mathcal{K}$). Entonces, para $y\in S$:
$$
y \text{ minimiza } \phi \text{ sobre } S \iff y \text{ minimiza } \|\cdot-x_*\|_A \text{ sobre } S.
$$

Consecuencia directa de la identidad (eq. 5.1): minimizar $\phi(y)$ y minimizar $\|y-x_*\|_A$ sobre el mismo conjunto son lo mismo, porque difieren solo en un factor y una constante que no dependen de $y$. Demostración completa: [[Lema 5.3]].

> [!idea]
> Este lema es la herramienta de **todo** el capítulo: cada iteración que minimiza $\phi$ sobre un afín está, en realidad, buscando el punto de ese afín más cercano a $x_*$ en norma $A$. Así se traducen preguntas de optimización en preguntas de aproximación.

## Ver también
- [[Descenso mas pronunciado]]
- [[Gradiente conjugado - definicion y propiedades]]
