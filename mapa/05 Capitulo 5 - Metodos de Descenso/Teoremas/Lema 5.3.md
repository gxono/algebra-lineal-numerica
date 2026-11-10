---
tags: [capitulo5, teorema]
---

# Lema 5.3 (lema de caracterización)

Sea $S\subset\mathbb{R}^n$ cualquier conjunto (en particular, un afín $x_0+\mathcal{K}$). Entonces, para $y\in S$:
$$
y \text{ minimiza } \phi \text{ sobre } S \iff y \text{ minimiza } \|\cdot-x_*\|_A \text{ sobre } S.
$$

> [!demostracion]-
> Por la identidad (eq. 5.1), para todo $y\in S$: $\phi(y)=\phi(x_*)+\tfrac12\|y-x_*\|_A^2$. Los dos términos de la derecha que dependen de $y$ difieren en un factor $\tfrac12$ y una constante aditiva $\phi(x_*)$ que **no depende de $y$** (solo de $x_*$, fijo). Minimizar $\phi(y)$ sobre $S$ es entonces equivalente a minimizar $\|y-x_*\|_A^2$ sobre $S$ (sumar una constante y multiplicar por un escalar positivo no cambia el minimizador), que a su vez es equivalente a minimizar $\|y-x_*\|_A$ (la raíz cuadrada es creciente).

> [!idea]
> Este lema es la herramienta de **todo** el capítulo: cada iteración que minimiza $\phi$ sobre un afín está, en realidad, buscando el punto de ese afín más cercano a $x_*$ en norma $A$. Así se traducen preguntas de optimización en preguntas de aproximación.

## Ver también
- [[Funcional de energia y caracterizacion]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
