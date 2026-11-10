---
tags: [capitulo6, lema]
---

# Lema 6.21 (perturbación de una matriz invertible)

Sea $A\in\mathbb{R}^{n\times n}$ invertible y $E\in\mathbb{R}^{n\times n}$ con $\|E\|<\dfrac{1}{\|A^{-1}\|}$. Entonces $A+E$ es invertible y
$$
\|(A+E)^{-1}\|\leq\frac{\|A^{-1}\|}{1-\|A^{-1}\|\|E\|}.
$$

> [!demostracion]-
> $A+E=A(I+A^{-1}E)$. Sea $M:=-A^{-1}E$; entonces $\|M\|\leq\|A^{-1}\|\|E\|<1$. Por el [[Lema 6.20]], $I-M$ es invertible con $\|(I-M)^{-1}\|\leq\frac{1}{1-\|M\|}$. Luego $A+E=A(I-M)$ es invertible (producto de invertibles) y
> $$\|(A+E)^{-1}\|=\|(I-M)^{-1}A^{-1}\|\leq\frac{\|A^{-1}\|}{1-\|M\|}\leq\frac{\|A^{-1}\|}{1-\|A^{-1}\|\|E\|},$$
> usando en el último paso que $\|M\|\leq\|A^{-1}\|\|E\|$ y que $t\mapsto\frac{1}{1-t}$ es creciente en $[0,1)$.

> [!idea]
> Dividiendo por $\|A\|$: si $\frac{\|E\|}{\|A\|}<\frac{1}{\kappa(A)}$, $A+E$ sigue siendo invertible. Las matrices mal condicionadas están **cerca** de ser singulares (ver [[Numero de condicion y error relativo]]).

## Ver también
- [[Corolario 6.22]]
- [[Metodo de Newton - definicion y buena definicion]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
