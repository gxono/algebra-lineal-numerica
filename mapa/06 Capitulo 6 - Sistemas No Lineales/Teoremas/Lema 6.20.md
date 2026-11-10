---
tags: [capitulo6, lema]
---

# Lema 6.20 (perturbación de la identidad)

Sea $M\in\mathbb{R}^{n\times n}$ con $\|M\|<1$ para alguna norma inducida. Entonces $I-M$ es invertible y
$$
\|(I-M)^{-1}\|\leq\frac{1}{1-\|M\|}.
$$

> [!demostracion]- (serie de Neumann)
> **Invertible.** Si $(I-M)x=0$ con $x\neq0$, entonces $Mx=x$ y $\|x\|=\|Mx\|\leq\|M\|\|x\|<\|x\|$, absurdo. Luego $I-M$ es invertible. (Equivalente: $1\notin\sigma(M)$ porque $\rho(M)\leq\|M\|<1$.)
>
> **Cota.** Las sumas parciales $S_N=\sum_{j=0}^{N}M^j$ forman una sucesión de Cauchy, pues para $N>m$
> $$\|S_N-S_m\|\leq\sum_{j=m+1}^N\|M\|^j\leq\frac{\|M\|^{m+1}}{1-\|M\|}\to0,$$
> así que convergen a cierta $S$. Además
> $$(I-M)S_N=S_N-MS_N=I-M^{N+1}\longrightarrow I\quad(\|M^{N+1}\|\leq\|M\|^{N+1}\to0),$$
> luego $(I-M)S=I$, es decir $S=(I-M)^{-1}=\sum_{j\geq0}M^j$. Finalmente,
> $$\|(I-M)^{-1}\|\leq\sum_{j=0}^\infty\|M\|^j=\frac{1}{1-\|M\|}.$$
>
> *(Alternativa sin series: de $(I-M)^{-1}=I+M(I-M)^{-1}$ sale $\|(I-M)^{-1}\|\leq1+\|M\|\|(I-M)^{-1}\|$, y despejando, la cota.)*

El apunte remite a la demostración del Teorema 4.25 (convergencia de métodos iterativos estacionarios), que usa exactamente esta serie.

## Ver también
- [[Lema 6.21]] (la generalización a $A+E$)
- [[Metodo de Newton - definicion y buena definicion]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
