---
tags: [capitulo6, teorema]
---

# Teorema 6.12 (Contraction Mapping Theorem)

Sea $D\subset\mathbb{R}^n$ cerrado y $G:D\to D$ una contracción en $D$ con constante $0<\gamma<1$. Entonces existe un único punto fijo $x_*\in D$, es punto de atracción, y para todo $x_0\in D$ la iteración $x_{k+1}=G(x_k)$ está bien definida y cumple
$$
\|x_{k+1}-x_*\|\leq\gamma\|x_k-x_*\|,\qquad\|x_k-x_*\|\leq\frac{\gamma^k}{1-\gamma}\|x_1-x_0\|.
$$

> [!demostracion]-
> **Bien definida.** Como $G(D)\subset D$ y $x_0\in D$, todos los $x_k$ están en $D$.
>
> **Cauchy.** Por ser contracción,
> $$\|x_{i+1}-x_i\|=\|G(x_i)-G(x_{i-1})\|\leq\gamma\|x_i-x_{i-1}\|\leq\cdots\leq\gamma^i\|x_1-x_0\|.$$
> Para $k\in\mathbb{N}_0$, $p\in\mathbb{N}$, por la desigualdad triangular y la suma geométrica:
> $$\|x_{k+p}-x_k\|\leq\sum_{i=k}^{k+p-1}\|x_{i+1}-x_i\|\leq\Big(\sum_{i=k}^{k+p-1}\gamma^i\Big)\|x_1-x_0\|=\frac{\gamma^k-\gamma^{k+p}}{1-\gamma}\|x_1-x_0\|\leq\frac{\gamma^k}{1-\gamma}\|x_1-x_0\|. \tag{$*$}$$
> Como $\gamma<1$, $\{x_k\}$ es de Cauchy en $\mathbb{R}^n$ (completo), así que $x_k\to x_*$ para algún $x_*\in\mathbb{R}^n$.
>
> **$x_*\in D$ y es punto fijo.** $D$ es cerrado y $\{x_k\}\subset D$, luego $x_*\in D$. $G$ es continua (es Lipschitz), así que $G(x_*)=\lim G(x_k)=\lim x_{k+1}=x_*$.
>
> **Unicidad.** Si $y_*\in D$ fuera otro punto fijo, $\|y_*-x_*\|=\|G(y_*)-G(x_*)\|\leq\gamma\|y_*-x_*\|<\|y_*-x_*\|$, absurdo.
>
> **Convergencia lineal.** $\|x_{k+1}-x_*\|=\|G(x_k)-G(x_*)\|\leq\gamma\|x_k-x_*\|$.
>
> **Estimación del error.** Tomando $p\to\infty$ en $(*)$: $\|x_*-x_k\|\leq\frac{\gamma^k}{1-\gamma}\|x_1-x_0\|$.

## Ver también
- [[Contracciones y teorema de punto fijo]]
- [[6.3 - Contraccion en acotados]] (una variante: contracción en cada acotado)
- [[Demostraciones - Capitulo 6|Demostraciones]]
