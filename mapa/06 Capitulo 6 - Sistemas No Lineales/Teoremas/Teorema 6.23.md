---
tags: [capitulo6, teorema]
---

# Teorema 6.23 (buena definición de Newton y convergencia local superlineal)

Sea $D\subset\mathbb{R}^n$ abierto, $F:D\to\mathbb{R}^n$ diferenciable en $D$ y $x_*\in D$ con $F(x_*)=0$. Si $F'$ es continua en $x_*$ y $F'(x_*)$ es invertible:
1. **Buena definición:** $G(x):=x-F'(x)^{-1}F(x)$ está bien definida en un entorno de $x_*$.
2. **Convergencia local:** $x_{k+1}=x_k-F'(x_k)^{-1}F(x_k)$ converge localmente a $x_*$.
3. **Superlineal:** si $x_0$ está suficientemente cerca de $x_*$, $\displaystyle\lim_{k\to\infty}\frac{\|x_{k+1}-x_*\|}{\|x_k-x_*\|}=0$.

> [!demostracion]-
> **1.** Es el [[Corolario 6.22]] con $C=F'$: existe $\delta>0$ con $F'(x)$ invertible en $B(x_*,\delta)$.
>
> **2.** Como $F(x_*)=0$, $G(x_*)=x_*$. Por el [[6.8 - Derivada de G en Newton es nula|Problema 6.8]] ([[Derivada de la iteracion de Newton en la raiz]]), $G$ es diferenciable en $x_*$ con $G'(x_*)=0$. Luego $\rho(G'(x_*))=0<1$ y por [[Teorema 6.8|Ostrowski]] $x_*$ es punto de atracción.
>
> **3.** $G$ diferenciable en $x_*$ con $G(x_*)=x_*$ y $G'(x_*)=0$ significa
> $$\lim_{x\to x_*}\frac{\|G(x)-x_*\|}{\|x-x_*\|}=0.$$
> Tomando $x_0$ en la bola de atracción, $x_k\to x_*$, y evaluando el límite anterior a lo largo de la sucesión ($x_{k+1}=G(x_k)$):
> $$\lim_{k\to\infty}\frac{\|x_{k+1}-x_*\|}{\|x_k-x_*\|}=0.$$
> (Si algún $x_k=x_*$ la sucesión es constante desde ahí y no hay nada que probar.)

## Ver también
- [[Metodo de Newton - definicion y buena definicion]]
- [[Teorema 6.24]] (con más regularidad: cuadrática)
- [[Demostraciones - Capitulo 6|Demostraciones]]
