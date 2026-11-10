---
tags: [capitulo6, teorema]
---

# Teorema 6.29 (Newton-Baluev: convergencia global de Newton)

Sea $F:\mathbb{R}^n\to\mathbb{R}^n$ convexa y continuamente diferenciable, con $F'(x)$ invertible y $F'(x)^{-1}\geq0$ (entradas no negativas) para todo $x\in\mathbb{R}^n$. Si existe $x_*$ con $F(x_*)=0$, este cero es único y Newton converge a $x_*$ desde cualquier $x_0\in\mathbb{R}^n$. Además, componente a componente,
$$
x_*^{(i)}\leq\cdots\leq x_{k+1}^{(i)}\leq x_k^{(i)}\leq\cdots\leq x_1^{(i)},\qquad i=1,\dots,n.
$$

Notación: $v\geq0$ significa $v^{(i)}\geq0$ para todo $i$ (y $v\geq w$ significa $v-w\geq0$). Si $B\geq0$ (entradas) y $v\geq0$, entonces $Bv\geq0$.

> [!demostracion]-
> Sea $x_0$ arbitrario y $x_{k+1}=x_k-F'(x_k)^{-1}F(x_k)$, es decir, $F'(x_k)(x_{k+1}-x_k)=-F(x_k)$. $\quad(6.10)$
>
> **1. $F(x_k)\geq0$ para $k\geq1$.** Cada $f_i$ es convexa; por la [[Proposicion 6.28]] y (6.10),
> $$f_i(x_{k+1})-f_i(x_k)\geq\nabla f_i(x_k)^T(x_{k+1}-x_k)=-f_i(x_k)\implies f_i(x_{k+1})\geq0.$$
>
> **2. Monotonía.** Para $k\geq1$, $x_k-x_{k+1}=F'(x_k)^{-1}F(x_k)\geq0$ (producto de una matriz $\geq0$ por un vector $\geq0$). Así $x_{k+1}\leq x_k$.
>
> **3. Cota inferior $x_k\geq x_*$ para $k\geq1$.** Por la Proposición 6.28 en $x_k$ y $x_*$:
> $$\underbrace{f_i(x_*)}_{=0}-f_i(x_k)\geq\nabla f_i(x_k)^T(x_*-x_k)\iff -F(x_k)-F'(x_k)(x_*-x_k)\geq0.$$
> Multiplicando por $F'(x_k)^{-1}\geq0$: $-F'(x_k)^{-1}F(x_k)+(x_k-x_*)\geq0$, o sea
> $$x_k-x_*\geq F'(x_k)^{-1}F(x_k)\geq0\qquad(k\geq1\text{, por el paso 1}).$$
>
> **4. Convergencia.** Cada sucesión $\{x_k^{(i)}\}_{k\geq1}$ es decreciente y acotada inferiormente por $x_*^{(i)}$, así que converge: $x_k\to y$. Por continuidad de $F$ y $F'$ y (6.10):
> $$F(y)=\lim F(x_k)=\lim F'(x_k)(x_k-x_{k+1})=F'(y)(y-y)=0.$$
>
> **5. Unicidad del cero** (el apunte lo deja como ejercicio). Sean $x_*$, $y_*$ ceros de $F$. Por la Proposición 6.28 aplicada a cada $f_i$ en $x_*$:
> $$0=f_i(y_*)-f_i(x_*)\geq\nabla f_i(x_*)^T(y_*-x_*)\implies F'(x_*)(y_*-x_*)\leq0.$$
> Multiplicando por $F'(x_*)^{-1}\geq0$: $y_*-x_*\leq0$. Intercambiando los roles de $x_*$ e $y_*$: $x_*-y_*\leq0$. Luego $x_*=y_*$.
>
> En particular el límite $y$ del paso 4 es **el** cero $x_*$, y las desigualdades de los pasos 2–3 dan la cadena monótona del enunciado.

## Ver también
- [[Convergencia global de Newton - Newton-Baluev]]
- [[6.18 - Newton-Baluev sin convexidad]] (qué pasa si falta la convexidad)
- [[Demostraciones - Capitulo 6|Demostraciones]]
