---
tags: [capitulo6, algoritmo, teorema, lema, corolario]
---

# Método de Newton - definición, buena definición y convergencia local

> [!idea]
> Reemplazar $F$ por su aproximación lineal en $x_k$ y anularla: $F'(x_k)h=F(x_k)$, $x_{k+1}=x_k-h$. **Nunca** se invierte $F'(x_k)$, se resuelve un sistema lineal por paso (con lo de los Capítulos 2–5). Si $F'$ es continua en $x_*$ y $F'(x_*)$ invertible, Newton está bien definido cerca de $x_*$ y converge localmente de forma **superlineal** (Teorema 6.23), porque es un punto fijo con $G'(x_*)=0$.

## La iteración (6.4)

Aproximación afín de $F$ en $x_k$: $P(x)=F(x_k)+F'(x_k)(x-x_k)$. Se define $x_{k+1}$ por $P(x_{k+1})=0$:
$$
\begin{cases}
\text{hallar } h\in\mathbb{R}^n \text{ tal que } F'(x_k)\,h=F(x_k),\\
x_{k+1}=x_k-h.
\end{cases}\tag{6.4}
$$
Es una iteración de punto fijo con $G(x)=x-F'(x)^{-1}F(x)$.

> [!algoritmo] Newton (versión del [[6.9 - Programar Newton|Problema 6.9]])
> ```julia
> function newton(F, DF, x0, miter, tol)
>     X = zeros(length(x0), miter + 1)
>     X[:, 1] = x0
>     x = copy(x0)
>     for k in 1:miter
>         h = DF(x) \ F(x)        # sistema lineal, NUNCA inv(DF(x))
>         x = x - h
>         X[:, k + 1] = x
>         norm(h) < tol && return x, k, X[:, 1:k+1]
>     end
>     @warn "Newton: no convergió en $miter iteraciones"
>     return x, miter, X
> end
> ```
> Costo por paso: evaluar $F$ y $F'$ ($n^2$ derivadas) + resolver un sistema $n\times n$ ($O(n^3)$ denso con LU, menos si $F'$ es rala/tridiagonal, como en el [[6.17 - Newton en problema discretizado|Problema 6.17]]).

## Lemas de perturbación

**Lema 6.20 (perturbación de la identidad).** Si $\|M\|<1$ (norma inducida), $I-M$ es invertible y $\|(I-M)^{-1}\|\leq\dfrac{1}{1-\|M\|}$. Ver [[Lema 6.20]] (serie de Neumann).

**Lema 6.21 (perturbación de una invertible).** $A$ invertible, $\|E\|<1/\|A^{-1}\|$ $\implies$ $A+E$ invertible y
$$
\|(A+E)^{-1}\|\leq\frac{\|A^{-1}\|}{1-\|A^{-1}\|\|E\|}.
$$
Ver [[Lema 6.21]]. Lectura: las matrices invertibles forman un **abierto**; si te movés menos que $1/\|A^{-1}\|$, seguís siendo invertible. (Relación con el [[Numero de condicion y error relativo|número de condición]]: la distancia relativa a la singularidad es $\geq 1/\kappa(A)$.)

**Corolario 6.22.** $C:D\to\mathbb{R}^{n\times n}$ continua en $x_*$ con $C(x_*)$ invertible $\implies$ existe $\delta>0$ con $C(x)$ invertible y $\|C(x)^{-1}\|\leq2\|C(x_*)^{-1}\|$ en $B(x_*,\delta)$. Ver [[Corolario 6.22]].

## Teorema 6.23 (buena definición y primer resultado de convergencia local)

> $D$ abierto, $F$ diferenciable en $D$, $F(x_*)=0$, $F'$ continua en $x_*$ y $F'(x_*)$ invertible. Entonces:
> 1. $G(x)=x-F'(x)^{-1}F(x)$ está **bien definida** en un entorno de $x_*$ (Corolario 6.22 con $C=F'$).
> 2. Newton **converge localmente** a $x_*$ (Ostrowski con $\rho(G'(x_*))=0$).
> 3. La convergencia es **superlineal**: $\displaystyle\lim_k\frac{\|x_{k+1}-x_*\|}{\|x_k-x_*\|}=0$.

Demostración: [[Teorema 6.23]]. El paso clave, $G'(x_*)=0$, es el [[6.8 - Derivada de G en Newton es nula|Problema 6.8]] ([[Derivada de la iteracion de Newton en la raiz]]).

> [!pendiente] Qué pasa si $F'(x_*)$ es singular
> El teorema no aplica y la convergencia típicamente se degrada a **lineal**: $F(x)=x^2$ da $x_{k+1}=x_k/2$ ([[6.14 - Newton para x al cuadrado|Problema 6.14]]). Raíces múltiples ⟹ Newton lineal.

## Variante: Newton simplificado

$x_{k+1}=x_k-F'(x_0)^{-1}F(x_k)$: se factoriza $F'(x_0)$ **una sola vez** y se reusa (cada paso cuesta $O(n^2)$ en vez de $O(n^3)$). Converge localmente, pero solo **linealmente** ([[6.16 - Newton simplificado|Problema 6.16]]): en el sistema del 6.10 la tasa medida es $\approx0.02$ por paso contra la duplicación de dígitos de Newton.

## Ver también
- [[Metodo de Newton - convergencia cuadratica]] (con $F'$ Lipschitz: orden 2)
- [[Convergencia global de Newton - Newton-Baluev]]
- [[Iteracion de punto fijo y teorema de Ostrowski]]
- [[Tablas Comparativas - Punto fijo vs Newton]]
- [[Capitulo 6 - MOC]]
