---
tags: [capitulo6, definicion, teorema, algoritmo]
---

# Iteración de punto fijo y teorema de Ostrowski

> [!idea]
> Reescribir $F(x)=0$ como $x=G(x)$ e iterar $x_{k+1}=G(x_k)$. Es la generalización no lineal de los [[Metodos iterativos estacionarios - marco general|métodos iterativos estacionarios]] ($G(x)=Mx+c$), y el rol de $M$ lo juega ahora el jacobiano $G'(x_*)$: si $\rho(G'(x_*))<1$, la iteración converge **localmente** (Ostrowski). Diferencia clave con el caso lineal: acá es solo **suficiente** y solo **local**.

## La iteración (6.3)

Dado $F(x)=0$, se lo transforma en $x=G(x)$ con $G:D\subset\mathbb{R}^n\to\mathbb{R}^n$ (la **función de iteración**), y se itera
$$
x_{k+1}=G(x_k),\qquad k=0,1,\dots \tag{6.3}
$$

> [!algoritmo] Punto fijo
> ```julia
> function punto_fijo(G, x0; tol=1e-10, maxit=1000)
>     x = x0
>     for k in 1:maxit
>         xn = G(x)
>         norm(xn - x) < tol && return xn, k
>         x = xn
>     end
>     @warn "no convergió en $maxit iteraciones"
>     return x, maxit
> end
> ```

La elección de $G$ **no es única** y decide todo: para el mismo $F$, una $G$ puede converger y otra divergir.

## Observación 6.5 (el caso lineal)

Si $F(x)=Ax-b$ y se reescribe como $x=Mx+c$, (6.3) es exactamente un método iterativo estacionario ([[Metodo de Jacobi|Jacobi]], [[Metodo de Gauss-Seidel|Gauss-Seidel]], [[Metodo de Richardson|Richardson]]...). El Capítulo 6 es la versión no lineal del Capítulo 4.

## Definiciones 6.6 y 6.7

- **Punto fijo:** $x_*=G(x_*)$.
- **Punto de atracción:** $x_*$ punto fijo tal que existe $\delta>0$ con: para todo $x_0\in B(x_*,\delta)$, los iterados quedan en $B(x_*,\delta)$ y $x_k\to x_*$. Equivalentemente: *la iteración converge localmente a $x_*$*.

## Teorema 6.8 (Ostrowski)

> Si $G:\mathbb{R}^n\to\mathbb{R}^n$ es diferenciable en un punto fijo $x_*$ y $\rho(G'(x_*))<1$, entonces $x_*$ es un punto de atracción.

Idea de la prueba: por el [[Teorema 4.20]] hay una norma inducida con $\|G'(x_*)\|<\sigma+\varepsilon$ ($\sigma=\rho(G'(x_*))$, $\sigma+2\varepsilon<1$); cerca de $x_*$, $G$ se parece a su linealización, así que $\|x_{k+1}-x_*\|\leq(\sigma+2\varepsilon)\|x_k-x_*\|$. Demostración completa: [[Teorema 6.8]].

De la prueba sale también la **tasa**: $\|x_k-x_*\|\leq(\sigma+2\varepsilon)^k\|x_0-x_*\|$, o sea convergencia lineal con factor asintótico $\approx\rho(G'(x_*))$.

## Observación 6.9 (comparación con el caso lineal)

| | Lineal ($G=Mx+c$, Cap. 4) | No lineal (Ostrowski) |
|---|---|---|
| Condición | $\rho(M)<1$ | $\rho(G'(x_*))<1$ |
| Tipo | **necesaria y suficiente** | solo **suficiente** |
| Alcance | **global** (todo $x_0$) | solo **local** ($x_0$ cerca de $x_*$) |
| Hace falta conocer $x_*$ | no | **sí** (hay que evaluar $G'$ en $x_*$) |

> [!pendiente] $\rho(G'(x_*))=1$ no decide nada
> Con $\rho=1$ el teorema no dice nada: puede converger (lentísimo) o no. El [[6.1 - Rho igual a 1 no decide|Problema 6.1]] lo muestra con $G(x)=x-x^3$ (converge, como $1/\sqrt{2k}$) y $G(x)=x+x^3$ (no converge), ambas con $G'(0)=1$.

## Una variante útil: contracción en una bola

Si además $G$ es $\mathcal C^1$ cerca de $x_*$, $\rho(G'(x_*))<1$ implica que $G$ es una **contracción** en alguna bola cerrada centrada en $x_*$ ([[6.2 - Contraccion local si rho menor que 1|Problema 6.2]]), y ahí vale todo lo del [[Contracciones y teorema de punto fijo|Teorema 6.12]] (unicidad, cota de error a priori).

## Ver también
- [[Contracciones y teorema de punto fijo]] (resultado global, sin conocer $x_*$ de antemano)
- [[Teorema de convergencia de metodos iterativos estacionarios]] (Capítulo 4, el caso lineal)
- [[Metodo de Newton - definicion y buena definicion]] (Newton es un punto fijo con $G'(x_*)=0$)
- [[Capitulo 6 - MOC]]
