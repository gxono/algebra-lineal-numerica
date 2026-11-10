---
tags: [capitulo6, teorema, definicion, observacion]
---

# Convergencia global de Newton - Teorema de Newton-Baluev

> [!idea]
> Los teoremas 6.23/6.24 son **locales**: hay que arrancar cerca de $x_*$. Para garantizar convergencia **desde cualquier $x_0$** hace falta estructura: si $F$ es **convexa** (cada componente), $F'(x)$ invertible con **$F'(x)^{-1}\geq0$ entrada a entrada**, y $F$ tiene un cero, entonces ese cero es único y Newton converge desde cualquier $x_0$, y además **en forma monótona decreciente** componente a componente desde $x_1$.

## Observación 6.26

Si $G$ es contracción en **todo** $\mathbb{R}^n$, el [[Contracciones y teorema de punto fijo|Teorema 6.12]] ya da convergencia global. Pero esa hipótesis casi nunca se cumple en problemas concretos.

## Definición 6.27 (función convexa)

$f:\mathbb{R}^n\to\mathbb{R}$ es convexa si $f(tx+(1-t)y)\leq tf(x)+(1-t)f(y)$ para todo $t\in[0,1]$. $F=(f_1,\dots,f_n)$ es convexa si cada $f_i$ lo es.

## Proposición 6.28 (caracterización con la derivada)

$f$ diferenciable es convexa $\iff$ $f(y)-f(x)\geq\nabla f(x)^T(y-x)$ para todo $x,y$. Para $n=1$: **el gráfico queda por encima de todas sus rectas tangentes**. Demostración: [[Proposicion 6.28]].

## Teorema 6.29 (Newton-Baluev)

> $F:\mathbb{R}^n\to\mathbb{R}^n$ convexa, $\mathcal C^1$, con $F'(x)$ invertible y $F'(x)^{-1}\geq0$ (entradas no negativas) para todo $x$. Si existe $x_*$ con $F(x_*)=0$, entonces es **único** y Newton converge a $x_*$ desde **cualquier** $x_0\in\mathbb{R}^n$. Además, componente a componente,
> $$x_*^{(i)}\leq\cdots\leq x_{k+1}^{(i)}\leq x_k^{(i)}\leq\cdots\leq x_1^{(i)}.$$

Demostración: [[Teorema 6.29]] (incluye la unicidad, que el apunte deja como ejercicio). Estructura:
1. Convexidad + definición de Newton ⟹ $F(x_k)\geq0$ para $k\geq1$ (Newton "se pasa" siempre del mismo lado).
2. $F'(x_k)^{-1}\geq0$ y $F(x_k)\geq0$ ⟹ $x_{k+1}\leq x_k$ (decreciente).
3. Convexidad en $x_*$ ⟹ $x_k\geq x_*$ (acotada inferiormente).
4. Monótona y acotada ⟹ converge; el límite es un cero por continuidad.

## Observación 6.30 (caso $n=1$)

Las hipótesis se reducen a: $F\in\mathcal C^1$, **convexa** y **$F'>0$**. Entonces Newton converge decrecientemente al único cero desde cualquier $x_0$ (desde $x_1$, que queda a la derecha de $x_*$).

## Las hipótesis no son decorativas

> [!pendiente] Sin convexidad, Newton puede no converger globalmente aunque todo lo demás valga
> [[6.18 - Newton-Baluev sin convexidad|Problema 6.18]]: $F(x)=\arctan x$ y $F(x)=2x+\sin x$ son $\mathcal C^1$, con $F'>0$ y un único cero en $0$, pero **no son convexas**, y Newton falla:
> - $\arctan$: diverge oscilando si $|x_0|>x_c\approx1.3917$ (raíz de $(1+x^2)\arctan x=2x$).
> - $2x+\sin x$: desde $x_0=\pi$, Newton salta $\pi\to-\pi\to\pi\to\cdots$ para siempre, y ese 2-ciclo es **superatractor**: un conjunto abierto de puntos iniciales cae en él.

Ejemplo típico donde **sí** aplica: sistemas $Ax+\Phi(x)=b$ con $A$ una M-matriz (como la tridiagonal $(-1,2,-1)$, cuya inversa tiene entradas positivas) y $\Phi$ convexa y creciente componente a componente (ej. $\phi(t)=e^t$), que salen de discretizar EDOs no lineales.

## Ver también
- [[Metodo de Newton - definicion y buena definicion]]
- [[Metodo de Newton - convergencia cuadratica]]
- [[6.18 - Newton-Baluev sin convexidad]]
- [[Capitulo 6 - MOC]]
