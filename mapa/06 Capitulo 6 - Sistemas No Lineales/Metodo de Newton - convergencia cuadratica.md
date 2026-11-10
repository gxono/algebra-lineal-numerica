---
tags: [capitulo6, teorema, observacion]
---

# Método de Newton - convergencia cuadrática

> [!idea]
> Si además de $F'(x_*)$ invertible pedimos que $F'$ sea **Lipschitz en $x_*$** ($\|F'(x)-F'(x_*)\|\leq L\|x-x_*\|$), Newton converge **cuadráticamente**: $\|x_{k+1}-x_*\|\leq c\|x_k-x_*\|^2$, con $c=3L\|F'(x_*)^{-1}\|$. En la práctica: **la cantidad de dígitos correctos se duplica en cada paso**.

## Teorema 6.24

> $D$ abierto, $F$ diferenciable en $D$, $F(x_*)=0$, $F'(x_*)$ invertible, y existen $\delta,L>0$ con
> $$\|F'(x)-F'(x_*)\|\leq L\|x-x_*\|,\qquad\forall x\in B(x_*,\delta). \tag{6.5}$$
> Entonces existe $c>0$ tal que, si $x_0$ está suficientemente cerca de $x_*$, $\|x_{k+1}-x_*\|\leq c\|x_k-x_*\|^2$.

Demostración: [[Teorema 6.24]]. Los dos ingredientes:
1. **Cota (6.6)**, de (6.5) + [[Criterio de contraccion via la derivada|Teoremas 6.16 y 6.17]]:
$$
\|F(x)-F(x_*)-F'(x_*)(x-x_*)\|\leq\tfrac{L}{2}\|x-x_*\|^2. \tag{6.6}
$$
2. Escribir $G(x)-x_*=-F'(x)^{-1}\big[F(x)-F(x_*)-F'(x_*)(x-x_*)-(F'(x)-F'(x_*))(x-x_*)\big]$ y acotar con (6.6), (6.5) y $\|F'(x)^{-1}\|\leq2\|F'(x_*)^{-1}\|$ ([[Corolario 6.22]]).

## Observación 6.25 (cuándo vale (6.5))

Cadena de condiciones suficientes, de la más débil a la más fuerte:
$$
F'\text{ Lipschitz cerca de }x_*\ \Longleftarrow\ \frac{\partial^2f_i}{\partial x_j\partial x_k}\text{ acotadas cerca de }x_*\ \Longleftarrow\ F\in\mathcal C^2\text{ cerca de }x_*.
$$
En la práctica: **$F$ dos veces continuamente diferenciable y $F'(x_*)$ invertible ⟹ Newton cuadrático**.

## Orden de convergencia y cómo medirlo

Una sucesión $x_k\to x_*$ tiene **orden $p$** si $\|x_{k+1}-x_*\|\approx C\|x_k-x_*\|^p$ (ver [[Orden de convergencia]]). Con tres errores consecutivos $e_k=\|x_k-x_*\|$:
$$
p\approx\frac{\log(e_{k+1}/e_k)}{\log(e_k/e_{k-1})}.
$$

```julia
orden(e) = [log(e[k+1]/e[k]) / log(e[k]/e[k-1]) for k in 2:length(e)-1]
```

> [!pendiente] En Float64 solo se ven 2 o 3 pasos útiles
> Con convergencia cuadrática, el error pasa de $10^{-4}$ a $10^{-8}$ a $10^{-16}$ y ya está en el piso de la aritmética: quedan poquísimos cocientes para estimar $p$, y el último suele salir basura. Para **ver** el orden con claridad, repetir la cuenta en `BigFloat` (`setprecision(BigFloat, 500)`), como en el [[6.10 - Orden experimental de Newton|Problema 6.10]].

## Cuando falla la hipótesis (6.5): ejemplos del capítulo

| $F$ | Qué falla | Orden observado |
|---|---|---|
| $x^2+y^2=4,\ x^2-y^2=1$ ([[6.10 - Orden experimental de Newton\|6.10]]) | nada | **2** |
| $x^2$ ([[6.14 - Newton para x al cuadrado\|6.14]]) | $F'(x_*)=0$ singular | **1** (tasa $1/2$) |
| $x+x^{1+\alpha}$, $0<\alpha<1$ ([[6.15 - Newton para x mas x a la 1 mas alfa\|6.15]]) | $F'$ no es Lipschitz en $0$ (solo Hölder-$\alpha$) | **$1+\alpha$** (superlineal, no cuadrático) |
| $Ax+h^2x^3$ ([[6.17 - Newton en problema discretizado\|6.17]]) | nada, y además $F''(x_*)=0$ | **3** |

La última fila muestra que el teorema da una **cota inferior** del orden: si además la parte cuadrática del error se anula, se puede tener más.

## Ver también
- [[Metodo de Newton - definicion y buena definicion]]
- [[Orden de convergencia]]
- [[Convergencia global de Newton - Newton-Baluev]]
- [[Capitulo 6 - MOC]]
