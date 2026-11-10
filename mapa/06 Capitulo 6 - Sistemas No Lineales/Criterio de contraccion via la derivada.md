---
tags: [capitulo6, teorema, definicion]
---

# Criterio de contracción vía la derivada (Teoremas 6.16–6.18, Corolario 6.19)

> [!idea]
> En varias variables **no hay** Teorema del Valor Medio con igualdad ([[6.6 - No hay teorema del valor medio vectorial|Problema 6.6]]), pero sí una versión integral (6.16) y una **desigualdad** (6.18): $\|G(y)-G(x)\|\leq\max_{t}\|G'(x+t(y-x))\|\,\|y-x\|$. Consecuencia práctica (6.19): en un **convexo**, si $\|G'\|\leq\alpha<1$, $G$ es una contracción. Es la herramienta para verificar las hipótesis del [[Contracciones y teorema de punto fijo|Teorema 6.12]].

## Definiciones 6.14 y 6.15

- **Integral de $H:[a,b]\to\mathbb{R}^n$:** componente a componente, $\int_a^bH=\big(\int_a^bh_1,\dots,\int_a^bh_n\big)^T$.
- **Convexo:** $x,y\in D,\ t\in[0,1]\implies tx+(1-t)y\in D$ (contiene los segmentos entre sus puntos). Bolas, rectángulos, $\mathbb{R}^n$: convexos.

## Teorema 6.16 (Teorema Fundamental del Cálculo vectorial)

$D$ convexo, $G\in\mathcal C^1(D)$:
$$
G(y)-G(x)=\int_0^1G'(x+t(y-x))(y-x)\,dt,\qquad\forall x,y\in D.
$$
Se prueba aplicando el TFC de una variable a $h_i(t)=g_i(x+t(y-x))$. Ver [[Teorema 6.16]].

## Teorema 6.17

Si $H:[a,b]\to\mathbb{R}^n$ es continua, para **cualquier** norma: $\left\|\int_a^bH(t)\,dt\right\|\leq\int_a^b\|H(t)\|\,dt$. Es el [[6.7 - Demostrar Teorema 6.17|Problema 6.7]] (por sumas de Riemann). Ver [[Teorema 6.17]].

## Teorema 6.18 (desigualdad del valor medio)

$D$ convexo, $G\in\mathcal C^1(D)$:
$$
\|G(y)-G(x)\|\leq\max_{t\in[0,1]}\|G'(x+t(y-x))\|\,\|y-x\|,\qquad\forall x,y\in D,
$$
con $\|\cdot\|$ una norma de $\mathbb{R}^n$ y su norma **inducida** en matrices. Ver [[Teorema 6.18]].

## Corolario 6.19

$D$ convexo, $G\in\mathcal C^1(D)$ y $\|G'(x)\|\leq\alpha$ para todo $x\in D$ $\implies$ $\|G(x)-G(y)\|\leq\alpha\|x-y\|$. Si $\alpha<1$, **$G$ es una contracción en $D$**. Ver [[Corolario 6.19]].

## Receta para los ejercicios

Para aplicar el [[Contracciones y teorema de punto fijo|Teorema 6.12]] a un sistema $x=G(x)$:

1. **Elegir $D$** cerrado y convexo (casi siempre un rectángulo $[a_1,b_1]\times\cdots\times[a_n,b_n]$), mirando qué valores **puede** tomar cada componente de $G$ (signos, máximos de las funciones involucradas).
2. **Verificar $G(D)\subset D$** acotando cada $g_i$ sobre $D$.
3. **Calcular $G'(x)$** y acotar $\|G'(x)\|_\infty=\max_i\sum_j|\partial g_i/\partial x_j|$ (máxima suma por filas, la más cómoda de acotar a mano) sobre $D$. Típicamente hay que maximizar funciones de una variable del tipo $\frac{s}{(1+s^2)^2}$, $xe^{-x^2}$, etc.
4. **Concluir** con el Corolario 6.19 que $\gamma=\sup_D\|G'\|_\infty<1$, y aplicar 6.12.

Ejemplos resueltos: [[6.4 - Region para el teorema de contraccion|Problema 6.4]] ($\gamma=\frac{9}{8\sqrt3}\approx0.65$), [[6.5 - Contraccion para h chico|Problema 6.5]] ($\gamma\leq2.28h$).

> [!pendiente] El dominio tiene que ser convexo
> Sin convexidad el segmento entre $x$ e $y$ puede salir de $D$ y la cota de $G'$ en $D$ no dice nada sobre ese segmento. Por eso los dominios de los ejercicios son siempre rectángulos o bolas.

## Ver también
- [[Contracciones y teorema de punto fijo]]
- [[Metodo de Newton - convergencia cuadratica]] (usa 6.16 y 6.17 para la cota (6.6))
- [[Normas de vectores y matrices]] (fórmulas de $\|\cdot\|_1$, $\|\cdot\|_\infty$)
- [[Capitulo 6 - MOC]]
