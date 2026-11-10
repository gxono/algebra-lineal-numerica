---
tags: [capitulo6, teorema]
---

# Teorema 6.16 (Teorema Fundamental del Cálculo vectorial)

Sea $D\subset\mathbb{R}^n$ convexo y $G:D\to\mathbb{R}^n$ continuamente diferenciable en $D$. Entonces
$$
G(y)-G(x)=\int_0^1G'(x+t(y-x))(y-x)\,dt,\qquad\forall x,y\in D.
$$

> [!demostracion]-
> Sea $G=(g_1,\dots,g_n)^T$ y $x,y\in D$ fijos. Por convexidad, $x+t(y-x)\in D$ para $t\in[0,1]$. Para cada $i$ definimos $h_i(t):=g_i(x+t(y-x))$. Por la regla de la cadena,
> $$h_i'(t)=\nabla g_i(x+t(y-x))^T(y-x),$$
> que es continua en $t$ (porque $G\in\mathcal C^1$). Por el Teorema Fundamental del Cálculo de una variable,
> $$g_i(y)-g_i(x)=h_i(1)-h_i(0)=\int_0^1h_i'(t)\,dt=\int_0^1\nabla g_i(x+t(y-x))^T(y-x)\,dt.$$
> La fila $i$ de $G'(z)$ es $\nabla g_i(z)^T$, así que juntando las $n$ componentes (la integral vectorial es componente a componente, Def. 6.14) se obtiene la identidad.

> [!idea]
> Esta es la versión **correcta** del valor medio en varias variables: un *promedio* de jacobianos a lo largo del segmento, no el jacobiano en *un* punto intermedio. Ese punto intermedio puede no existir ([[6.6 - No hay teorema del valor medio vectorial|Problema 6.6]]).

## Ver también
- [[Criterio de contraccion via la derivada]]
- [[Teorema 6.18]], [[Teorema 6.24]] (lo usan)
- [[Demostraciones - Capitulo 6|Demostraciones]]
