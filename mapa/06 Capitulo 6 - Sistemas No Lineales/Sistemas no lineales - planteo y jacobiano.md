---
tags: [capitulo6, definicion]
---

# Sistemas no lineales - planteo y matriz jacobiana

> [!idea]
> Todo sistema de $n$ ecuaciones con $n$ incógnitas se escribe $F(x)=0$ con $F:D\subset\mathbb{R}^n\to\mathbb{R}^n$. A diferencia del caso lineal (0, 1 o infinitas soluciones), acá **cualquier cantidad** de soluciones es posible. La herramienta central de todo el capítulo es la **matriz jacobiana** $F'(x)$: la "mejor aproximación lineal" de $F$ cerca de un punto.

## El problema

$$
\begin{cases}
f_1(x_1,\dots,x_n)=0\\
\quad\vdots\\
f_n(x_1,\dots,x_n)=0
\end{cases}
\qquad\Longleftrightarrow\qquad F(x)=0,\quad F(x)=\begin{pmatrix}f_1(x)\\\vdots\\f_n(x)\end{pmatrix}. \tag{6.1–6.2}
$$

- **Ejemplo 6.1** ($n=1$): $x-\tan x+1=0$.
- **Ejemplo 6.2** ($n=2$): intersección de la circunferencia $x_1^2+x_2^2=1$ con la elipse $\frac{x_1^2}{9}+4x_2^2=1$.

Geométricamente, en el caso lineal las soluciones son la intersección de $n$ hiperplanos (única, infinitas o ninguna); en el no lineal, de $n$ **hipersuperficies**: puede no haber ninguna, haber finitas o infinitas.

**Dónde aparecen:**
- Minimizar $f:\mathbb{R}^n\to\mathbb{R}$ diferenciable ⟹ el mínimo $x_*$ cumple $\nabla f(x_*)=0$, un sistema no lineal con $F=\nabla f$ (ver [[6.13 - Minimo con Newton|Problema 6.13]]).
- Métodos implícitos para problemas de valores iniciales (cada paso resuelve un sistema no lineal).

## Definición 6.3 (función diferenciable, matriz jacobiana)

$D\subset\mathbb{R}^n$ abierto, $F:D\to\mathbb{R}^n$. $F$ es **diferenciable en $x_0$** si existe $J\in\mathbb{R}^{n\times n}$ con
$$
\lim_{x\to x_0}\frac{\|F(x)-F(x_0)-J(x-x_0)\|}{\|x-x_0\|}=0.
$$
$J$ es única, se llama **matriz jacobiana** y se denota $F'(x_0)$:
$$
F'(x_0)=\begin{pmatrix}
\frac{\partial f_1}{\partial x_1}(x_0) & \cdots & \frac{\partial f_1}{\partial x_n}(x_0)\\
\vdots & \ddots & \vdots\\
\frac{\partial f_n}{\partial x_1}(x_0) & \cdots & \frac{\partial f_n}{\partial x_n}(x_0)
\end{pmatrix}
=\begin{pmatrix}\nabla f_1(x_0)^T\\\vdots\\\nabla f_n(x_0)^T\end{pmatrix}.
$$
Diferenciable en $x_0$ ⟹ existen todas las derivadas parciales en $x_0$ (la recíproca es falsa).

> [!idea] Lectura geométrica (nota propia del cuaderno)
> $x\mapsto F(x_0)+J(x-x_0)$ es el "plano tangente" a $F$ en $x_0$: la definición dice que el error de aproximar $F$ por ese plano es $o(\|x-x_0\|)$. El [[Metodo de Newton - definicion y buena definicion|método de Newton]] consiste exactamente en reemplazar $F$ por este plano tangente y resolver el sistema **lineal** resultante.

## Definición 6.4 (continuamente diferenciable)

$F$ es **continuamente diferenciable en $x_0$** si todas las $\frac{\partial f_i}{\partial x_j}$ son continuas en $x_0$. Es la hipótesis que se usa en casi todos los teoremas del capítulo (Ostrowski, [[Criterio de contraccion via la derivada|Corolario 6.19]], Newton).

## Cómo calcular el jacobiano en Julia

```julia
# A mano (lo más claro para sistemas chicos)
F(x, y)  = [x^2 + y^2 - 4, x^2 - y^2 - 1]
DF(x, y) = [2x 2y; 2x -2y]

# Simbólico, para no equivocarse derivando
using Symbolics
@variables x y
J = Symbolics.jacobian([x^2 + y^2 - 4, x^2 - y^2 - 1], [x, y])
DFnum = eval(build_function(J, x, y)[1])   # función numérica
```

## Ver también
- [[Iteracion de punto fijo y teorema de Ostrowski]] (siguiente tema)
- [[Metodo de Newton - definicion y buena definicion]]
- [[Capitulo 6 - MOC]]
