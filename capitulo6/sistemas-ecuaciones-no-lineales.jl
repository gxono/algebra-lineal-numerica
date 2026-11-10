### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000001
using LinearAlgebra, PlutoUI, CairoMakie

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000002
TableOfContents()

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000003
md"""
# Capítulo 6: Sistemas de ecuaciones no lineales
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000004
md"""
En este capítulo estudiaremos métodos para resolver sistemas de ecuaciones no lineales. En general, para ``n\in\mathbb{N}``, un sistema de ``n`` ecuaciones con ``n`` incógnitas ``(x_1,x_2,\dots,x_n)`` viene dado por
```math
\begin{cases}
f_1(x_1,x_2,\dots,x_n)=0\\
f_2(x_1,x_2,\dots,x_n)=0\\
\qquad\vdots\\
f_n(x_1,x_2,\dots,x_n)=0,
\end{cases} \tag{6.1}
```
donde ``f_i:D\subset\mathbb{R}^n\to\mathbb{R}``, para ``i=1,2,\dots,n``. El sistema (6.1) siempre puede escribirse como
```math
F(x) = 0, \tag{6.2}
```
donde ``F:D\subset\mathbb{R}^n\to\mathbb{R}^n`` está dada por ``F(x)=\begin{pmatrix}f_1(x)\\f_2(x)\\\vdots\\f_n(x)\end{pmatrix}``, para ``x\in D``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000005
md"""
**Ejemplo 6.1.** Para ``n=1``, un sistema de una ecuación con una incógnita es simplemente *una ecuación no lineal*, por ejemplo, hallar ``x\in\mathbb{R}`` tal que
```math
x-\tan x+1 = 0.
```

**Ejemplo 6.2.** Para ``n=2``, se tiene un sistema de dos ecuaciones con dos incógnitas, cuya solución en general representa la intersección de dos curvas en el plano, por ejemplo, hallar ``(x_1,x_2)\in\mathbb{R}^2`` tales que
```math
\begin{cases}
x_1^2+x_2^2=1\\
\frac{x_1^2}{9}+4x_2^2=1.
\end{cases}
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000006
md"""
Cuando un sistema de ``n`` ecuaciones con ``n`` incógnitas es lineal, la solución o las soluciones representan la intersección de ``n`` hiperplanos en ``\mathbb{R}^n``, y en este caso sabemos que existen sólo tres posibilidades, a saber: tiene una única solución, tiene infinitas soluciones, o no tiene ninguna. En el caso de un sistema *no lineal* de ``n`` ecuaciones con ``n`` incógnitas la solución o las soluciones (si existen) corresponden a la intersección de ``n`` hipersuperficies en ``\mathbb{R}^n`` y pueden darse todas las posibilidades, es decir, puede no tener solución, o puede tener una cantidad finita o infinita de soluciones.

La necesidad de resolver sistemas de ecuaciones no lineales surge, por ejemplo, cuando se quiere hallar el mínimo de una función de varias variables ``f:\mathbb{R}^n\to\mathbb{R}``. Sabemos que si ``f`` es diferenciable y alcanza su mínimo en un valor ``x_*\in\mathbb{R}^n`` entonces
```math
\nabla f(x_*) = 0,
```
es decir, ``x_*`` debe ser solución del sistema no lineal (6.2) con ``F(x):=\nabla f(x)``.

Otra situación en donde surge la necesidad de resolver un sistema de ecuaciones no lineales, es al aplicar métodos implícitos para resolver problemas de valores iniciales (PVI).

En este capítulo estudiaremos métodos numéricos para resolver sistemas de ecuaciones no lineales. Antes de presentar estos métodos, recordamos algunos conceptos de cálculo en varias variables. En lo que sigue, ``\|\cdot\|`` denota una norma en ``\mathbb{R}^n``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000007
md"""
**Definición 6.3** (Función diferenciable). Sea ``D\subset\mathbb{R}^n`` abierto y sea ``F:D\to\mathbb{R}^n``. Decimos que ``F`` *es diferenciable en* ``x_0\in D`` si existe ``J\in\mathbb{R}^{n\times n}`` tal que
```math
\lim_{x\to x_0}\frac{\|F(x)-F(x_0)-J(x-x_0)\|}{\|x-x_0\|} = 0.
```
En este caso, la matriz ``J`` se llama *matriz jacobiana* de ``F`` en ``x_0`` y se denota por ``F'(x_0)``.

Se puede demostrar que si ``F`` es diferenciable en ``x_0``, la matriz jacobiana es única y viene dada por
```math
F'(x_0) = \begin{pmatrix}
\frac{\partial f_1}{\partial x_1}(x_0) & \frac{\partial f_1}{\partial x_2}(x_0) & \cdots & \frac{\partial f_1}{\partial x_n}(x_0)\\
\frac{\partial f_2}{\partial x_1}(x_0) & \frac{\partial f_2}{\partial x_2}(x_0) & \cdots & \frac{\partial f_2}{\partial x_n}(x_0)\\
\vdots & \vdots & \ddots & \vdots\\
\frac{\partial f_n}{\partial x_1}(x_0) & \frac{\partial f_n}{\partial x_2}(x_0) & \cdots & \frac{\partial f_n}{\partial x_n}(x_0)
\end{pmatrix}
= \begin{pmatrix}\nabla f_1(x_0)^T\\\nabla f_2(x_0)^T\\\vdots\\\nabla f_n(x_0)^T\end{pmatrix},
\quad\text{si } F(x)=\begin{pmatrix}f_1(x)\\f_2(x)\\\vdots\\f_n(x)\end{pmatrix};
```
es decir, diferenciabilidad en ``x_0`` implica existencia de todas las derivadas parciales en ``x_0``.

**Definición 6.4** (Función continuamente diferenciable). Sea ``D\subset\mathbb{R}^n`` abierto y sea ``F:D\to\mathbb{R}^n`` con ``F(x)=\begin{pmatrix}f_1(x)\\f_2(x)\\\vdots\\f_n(x)\end{pmatrix}``. Decimos que ``F`` *es continuamente diferenciable en* ``x_0\in D`` si ``\frac{\partial f_i}{\partial x_j}`` es continua en ``x_0``, para ``i,j=1,2,\dots,n``.
"""

# ╔═╡ b8cd7866-eaba-4c17-901d-e7ce33fa3ade
@info "F(x₀) + J ⋅ (x - x₀) el conjunto tangente en la definicion"

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000008
md"""
## 6.1 Métodos de punto fijo

Dado el sistema ``F(x)=0``, queremos transformarlo en una ecuación (vectorial) de la forma ``x=G(x)``, donde ``G:D\subset\mathbb{R}^n\to\mathbb{R}^n`` es, en general, una función no lineal. La *iteración de punto fijo* asociada a la función ``G`` se define por:
```math
x_{k+1} = G(x_k), \qquad k=0,1,\dots \tag{6.3}
```
donde ``x_0\in D`` es un valor inicial dado. En este caso, la función ``G`` se llama *función de iteración*.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000009
md"""
**Observación 6.5.** Sea ``F(x)=Ax-b``, donde ``A\in\mathbb{R}^{n\times n}`` y ``b\in\mathbb{R}^n``. Si reescribimos el sistema ``F(x)=0`` como ``x=G(x)`` donde ``G(x)=Mx+c``, para ``M\in\mathbb{R}^{n\times n}`` y ``c\in\mathbb{R}^n``, entonces la iteración de punto fijo (6.3) es un método iterativo estacionario para el sistema lineal ``F(x)=0`` (``Ax=b``). Así, podemos pensar la iteración de punto fijo (6.3) como una generalización al caso de sistemas no lineales de los métodos iterativos estacionarios definidos para sistemas lineales.

**Definición 6.6** (Punto fijo). Decimos que ``x_*\in D`` es *un punto fijo* de ``G:D\subset\mathbb{R}^n\to\mathbb{R}^n`` si
```math
x_* = G(x_*).
```

**Definición 6.7** (Punto de atracción). Sea ``G:D\subset\mathbb{R}^n\to\mathbb{R}^n`` y sea ``x_*\in D`` un punto fijo de ``G``. Decimos que ``x_*`` es *un punto de atracción de* ``G``, o equivalentemente, que *la iteración (6.3) converge localmente a* ``x_*``, si existe ``\delta>0`` tal que para todo ``x_0\in B(x_*,\delta)`` se cumple que ``x_{k+1}:=G(x_k)\in B(x_*,\delta)``, para todo ``k=0,1,\dots`` y
```math
\lim_{k\to\infty}x_k = x_*.
```

El siguiente resultado nos da una condición suficiente para determinar cuándo un punto fijo de ``G`` es un punto de atracción.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000000a
md"""
**Teorema 6.8** (Teorema de Ostrowski). Sea ``G:\mathbb{R}^n\to\mathbb{R}^n`` diferenciable en un punto ``x_*\in\mathbb{R}^n`` tal que ``x_*=G(x_*)``. Si ``\rho(G'(x_*))<1`` entonces ``x_*`` es un punto de atracción.

**Observación 6.9.** Notemos que este resultado es una generalización del teorema de convergencia de métodos iterativos estacionarios para sistemas lineales (ver Teorema 4.25), ya que si ``G(x)=Mx+c``, se tiene que ``G'(x)=M`` para todo ``x\in\mathbb{R}^n``. Notemos en particular que en el caso de que el sistema es no lineal este teorema es sólo una condición suficiente para la convergencia y por otro lado, que ésta se garantiza sólo localmente, es decir, partiendo de un valor inicial ``x_0`` suficientemente cerca de ``x_*``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000000b
md"""
*Demostración del Teorema 6.8.* Sea ``\sigma:=\rho(G'(x_*))<1`` y elijamos ``\varepsilon>0`` (fijo) tal que ``\sigma+2\varepsilon<1``. Por el Teorema 4.20, sabemos que existe una norma inducida ``\|\cdot\|`` tal que
```math
\|G'(x_*)\| < \sigma+\varepsilon.
```
Por otro lado, puesto que ``G`` es diferenciable en ``x_*``, existe ``\delta>0`` tal que si ``x\in B(x_*,\delta)`` entonces
```math
\|G(x)-G(x_*)-G'(x_*)(x-x_*)\| < \varepsilon\|x-x_*\|.
```
Sea ``x_0\in B(x_*,\delta)`` y ``x_{k+1}:=G(x_k)``, para ``k=0,1,2,\dots``. Notemos que si ``x_k\in B(x_*,\delta)``, para algún ``k\in\mathbb{N}_0``,
```math
\begin{align*}
\|x_{k+1}-x_*\| = \|G(x_k)-G(x_*)\| &\leq \|G(x_k)-G(x_*)-G'(x_*)(x_k-x_*)\|+\|G'(x_*)(x_k-x_*)\|\\
&< \varepsilon\|x_k-x_*\|+(\sigma+\varepsilon)\|x_k-x_*\| = (\sigma+2\varepsilon)\|x_k-x_*\| < \|x_k-x_*\|,
\end{align*}
```
y en consecuencia, ``x_{k+1}\in B(x_*,\delta)``. Así, tenemos que ``x_k\in B(x_*,\delta)``, para todo ``k\in\mathbb{N}_0``. Por otro lado, de la última estimación se sigue también que
```math
\|x_{k+1}-x_*\| \leq (\sigma+2\varepsilon)^{k+1}\|x_0-x_*\|, \qquad k=0,1,\dots
```
Puesto que ``\sigma+2\varepsilon<1``, tenemos que ``x_k\to x_*`` cuando ``k\to\infty``, lo que concluye la demostración.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000000c
md"""
En el resultado anterior se necesita conocer de antemano la existencia de un punto fijo de la función de iteración, y más aún, se pide que ésta sea diferenciable en dicho punto. A continuación establecemos un resultado que garantiza la existencia de un único punto fijo y la convergencia de la iteración correspondiente sin pedir diferenciabilidad. Recordamos antes dos conceptos que serán importantes.

**Definición 6.10** (Función Lipschitz). Decimos que ``G:D\subset\mathbb{R}^n\to\mathbb{R}^n`` es *Lipschitz continua* sobre ``D`` con constante ``\gamma>0`` si
```math
\|G(x)-G(y)\| \leq \gamma\|x-y\|, \qquad \forall x,y\in D.
```

**Definición 6.11** (Contracción). Decimos que ``G:D\subset\mathbb{R}^n\to\mathbb{R}^n`` es *una contracción* en ``D`` con la norma ``\|\cdot\|`` si ``G`` es Lipschitz continua con constante ``0<\gamma<1``, es decir, si existe ``\gamma\in(0,1)`` tal que
```math
\|G(x)-G(y)\| \leq \gamma\|x-y\|, \qquad \forall x,y\in D.
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000000d
md"""
**Teorema 6.12** (Contraction Mapping Theorem). Sea ``D\subset\mathbb{R}^n`` cerrado. Sea ``G:D\to D`` una contracción en ``D`` (con alguna norma ``\|\cdot\|``) con constante ``0<\gamma<1``. Entonces, existe un único punto fijo ``x_*`` de ``G`` en ``D`` y ``x_*`` es un punto de atracción de ``G``. Más aún, para todo ``x_0\in D``, la iteración
```math
x_{k+1} = G(x_k), \qquad k=0,1,\dots
```
está bien definida y se cumplen:

- Convergencia lineal:
```math
\|x_{k+1}-x_*\| \leq \gamma\|x_k-x_*\|, \qquad k=0,1,\dots
```

- Estimación del error:
```math
\|x_k-x_*\| \leq \frac{\gamma^k}{1-\gamma}\|x_1-x_0\|, \qquad k=0,1,\dots
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000000e
md"""
*Demostración.* Sea ``x_0\in D``. Puesto que ``G:D\to D``, la iteración
```math
x_{k+1} := G(x_k), \qquad k=0,1,2,\dots
```
está bien definida. Demostraremos ahora que ``\{x_k\}_{k\in\mathbb{N}_0}`` es una sucesión de Cauchy.

Puesto que ``G`` es una contracción, tenemos que
```math
\|x_{i+1}-x_i\| = \|G(x_i)-G(x_{i-1})\| \leq \gamma\|x_i-x_{i-1}\| \leq \dots \leq \gamma^i\|x_1-x_0\|,
```
y luego, para ``k\in\mathbb{N}_0`` y ``p\in\mathbb{N}``,
```math
\begin{align*}
\|x_{k+p}-x_k\| &\leq \sum_{i=k}^{k+p-1}\|x_{i+1}-x_i\| \leq \left(\sum_{i=k}^{k+p-1}\gamma^i\right)\|x_1-x_0\|\\
&= \frac{\gamma^k-\gamma^{k+p}}{1-\gamma}\|x_1-x_0\| \leq \frac{\gamma^k}{1-\gamma}\|x_1-x_0\|.
\end{align*}
```
Teniendo en cuenta que ``\gamma<1``, concluimos que ``\{x_k\}_{k\in\mathbb{N}}`` es una sucesión de Cauchy en ``\mathbb{R}^n``. En consecuencia, existe ``x_*\in\mathbb{R}^n`` tal que
```math
\lim_{k\to\infty}x_k = x_*.
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000000f
md"""
Puesto que ``D`` es cerrado y que ``\{x_k\}_{k\in\mathbb{N}}\subset D``, se tiene que ``x_*\in D``. Por otro lado, puesto que ``G`` es continua se tiene que
```math
\lim_{k\to\infty}G(x_k) = G(x_*),
```
y como ``G(x_k)=x_{k+1}``, tenemos que ``x_*=G(x_*)``, es decir, que ``x_*`` es un punto fijo de ``G``.

Es fácil verificar que ``x_*`` es el único punto fijo de ``G``. En efecto, si ``y_*\in D`` es un punto fijo de ``G`` distinto de ``x_*``, entonces
```math
\|y_*-x_*\| = \|G(y_*)-G(x_*)\| \leq \gamma\|y_*-x_*\| < \|y_*-x_*\|,
```
lo cual no puede ocurrir.

Hasta aquí hemos demostrado que ``G`` tiene un único punto fijo ``x_*\in D``, y que para cualquier ``x_0\in D``, la iteración de punto fijo ``x_{k+1}=G(x_k)``, para ``k\in\mathbb{N}_0``, está bien definida y converge a ``x_*``. Sólo resta ver que la convergencia es lineal y la estimación de error.

Usando nuevamente que ``G`` es una contracción, tenemos que
```math
\|x_{k+1}-x_*\| = \|G(x_k)-G(x_*)\| \leq \gamma\|x_k-x_*\|, \qquad k=0,1,2,\dots,
```
y por lo tanto, ``\{x_k\}_{k\in\mathbb{N}}`` converge linealmente. Finalmente, recordando que
```math
\|x_{k+p}-x_k\| \leq \frac{\gamma^k}{1-\gamma}\|x_1-x_0\|, \qquad \forall k\in\mathbb{N}_0,\ p\in\mathbb{N},
```
y tomando el límite cuando ``p`` tiende a infinito, obtenemos la estimación de error deseada.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000010
md"""
**Observación 6.13.** Si ``D=\mathbb{R}^n`` y ``G(x):=Mx+c``, notemos que ``G`` es una contracción en alguna norma si y sólo si ``\rho(M)<1``. En efecto, si para alguna norma ``\|\cdot\|`` y algún ``0<\gamma<1`` se cumple que
```math
\|G(x)-G(y)\| \leq \gamma\|x-y\|, \qquad \forall x,y\in\mathbb{R}^n,
```
se tiene que
```math
\rho(M) \leq \|M\| = \sup_{z\in\mathbb{R}^n}\frac{\|Mz\|}{\|z\|} = \sup_{z\in\mathbb{R}^n}\frac{\|G(z)-G(0)\|}{\|z\|} \leq \gamma < 1.
```
Por otro lado, si ``\rho(M)<1``, por el Teorema 4.20 sabemos que existe una norma inducida ``\|\cdot\|`` tal que ``\|M\|<1``, y por lo tanto,
```math
\|G(x)-G(y)\| = \|Mx+c-(My+c)\| = \|M(x-y)\| \leq \|M\|\|x-y\|, \qquad \forall x,y\in\mathbb{R}^n,
```
lo que demuestra que ``G`` es una contracción en la norma ``\|\cdot\|``.

Así, podemos ver al Teorema 6.12 como una generalización al caso no lineal del Teorema 4.25, que establece una caracterización de la convergencia de métodos de punto fijo asociados a sistemas lineales.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000011
md"""
A continuación, y para concluir esta sección, nos proponemos establecer un criterio para determinar cuándo una función diferenciable ``G`` es una contracción en alguna norma (ver Corolario 6.19 debajo) y así poder utilizar el Teorema 6.12.

**Definición 6.14** (Integral de funciones vectoriales). Sea ``H:[a,b]\to\mathbb{R}^n``. Si ``H(t)=\begin{pmatrix}h_1(t)\\h_2(t)\\\vdots\\h_n(t)\end{pmatrix}``, entonces ``\displaystyle\int_a^bH(t)\,dt:=\begin{pmatrix}\int_a^bh_1(t)\,dt\\\int_a^bh_2(t)\,dt\\\vdots\\\int_a^bh_n(t)\,dt\end{pmatrix}``.

**Definición 6.15** (Conjunto convexo en ``\mathbb{R}^n``). Un conjunto ``D\subset\mathbb{R}^n`` es *convexo* si para todo par de puntos de ``D``, el segmento que los une está contenido en ``D``, esto es, si
```math
\forall x,y\in D,\quad t\in[0,1], \qquad tx+(1-t)y\in D.
```

El siguiente teorema es una generalización del Teorema Fundamental del Cálculo para funciones de varias variables.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000012
md"""
**Teorema 6.16** (Teorema Fundamental del Cálculo). Sea ``D\subset\mathbb{R}^n`` convexo. Si ``G:D\to\mathbb{R}^n`` es continuamente diferenciable en ``D`` entonces
```math
G(y)-G(x) = \int_0^1G'(x+t(y-x))(y-x)\,dt, \qquad \forall x,y\in D.
```

*Demostración.* Sea ``G=\begin{pmatrix}g_1\\g_2\\\vdots\\g_n\end{pmatrix}`` y sean ``x,y\in D`` fijos. Para ``i=1,2,\dots,n``, definimos ``h_i(t):=g_i(x+t(y-x))``, para ``t\in[0,1]``. Entonces, ``h_i'(t)=\nabla g_i(x+t(y-x))^T(y-x)`` y por el Teorema Fundamental del Cálculo tenemos que
```math
g_i(y)-g_i(x) = h_i(1)-h_i(0) = \int_0^1h_i'(t)\,dt = \int_0^1\nabla g_i(x+t(y-x))^T(y-x)\,dt,
```
lo que concluye la demostración.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000013
md"""
**Teorema 6.17.** Si ``H:[a,b]\to\mathbb{R}^n`` es continua entonces, para cualquier norma de ``\mathbb{R}^n``,
```math
\left\|\int_a^bH(t)\,dt\right\| \leq \int_a^b\|H(t)\|\,dt.
```

*Demostración.* La demostración se deja como ejercicio (ver Ejercicio 6.7).

Como consecuencia de los dos últimos teoremas tenemos el siguiente resultado.

**Teorema 6.18.** Sea ``D\subset\mathbb{R}^n`` convexo. Si ``G:D\to\mathbb{R}^n`` es continuamente diferenciable en ``D`` entonces
```math
\|G(y)-G(x)\| \leq \max_{t\in[0,1]}\|G'(x+t(y-x))\|\|y-x\|, \qquad \forall x,y\in D;
```
aquí, ``\|\cdot\|`` denota tanto una norma de ``\mathbb{R}^n`` como la norma inducida por ella en el espacio de las matrices.

*Demostración.* Por el Teorema 6.16 se tiene que
```math
G(y)-G(x) = \int_0^1G'(x+t(y-x))(y-x)\,dt, \qquad \forall x,y\in D.
```
Teniendo en cuenta ahora el Teorema 6.17 obtenemos
```math
\begin{align*}
\|G(y)-G(x)\| &\leq \int_0^1\|G'(x+t(y-x))(y-x)\|\,dt \leq \int_0^1\|G'(x+t(y-x))\|\|y-x\|\,dt\\
&\leq \max_{t\in[0,1]}\|G'(x+t(y-x))\|\|y-x\|.
\end{align*}
```

Finalmente, una aplicación inmediata del teorema anterior nos permite establecer el siguiente corolario.

**Corolario 6.19.** Sea ``D\subset\mathbb{R}^n`` convexo. Si ``G:D\to\mathbb{R}^n`` es continuamente diferenciable en ``D`` y existe ``\alpha>0`` tal que
```math
\|G'(x)\| \leq \alpha, \qquad \forall x\in D,
```
entonces
```math
\|G(x)-G(y)\| \leq \alpha\|x-y\|, \qquad \forall x,y\in D.
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000014
md"""
## 6.2 El método de Newton

Consideremos el sistema ``F(x)=0``, y supongamos que ``F`` es continuamente diferenciable en un entorno de una solución ``x_*`` del sistema. Dada una aproximación ``x_k`` de la solución, consideramos la aproximación lineal afín de ``F`` dada por
```math
P(x) := F(x_k)+F'(x_k)(x-x_k),
```
y definimos ``x_{k+1}`` tal que ``P(x_{k+1})=0``, es decir,
```math
F'(x_k)(x_{k+1}-x_k) = -F(x_k),
```
o, suponiendo que ``F'(x_k)`` es invertible,
```math
x_{k+1} = x_k-(F'(x_k))^{-1}F(x_k).
```

La iteración anterior corresponde al *método de Newton*. En la práctica no se calcula nunca explícitamente la matriz ``(F'(x_k))^{-1}``, sino que se resuelve un sistema lineal asociado a la matriz ``F'(x_k)``. Más precisamente, dado ``x_k`` se obtiene ``x_{k+1}`` de la siguiente manera:
```math
\begin{cases}
\text{Hallar } h\in\mathbb{R}^n \text{ tal que } F'(x_k)h=F(x_k).\\
x_{k+1}=x_k-h.
\end{cases} \tag{6.4}
```

Notemos que el método de Newton es una iteración de punto fijo cuya función de iteración es
```math
G(x) := x-(F'(x))^{-1}F(x).
```

La primera pregunta importante que surge es bajo qué hipótesis la iteración de Newton está bien definida. Notemos que en particular se necesita que la matriz jacobiana de ``F`` sea invertible. Los siguientes resultados nos llevaran a saber en qué circunstancias ésto ocurre.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000015
md"""
**Lema 6.20** (Perturbación de la matriz Identidad). Sea ``M\in\mathbb{R}^{n\times n}``. Si ``\|M\|<1`` para alguna norma inducida, entonces ``I-M`` es invertible y
```math
\|(I-M)^{-1}\| \leq \frac{1}{1-\|M\|}.
```

*Demostración.* Ver demostración del Teorema 4.25.

**Lema 6.21** (Perturbación de una matriz invertible). Sea ``A\in\mathbb{R}^{n\times n}`` invertible y ``E\in\mathbb{R}^{n\times n}`` tal que ``\|E\|<\frac{1}{\|A^{-1}\|}``. Entonces, ``A+E`` es invertible y
```math
\|(A+E)^{-1}\| \leq \frac{\|A^{-1}\|}{1-\|A^{-1}\|\|E\|}.
```

*Demostración.* Puesto que ``A`` es invertible, ``A+E=A(I+A^{-1}E)``. Si llamamos ``M:=-A^{-1}E``, tenemos que ``\|M\|=\|A^{-1}E\|\leq\|A^{-1}\|\|E\|<1``, y por el lema anterior tenemos que ``I-M`` es invertible y que ``\|(I-M)^{-1}\|\leq\frac{1}{1-\|M\|}``. En consecuencia, ``A+E=A(I-M)`` es invertible y
```math
\|(A+E)^{-1}\| = \|(A(I-M))^{-1}\| = \|(I-M)^{-1}A^{-1}\| \leq \frac{\|A^{-1}\|}{1-\|M\|} \leq \frac{\|A^{-1}\|}{1-\|A^{-1}\|\|E\|}.
```

Como consecuencia del último lema se tiene el siguiente corolario.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000016
md"""
**Corolario 6.22.** Sea ``C:D\subset\mathbb{R}^n\to\mathbb{R}^{n\times n}`` continua en ``x_*\in D``. Si ``C(x_*)`` es invertible, entonces existe ``\delta>0`` tal que para todo ``x\in B(x_*,\delta)``, ``C(x)`` es invertible y más aún, ``\|C(x)^{-1}\|\leq2\|C(x_*)^{-1}\|``.

*Demostración.* Notemos que para cada ``x\in D``, la matriz ``C(x)`` puede verse como una perturbación de ``C(x_*)``, esto es,
```math
C(x) = \underbrace{C(x_*)}_{=:A}+\underbrace{(C(x)-C(x_*))}_{=:E}.
```
Puesto que ``C`` es continua en ``x_*``, sabemos que existe ``\delta>0`` tal que
```math
\|C(x)-C(x_*)\| < \frac{1}{2\|C(x_*)^{-1}\|}, \qquad \text{si } x\in B(x_*,\delta).
```
Por el lema anterior, tenemos que ``C(x)`` es invertible para todo ``x\in B(x_*,\delta)`` y que
```math
\|C(x)^{-1}\| \leq \frac{\|C(x_*)^{-1}\|}{1-\|C(x_*)^{-1}\|\|C(x)-C(x_*)\|} < \frac{\|C(x_*)^{-1}\|}{1-\|C(x_*)^{-1}\|\frac{1}{2\|C(x_*)^{-1}\|}} = 2\|C(x_*)^{-1}\|.
```

Ahora estamos en condiciones de establecer cuándo la iteración de Newton está bien definida.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000017
md"""
**Teorema 6.23** (Buena definición del método de Newton y primer resultado de convergencia local). Sea ``D\subset\mathbb{R}^n`` abierto y supongamos que ``F:D\to\mathbb{R}^n`` es diferenciable en ``D``. Sea ``x_*\in D`` un cero de ``F`` (i.e., ``F(x_*)=0``). Si ``F'`` es continua en ``x_*`` y ``F'(x_*)`` es invertible entonces se cumple lo siguiente:

- Buena definición: ``G(x):=x-(F'(x))^{-1}F(x)`` está bien definida en un entorno de ``x_*``.

- Convergencia local: La iteración
```math
x_{k+1} = x_k-(F'(x_k))^{-1}F(x_k), \qquad k=0,1,2,\dots
```
converge localmente a ``x_*``.

- Velocidad de convergencia superlineal: Si ``x_0`` se elige suficientemente cerca de ``x_*`` entonces
```math
\lim_{k\to\infty}\frac{\|x_{k+1}-x_*\|}{\|x_k-x_*\|} = 0.
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000018
md"""
*Demostración.*

- Tomando ``C=F'`` en el Corolario 6.22 tenemos la buena definición del método de Newton.

- Puesto que ``F(x_*)=0``, se tiene que ``x_*=G(x_*)``. Además, se puede demostrar (ver Ejercicio 6.8) que ``G`` es diferenciable en ``x_*`` y que ``G'(x_*)=0``. Por el Teorema de Ostrowski (Teorema 6.8), ya que ``\rho(G'(x_*))=0<1``, se tiene que ``x_*`` es un punto de atracción de ``G``.

- Puesto que ``G`` es diferenciable en ``x_*`` tenemos que
```math
\lim_{x\to x_*}\frac{\|G(x)-G(x_*)-G'(x_*)(x-x_*)\|}{\|x-x_*\|} = 0.
```
Considerando que ``G(x_*)=x_*`` y que ``G'(x_*)=0``, la igualdad anterior puede reescribirse como
```math
\lim_{x\to x_*}\frac{\|G(x)-x_*\|}{\|x-x_*\|} = 0.
```
Puesto que ``x_*`` es un punto de atracción de ``G``, podemos tomar ``x_0`` suficientemente cerca de ``x_*`` de modo que
```math
x_{k+1} := G(x_k)\to x_*, \qquad \text{cuando } k\to\infty.
```
Por lo tanto,
```math
\lim_{k\to\infty}\frac{\|x_{k+1}-x_*\|}{\|x_k-x_*\|} = 0,
```
lo que concluye la demostración.

Si en el teorema anterior reemplazamos la hipótesis de continuidad de ``F'`` en ``x_*``, por una hipótesis un poco más fuerte (ver (6.5)) tenemos que el orden de convergencia de la iteración de Newton es cuadrático, como se muestra en el siguiente resultado.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000019
md"""
**Teorema 6.24** (Segundo resultado de convergencia local: Convergencia cuadrática). Sea ``D\subset\mathbb{R}^n`` abierto y supongamos que ``F:D\to\mathbb{R}^n`` es diferenciable en ``D``. Sea ``x_*\in D`` un cero de ``F`` (i.e., ``F(x_*)=0``). Supongamos que existen ``\delta>0`` y ``L>0`` tales que
```math
\|F'(x)-F'(x_*)\| \leq L\|x-x_*\|, \qquad \forall x\in B(x_*,\delta). \tag{6.5}
```
Entonces, si ``F'(x_*)`` es invertible, existe ``c>0`` tal que si ``x_0`` se elige suficientemente cerca de ``x_*``, la iteración
```math
x_{k+1} = x_k-(F'(x_k))^{-1}F(x_k), \qquad k=0,1,2,\dots
```
converge cuadráticamente, es decir,
```math
\|x_{k+1}-x_*\| \leq c\|x_k-x_*\|^2, \qquad k=0,1,2,\dots
```

**Observación 6.25** (Condiciones suficientes para que se cumpla (6.5)). Si ``F'`` es Lipschitz continua en un entorno de ``x_*`` entonces en particular se cumple (6.5). Por otro lado, ``F'`` es Lipschitz continua si y sólo si ``\frac{\partial f_i}{\partial x_j}`` es Lipschitz continua, para ``i,j=1,2,\dots,n``. Notemos además que para esto es suficiente que ``\frac{\partial^2f_i}{\partial x_j\partial x_k}`` estén acotadas en un entorno de ``x_*``, para ``i,j,k=1,2,\dots,n``. Finalmente, notemos también que para que se cumpla la última condición, es suficiente que ``F`` sea dos veces continuamente diferenciable en un entorno de ``x_*``, es decir, que ``\frac{\partial^2f_i}{\partial x_j\partial x_k}`` sean continuas en un entorno de ``x_*``, para ``i,j,k=1,2,\dots,n``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000001a
md"""
*Demostración del Teorema 6.24.* Veamos primero que (6.5) implica que
```math
\|F(x)-F(x_*)-F'(x_*)(x-x_*)\| \leq \frac{L}{2}\|x-x_*\|^2, \qquad \forall x\in B(x_*,\delta). \tag{6.6}
```
Por el Teorema 6.16 se tiene que
```math
F(x)-F(x_*) = \int_0^1F'(x_*+t(x-x_*))(x-x_*)\,dt, \qquad \forall x\in B(x_*,\delta).
```
En consecuencia,
```math
\begin{align*}
\|F(x)-F(x_*)-F'(x_*)(x-x_*)\| &= \left\|\int_0^1(F'(x_*+t(x-x_*))-F'(x_*))(x-x_*)\,dt\right\|\\
&\leq \int_0^1\|(F'(x_*+t(x-x_*))-F'(x_*))(x-x_*)\|\,dt\\
&\leq \int_0^1\|F'(x_*+t(x-x_*))-F'(x_*)\|\|x-x_*\|\,dt\\
&\leq \int_0^1L\|t(x-x_*)\|\|x-x_*\|\,dt\\
&= L\|x-x_*\|^2\int_0^1t\,dt\\
&= \frac{L}{2}\|x-x_*\|^2,
\end{align*}
```
para todo ``x\in B(x_*,\delta)``, y así, (6.6) queda demostrada.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000001b
md"""
Ahora bien, si ``x\in B(x_*,\delta)``, tenemos que
```math
\begin{align*}
\|G(x)-G(x_*)\| &= \|x-F'(x)^{-1}F(x)-x_*\|\\
&= \|-F'(x)^{-1}[F(x)-\underbrace{F(x_*)}_{=0}-F'(x)(x-x_*)]\|\\
&= \|-F'(x)^{-1}[F(x)-F(x_*)-F'(x_*)(x-x_*)-(F'(x)-F'(x_*))(x-x_*)]\|\\
&\leq \|F'(x)^{-1}\|\Big[\|F(x)-F(x_*)-F'(x_*)(x-x_*)\|+\|F'(x)-F'(x_*)\|\|x-x_*\|\Big]\\
&\leq \|F'(x)^{-1}\|\Big[\frac{L}{2}\|x-x_*\|^2+L\|x-x_*\|^2\Big]\\
&= \|F'(x)^{-1}\|\frac{3L}{2}\|x-x_*\|^2.
\end{align*}
```
Puesto que ``C:=F'`` es continua en ``x_*`` y que ``C(x_*)=F'(x_*)`` es invertible, por el Corolario 6.22, tenemos que existe ``\delta'>0`` tal que ``C(x)=F'(x)`` es invertible para ``x\in B(x_*,\delta')`` y
```math
\|F'(x)^{-1}\| \leq 2\|F'(x_*)^{-1}\|, \qquad \forall x\in B(x_*,\delta').
```
Por lo tanto,
```math
\|G(x)-G(x_*)\| \leq \underbrace{3L\|F'(x_*)^{-1}\|}_{=:c}\|x-x_*\|^2, \qquad \forall x\in B(x_*,\delta)\cap B(x_*,\delta').
```
Finalmente, puesto que ``x_*`` es un punto de atracción de ``G``, podemos elegir ``\delta''>0`` con ``\delta''\leq\min\{\delta,\delta'\}`` tal que para todo ``x_0\in B(x_*,\delta'')``, ``x_{k+1}=G(x_k)\in B(x_*,\delta'')`` y así, ``\|x_{k+1}-x_*\|=\|G(x_k)-G(x_*)\|\leq c\|x_k-x_*\|^2``, ``k=0,1,2,\dots``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000001c
md"""
### 6.2.1 Convergencia global del método de Newton

En esta sección nos proponemos estudiar en qué casos el método de Newton converge partiendo de cualquier valor inicial.

**Observación 6.26** (Convergencia global de métodos de punto fijo). Si ``G:\mathbb{R}^n\to\mathbb{R}^n`` es una contracción, entonces por el Teorema 6.12 tenemos que partiendo de cualquier ``x_0\in\mathbb{R}^n``, la iteración
```math
x_{k+1} := G(x_k), \qquad k=0,1,2,\dots
```
converge a un punto fijo ``x_*`` de ``G`` (que existe y es único).

La observación anterior nos da un resultado de convergencia global, pero es difícil que la hipótesis se cumpla en problemas concretos. A continuación veremos un resultado de convergencia global para el método de Newton, bajo ciertas hipótesis sobre la función ``F``, que incluyen el concepto de *convexidad*.

**Definición 6.27** (Función convexa). Decimos que una función ``f:\mathbb{R}^n\to\mathbb{R}`` es *convexa* si
```math
f(tx+(1-t)y) \leq tf(x)+(1-t)f(y), \qquad \forall t\in[0,1],\ x,y\in\mathbb{R}^n.
```
Además, una función ``F:\mathbb{R}^n\to\mathbb{R}^n`` con ``F=\begin{pmatrix}f_1\\f_2\\\vdots\\f_n\end{pmatrix}`` es *convexa* si ``f_i`` es convexa para todo ``i=1,2,\dots,n``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000001d
md"""
**Proposición 6.28** (Caracterización de las funciones convexas). Sea ``f:\mathbb{R}^n\to\mathbb{R}`` diferenciable. Entonces, ``f`` es convexa si y sólo si
```math
f(y)-f(x) \geq \nabla f(x)^T(y-x), \qquad \forall x,y\in\mathbb{R}^n.
```

Si ``n=1``, tenemos que ``\nabla f(x)=f'(x)``, y la caracterización dice que una función es convexa si y sólo si el gráfico de ``f`` está por arriba de su recta tangente en todos los puntos, es decir,
```math
f(y) \geq f(x)+f'(x)(y-x), \qquad \forall x,y\in\mathbb{R}.
```

*Demostración.* (``\Rightarrow``) Sean ``x,y\in\mathbb{R}^n``, y sea ``0<h<1``. Puesto que ``f`` es convexa, tenemos que
```math
f((1-h)x+hy) \leq (1-h)f(x)+hf(y).
```
Teniendo en cuenta que ``(1-h)x+hy=x+h(y-x)``, obtenemos
```math
\frac{f(x+h(y-x))-f(x)}{h} \leq f(y)-f(x),
```
y tomando el límite para ``h\to0``,
```math
\nabla f(x)^T(y-x) \leq f(y)-f(x),
```
como queríamos demostrar.

(``\Leftarrow``) Sean ``x,y\in\mathbb{R}^n``, y sea ``0<t<1``. Queremos demostrar que
```math
f(tx+(1-t)y) \leq tf(x)+(1-t)f(y). \tag{6.7}
```
Sea ``z:=tx+(1-t)y``. Por hipótesis,
```math
f(x)-f(z) \geq \nabla f(z)^T(x-z), \tag{6.8}
```
y
```math
f(y)-f(z) \geq \nabla f(z)^T(y-z). \tag{6.9}
```
Multiplicando (6.8) por ``t`` y (6.9) por ``1-t`` y sumando miembro a miembro, obtenemos
```math
tf(x)+(1-t)f(y)-f(z) \geq \nabla f(z)^T(tx-tz+(1-t)y-z+tz) = \nabla f(z)^T0,
```
de donde se sigue (6.7) y por lo tanto, que ``f`` es convexa.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000001e
md"""
**Teorema 6.29** (Teorema de Newton-Baluev: Convergencia global del método de Newton). Sea ``F:\mathbb{R}^n\to\mathbb{R}^n`` convexa. Supongamos que ``F`` es continuamente diferenciable en ``\mathbb{R}^n``, que ``F'(x)`` es invertible y que las entradas de ``F'(x)^{-1}`` son no negativas, para todo ``x\in\mathbb{R}^n``. Entonces, si existe ``x_*\in\mathbb{R}^n`` tal que ``F(x_*)=0``, este cero es único y la iteración
```math
x_{k+1} := x_k-F'(x_k)^{-1}F(x_k), \qquad k=0,1,2,\dots
```
converge a ``x_*``, partiendo de cualquier valor inicial ``x_0\in\mathbb{R}^n``. Más aún, si ``x_*=(x_*^{(1)},x_*^{(2)},\dots,x_*^{(n)})^T`` y ``x_k=(x_k^{(1)},x_k^{(2)},\dots,x_k^{(n)})^T``, para ``k\in\mathbb{N}``, entonces
```math
x_*^{(i)} \leq \dots \leq x_{k+1}^{(i)} \leq x_k^{(i)} \leq \dots \leq x_1^{(i)}, \qquad i=1,2,\dots,n.
```

**Observación 6.30.** Notemos que si ``n=1``, es decir, si ``F:\mathbb{R}\to\mathbb{R}`` entonces la matriz jacobiana de ``F`` en ``x`` es la derivada de ``F`` en ``x``, y las hipótesis del teorema anterior equivalen a pedir que ``F\in\mathcal{C}^1``, ``F`` convexa y que ``F'(x)>0`` para todo ``x\in\mathbb{R}``. En este caso, si existe un cero ``x_*`` de ``F``, éste es único y la sucesión del método de Newton converge decrecientemente a ``x_*``, partiendo de cualquier ``x_0\in\mathbb{R}``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000001f
md"""
*Demostración del Teorema 6.29.* En esta demostración utilizaremos la notación ``x\geq0`` para un vector ``x=(x^{(1)},x^{(2)},\dots,x^{(n)})^T\in\mathbb{R}^n`` si satisface
```math
x^{(i)} \geq 0, \qquad \forall i=1,2,\dots,n.
```

Sea ``x_0\in\mathbb{R}^n`` arbitrario, y consideremos la iteración de Newton
```math
x_{k+1} := x_k-F'(x_k)^{-1}F(x_k), \qquad k=0,1,2,\dots \tag{6.10}
```
Sea ``F=\begin{pmatrix}f_1\\f_2\\\vdots\\f_n\end{pmatrix}``. Puesto que ``F`` es convexa, tenemos que ``f_i`` es convexa, para ``i=1,2,\dots,n``. Para cada ``k\in\mathbb{N}_0``, de la Proposición 6.28, y teniendo en cuenta (6.10), se sigue que
```math
f_i(x_{k+1})-f_i(x_k) \geq \nabla f_i(x_k)^T(x_{k+1}-x_k) = -f_i(x_k), \qquad i=1,2,\dots,n,
```
y así,
```math
f_i(x_{k+1}) \geq 0, \qquad i=1,2,\dots,n.
```
Considerando la notación introducida al principio de esta demostración, la última desigualdad implica que
```math
F(x_k) \geq 0, \qquad k=1,2,\dots
```

Ahora vamos a demostrar que ``\{x_k^{(i)}\}_{k\in\mathbb{N}}`` es una sucesión decreciente, para cada ``i=1,2,\dots,n``. En efecto, puesto que las entradas de ``F'(x_k)^{-1}`` son no negativas y que ``F(x_k)\geq0``, de (6.10) se sigue que
```math
x_k-x_{k+1} = F'(x_k)^{-1}F(x_k) \geq 0,
```
lo que significa que
```math
x_{k+1}^{(i)} \leq x_k^{(i)}, \qquad i=1,2,\dots,n.
```
Por otro lado, veamos que ``\{x_k^{(i)}\}_{k\in\mathbb{N}}`` está acotada inferiormente por ``x_*^{(i)}``, para cada ``i=1,2,\dots,n``; es decir, veamos que para ``k=1,2,\dots`` se tiene que
```math
x_*^{(i)} \leq x_k^{(i)}, \qquad i=1,2,\dots,n. \tag{6.11}
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000020
md"""
Sea ``k\in\mathbb{N}``. Haciendo uso de la Proposición 6.28 nuevamente, tenemos que
```math
\underbrace{f_i(x_*)}_{=0}-f_i(x_k) \geq \nabla f_i(x_k)^T(x_*-x_k), \qquad i=1,2,\dots,n,
```
o equivalentemente,
```math
-f_i(x_k)-\nabla f_i(x_k)^T(x_*-x_k) \geq 0, \qquad i=1,2,\dots,n.
```
La última desigualdad puede reescribirse en la nueva notación como
```math
-F(x_k)-F'(x_k)(x_*-x_k) \geq 0.
```
Ya que las entradas de ``F'(x_k)^{-1}`` son no negativas, tenemos que
```math
F'(x_k)^{-1}[-F(x_k)-F'(x_k)(x_*-x_k)] \geq 0.
```
En consecuencia,
```math
-F'(x_k)^{-1}F(x_k)+(x_k-x_*) \geq 0,
```
y luego
```math
x_k-x_* \geq \underbrace{F'(x_k)^{-1}\overbrace{F(x_k)}^{\geq0}}_{\geq0},
```
y finalmente,
```math
x_k-x_* \geq 0,
```
de donde se sigue (6.11). Puesto que toda sucesión monótona decreciente y acotada inferiormente converge, tenemos que existe ``y^{(i)}\in\mathbb{R}`` tal que
```math
\lim_{k\to\infty}x_k^{(i)} = y^{(i)}, \qquad i=1,2,\dots,n.
```
Definiendo ``y:=\begin{pmatrix}y^{(1)}\\y^{(2)}\\\vdots\\y^{(n)}\end{pmatrix}``, tenemos entonces que
```math
\lim_{k\to\infty}x_k = y.
```
Sólo resta ver que ``y`` es un cero de ``F`` y que ``F`` tiene un único cero (en particular, ``y=x_*``). Utilizando la continuidad de ``F`` y de ``F'``, y la definición de la iteración (6.10) tenemos que
```math
F(y) = \lim_{k\to\infty}F(x_k) = \lim_{k\to\infty}F'(x_k)(x_k-x_{k+1}) = F'(y)(y-y) = 0.
```
Así, ``y`` es un cero de ``F``. Se deja como ejercicio demostrar que ``F`` tiene un único cero. Como ayuda, suponer que existen dos ceros ``x_*`` y ``y_*`` y utilizar la condición de convexidad de ``F`` (aplicada a ``x_*`` e ``y_*``) dada por la Proposición 6.28.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000021
md"""
## 6.3 Problemas
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000022
md"""
### 6.1)

Consideremos la iteración de punto fijo asociada a una función ``G:\mathbb{R}\to\mathbb{R}``.

#### 6.1.a)
Si ``G(x):=x-x^3``, demostrar que aunque ``G'(0)=1\ (\not<1)``, el punto 0 es un punto de atracción de ``G``.

#### 6.1.b)
Si ``G(x):=x+x^3``, demostrar que ``G'(0)=1``, y que el punto 0 *no* es un punto de atracción de ``G``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000023
md"""
### 6.2)

Sea ``G:D\subset\mathbb{R}^n\to\mathbb{R}^n``. Sea ``x_*\in D`` un punto fijo de ``G``. Demostrar que si ``G`` es continuamente diferenciable en un entorno de ``x_*`` y ``\rho(G'(x_*))<1``, entonces ``G`` es una contracción en alguna bola cerrada centrada en ``x_*``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000024
md"""
### 6.3)

Supongamos que ``G:\mathbb{R}^n\to\mathbb{R}^n`` tiene la propiedad que para todo conjunto cerrado y acotado ``C\subset\mathbb{R}^n`` existe una constante ``\alpha_C<1`` tal que
```math
\|G(x)-G(y)\| \leq \alpha_C\|x-y\|, \qquad \forall x,y\in C.
```
Supongamos además que ``G`` tiene un punto fijo ``x_*``. Mostrar que ``x_*`` es único y que las iteraciones ``x_{k+1}=G(x_k)``, ``k=0,1,\dots`` convergen a ``x_*`` para cualquier ``x_0``. Dar un ejemplo en ``\mathbb{R}^1`` que muestre que la condición sobre ``G`` no es suficiente para garantizar la existencia de un punto fijo.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000025
md"""
### 6.4)

Considerar el sistema
```math
x = \frac{\frac12}{1+(x+y)^2} \qquad y = \frac{\frac12}{1+(x-y)^2}.
```
Encontrar un conjunto acotado ``D\subset\mathbb{R}^2`` donde se pueda aplicar el Teorema 6.12.

*Ayuda: ¿Cuál es el signo de las componentes de la solución del sistema? ¿Cuáles son los máximos valores posibles para dichas componentes según las ecuaciones del sistema?*
"""

# ╔═╡ f4a16bce-53e9-4382-8dcc-0774179f747a
@warn "Hecho hasta el miercoles"

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000026
md"""
### 6.5)

Considerar el sistema
```math
x = 1+h\frac{e^{-x^2}}{1+y^2} \qquad y = \frac12+h\tan^{-1}(x^2+y^2).
```
Mostrar que si ``h`` se elige suficientemente chico, entonces este sistema tiene una única solución ``x_*`` dentro de alguna región rectangular. Más aún, mostrar que la iteración de punto fijo (propuesta en la forma que viene dado el sistema) converge a dicha solución.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000027
md"""
### 6.6)

El Teorema 6.16 puede interpretarse también como la generalización del Teorema del Valor Medio del Cálculo Diferencial a funciones de valores vectoriales. En el siguiente ejemplo se muestra que una generalización, quizás más natural que la del mencionado teorema, no es posible.

Sea ``G:\mathbb{R}^2\to\mathbb{R}^2`` definida por
```math
G\begin{pmatrix}x_1\\x_2\end{pmatrix} = \begin{pmatrix}x_1^3\\x_1^2\end{pmatrix}.
```
Sean ``x=(0,0)^T`` e ``y=(1,1)^T``. Demostrar que no existe ningún ``z\in\mathbb{R}^2`` tal que
```math
G(x)-G(y) = G'(z)(x-y).
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000028
md"""
### 6.7)

Demostrar el Teorema 6.17.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000029
md"""
### 6.8)

Sea ``D\subset\mathbb{R}^n`` abierto. Sea ``F:D\to\mathbb{R}^n`` diferenciable en ``D``. Sea ``x_*\in D`` tal que ``F(x_*)=0``. En la primera parte del Teorema 6.23 vimos que si ``F'`` es continua en ``x_*`` y ``F'(x_*)`` es invertible, entonces
```math
G(x) := x-F'(x)^{-1}F(x),
```
está bien definida en un entorno de ``x_*``. Demostrar que, además, ``G`` es diferenciable en ``x_*`` y que ``G'(x_*)=0``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000002a
md"""
### 6.9)

Considerar un sistema no lineal
```math
F(x) = 0,
```
donde ``F:D\subset\mathbb{R}^n\to\mathbb{R}^n``. Escribir una función `[x,iter,X]=newton(F, DF, x0, miter,tol)` para resolver el sistema por el método de Newton. Los parámetros de entrada son una referencia `F` a la función ``F(x)``, una referencia `DF` a la matriz jacobiana ``F'(x)``, un punto inicial `x0`, el máximo número de iteraciones permitido `miter` y una tolerancia `tol` permitida. Los parámetros de salida son la aproximación `x` de una solución del sistema, la cantidad de iteraciones realizadas `iter` y la matriz `X` cuyas columnas son las sucesivas aproximaciones que se calcularon.

*Sugerencia:* Para resolver el sistema lineal necesario en cada iteración utilizar el comando "``\backslash``" de MATLAB/OCTAVE.
"""

# ╔═╡ 10e8fced-4b9f-4544-935b-616a446bec68
"""
	newton(F, DF, x0, miter, tol)
"""
function newton(F, DF, x0, miter, tol)
	Ⅎ(v) = DF(v...) \ F(v...)
	X = Matrix{Float64}(undef, length(x0), miter)
	
	X[:, 1] = x0
	k = 1
	
	while k < miter 
		y = Ⅎ(X[:, k])
		X[:, k+1] = X[:, k] - y
		k += 1
		
		norm(y) < tol && break
	end

	return X[:,k], k-1, X[:, 1:k]
end

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000002b
md"""
### 6.10)

Estudiar la convergencia y el orden experimental de la misma (si corresponde) al aplicar el método de Newton al sistema
```math
x^2+y^2 = 4 \qquad x^2-y^2 = 1,
```
utilizando como valor inicial ``\begin{pmatrix}1.6\\1.2\end{pmatrix}``.
"""

# ╔═╡ 5a291aef-d9b7-4c18-838e-a8d828dd092c
let
	f1(x, y) = x^2 + y^2 - 4
	f2(x, y) = x^2 - y^2 -1
	x0 = [0.5, 0.5]
	
	F(x, y) = [f1(x, y), f2(x, y)]
	DF(x, y) = [2x 2y; 2x -2y]
	
	x = newton(F, DF, x0, 100, 0.0001)
	sol = x[1]


	fig = Figure(size = (700,400))

	ax1 = Axis(fig[1,1], aspect = DataAspect())

	xs = LinRange(-2, 2, 100)
	ys = LinRange(-2, 2, 100)
	zs1 = [f1(x, y) for x in xs, y in ys]
	zs2 = [f2(x, y) for x in xs, y in ys]

	contour!(ax1, xs, ys, zs1, levels = [0], linewidth=2, color = :firebrick)
	contour!(ax1, xs, ys, zs2, levels = [0], linewidth=2, color = :steelblue)

	xv1 = [x[3][:,i][1] for i in 1:x[2]]
	yv1 = [x[3][:,i][2] for i in 1:x[2]]

	
	scatter!(ax1, Point(sol...))
	scatterlines!(ax1, xv1, yv1, linewidth=1, markersize=5, color = :black)
	
	fig
end

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000002c
md"""
### 6.11)

#### 6.11.a)
El siguiente sistema tiene cuatro ceros en el dominio ``(-4,4)\times(-4,4)``
```math
f(x,y) = x^2+xy^3-9 = 0, \qquad g(x,y) = 3x^2y-y^3-4 = 0
```
Utilizar `meshgrid` y `contour` para graficar las curvas de nivel cero de ``f`` y ``g`` en el mismo dibujo y así determinar puntos iniciales para hallar cada una de las cuatro soluciones.

#### 6.11.b)
Utilizar la función `newton` para calcular las cuatro soluciones del sistema del item anterior. ¿Cuál fue el número de iteraciones realizadas en cada caso? Estimar la velocidad de convergencia.

#### 6.11.c)
Definir la función ``h(x,y)=f(x,y)^2+g(x,y)^2``. Luego, los ceros del sistema corresponden a los mínimos de ``h``. Utilizar la función `fminsearch` de MATLAB/OCTAVE para aproximar los cuatro ceros.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000002d
md"""
### 6.12)

Utilizar el método de Newton para hallar todas las raíces de cada uno de los siguientes sistemas no lineales. Utilizar gráficos para escoger valores iniciales.

#### 6.12.a)
```math
x^2+y^2-2x-2y+1 = 0 \qquad x+y-2xy = 0.
```

#### 6.12.b)
```math
x^2+2xy+y^2-x+y-4 = 0 \qquad 5x^2-6xy+5y^2+16x-16y+12 = 0.
```
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000002e
md"""
### 6.13)

Usar el método de Newton para hallar el mínimo valor de
```math
f(x,y) = x^4+xy+(1+y)^2.
```
Experimentar con varios valores iniciales y observar los patrones de convergencia.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-00000000002f
md"""
### 6.14)

Considerar la iteración de Newton en una dimensión para la función ``F(x)=x^2``.

#### 6.14.a)
Demostrar que la iteración converge localmente a ``x_*=0`` pero que la rapidez de convergencia no es cuadrática.

#### 6.14.b)
Mostrar qué hipótesis del Teorema 6.24 no se cumple.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000030
md"""
### 6.15)

Repetir el ejercicio anterior considerando la función ``F(x)=x+x^{1+\alpha}``, donde ``0<\alpha<1``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000031
md"""
### 6.16)

La iteración
```math
x_{k+1} = x_k-F'(x_0)^{-1}F(x_k), \qquad k=0,1,\dots
```
se conoce como *iteración simplificada de Newton*. Demostrar que si ``F`` es diferenciable en un entorno de un punto ``x_*`` tal que ``F(x_*)=0``, si ``F'`` es continua en ``x_*``, y ``F'(x_*)`` es invertible, entonces esta iteración converge localmente a ``x_*``.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000032
md"""
### 6.17)

Sea ``h>0`` un parámetro fijo. Sea ``F:\mathbb{R}^n\to\mathbb{R}^n`` dada por ``F(x)=Ax+\Phi(x)`` donde
```math
A = \begin{bmatrix}
2&-1&&&&\\
-1&2&-1&&0&\\
&-1&2&-1&&\\
&&\ddots&\ddots&\ddots&\\
&0&&-1&2&-1\\
&&&&-1&2
\end{bmatrix}, \qquad
\Phi(x) = h^2\begin{bmatrix}f(x^{(1)})\\f(x^{(2)})\\f(x^{(3)})\\\vdots\\f(x^{(n-1)})\\f(x^{(n)})\end{bmatrix},
```
con ``x=(x^{(1)},x^{(2)},x^{(3)},\dots,x^{(n-1)},x^{(n)})^T`` y ``f(t)=t^3``.

#### 6.17.a)
Usando que ``A`` es invertible, demostrar que el método de Newton aplicado al sistema ``F(x)=0`` converge localmente a ``x_*=0``.

#### 6.17.b)
Estudiar el orden de convergencia del método de Newton para este caso.
"""

# ╔═╡ f6a7b8c9-eeee-4a1b-8c2d-000000000033
md"""
### 6.18)

Considerar las funciones ``F:\mathbb{R}\to\mathbb{R}`` dadas por ``F(x)=\tan^{-1}x`` y ``F(x)=2x+\text{sen}\,x``. Mostrar que para ambas funciones se satisfacen todas las hipótesis del Teorema 6.29 (de Newton-Baluev), excepto la convexidad, y que la iteración de Newton correspondiente no converge globalmente.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CairoMakie = "~0.15.15"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "87a04fcd9461c40396f39e6a5f41080733bb0526"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "dee8388a9bec89c8713ad176f96960ae2f98598a"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.4"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
registries = "General"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
registries = "General"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
registries = "General"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
registries = "General"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
registries = "General"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "1cda0b7d5abfc95357dae18aca934d401f7869ad"
registries = "General"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "7b841680738c19948120f6e4cf8d0200518bd564"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.8+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
registries = "General"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
registries = "General"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "61761f58648aa7217445f24f841839b78c712232"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.3"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "291665b547f137df070e4dd83e432b5fee8cc4a0"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.2"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "7bc84b769c1d384315e7b5c4ac03a6c303e6cf35"
registries = "General"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.8"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
registries = "General"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "4ac548adcad90c1d5d677af13568a748af4c952b"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a958ab3a40c755563f5e1405c0846cb0446bf19d"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.131"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
registries = "General"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
registries = "General"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2bfb1e047e2ad0a5ca94365340bde8005d637568"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.4+0"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d9d3cd382f2c999f684e6bfa4039fcc16c153b6f"
registries = "General"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "9.0.2+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
registries = "General"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
registries = "General"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
registries = "General"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
registries = "General"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "086b5fbd032baf544678cc15b52b015e4c4aceb8"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.1"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
registries = "General"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
registries = "General"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
registries = "General"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "ec46c5825710fa1a15d468acb2d93cc939a7a5fe"
registries = "General"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.13"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
registries = "General"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3efbbc027441271444dcd0c0a46f2d119dc4329"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "6.1.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "090526e65de8f69648ac156daae153de8b56df62"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.88.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "ef70da5e123a06a29e2d6ddff0f09985bc226491"
registries = "General"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "9d9531a9cb63a9edc33836414e82a07e81710de2"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.14004.0+0"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "31bb6c92405c084617facc1d7ed9eb6c402d061e"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.30"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "f0f005f997dfb8c5fe23920d99458a9619873893"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.10"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "1c531bf0f8a5c60a340926e058fd3f209b5eef5d"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.12"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticMakieExt = "Makie"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "0db6aea5b64caa1c47e9fdd27394f0296f05e8bf"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.15"

    [deps.IntervalSets.extensions]
    IntervalSetsMakieExt = "Makie"
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
registries = "General"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "4657a834b01ce00f2b94546095c4350617b9af2c"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.10.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
registries = "General"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
registries = "General"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39bca05343661c347aae0bca57a5994a0bf4f08d"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.2.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "b85e2797b2409570e84c4de46238c0ed5f6476ae"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.2"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "5f6f5d1b1fb7ff98c9a083bbfd4c661a9808e758"
registries = "General"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.15"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
registries = "General"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
registries = "General"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
registries = "General"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e9ae72527609bb66882d38ad009ffd47ca5713fb"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.16+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
registries = "General"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
registries = "General"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
registries = "General"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "663e8b48b789916221e0765393b289ca6c88f24e"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "ddcc36fd83fb4cd459a53c891981697bde0d3de6"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
registries = "General"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "5005266de4bfe50e53ff44a5cb5c540b6e47a254"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
registries = "General"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "e6ffa7bd4f8d0b8744d18f219304738d952f9a37"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.11"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "57aa595158717ef165e6f5ab639fe2e3178c0a2b"
registries = "General"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
registries = "General"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2c79185e5261b159474903598c571ddd9f12ff7b"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "429071b23f4c9a13fb6582f807cc2ef454082408"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "de349d5585559b692f518b6d3fc7c179d64b2364"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.23"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "ad8002667372439f2e3611cfd14097e03fa4bccd"
registries = "General"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.3"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "b814d5005d6a529d740ffe06f8a86396f6501138"
registries = "General"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.9.2"

    [deps.StructUtils.extensions]
    StructUtilsLazilyInitializedFieldsExt = ["LazilyInitializedFields"]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    LazilyInitializedFields = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
registries = "General"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
registries = "General"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
registries = "General"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
registries = "General"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
registries = "General"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1210ba774d3427387d307bf1f416d699b7c39417"
registries = "General"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.15.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cb007192783c56d8249db4cf0e3495001edfe414"
registries = "General"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.5+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "28e57478e8a160d346a19c28b3fffb9273bcc9c2"
registries = "General"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.134+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
registries = "General"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1fcd6cf05d23adab0b5ffc496d020537215751dc"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "e067c8bae65bb40866552296a2d5b4dc65b8928f"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.80.702+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
registries = "General"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
registries = "General"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "52d3b9475133c3bc8c0a7f90f18bfc8cdb5a443a"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
registries = "General"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
registries = "General"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═f6a7b8c9-eeee-4a1b-8c2d-000000000001
# ╠═f6a7b8c9-eeee-4a1b-8c2d-000000000002
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000003
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000004
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000005
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000006
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000007
# ╟─b8cd7866-eaba-4c17-901d-e7ce33fa3ade
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000008
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000009
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000000a
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000000b
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000000c
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000000d
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000000e
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000000f
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000010
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000011
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000012
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000013
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000014
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000015
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000016
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000017
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000018
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000019
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000001a
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000001b
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000001c
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000001d
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000001e
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000001f
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000020
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000021
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000022
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000023
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000024
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000025
# ╟─f4a16bce-53e9-4382-8dcc-0774179f747a
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000026
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000027
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000028
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000029
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000002a
# ╠═10e8fced-4b9f-4544-935b-616a446bec68
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000002b
# ╠═5a291aef-d9b7-4c18-838e-a8d828dd092c
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000002c
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000002d
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000002e
# ╟─f6a7b8c9-eeee-4a1b-8c2d-00000000002f
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000030
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000031
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000032
# ╟─f6a7b8c9-eeee-4a1b-8c2d-000000000033
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
