### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000001
begin 
	using LinearAlgebra, PlutoUI, CairoMakie
	import Plots
end

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000002
TableOfContents()

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000003
md"""
# Capítulo 5: Métodos de descenso para sistemas lineales
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000004
md"""
## 5.1 Sistemas sdp y minimización de un funcional

En esta sección estudiaremos métodos iterativos para resolver sistemas de ecuaciones lineales asociados a matrices simétricas y definidas positivas. Estos métodos están basados en el hecho de que la solución de un sistema asociado a una matriz ``A`` sdp está caracterizada por ser el punto donde cierto funcional de energía asociado a ``A`` alcanza su mínimo.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000005
md"""
**Teorema 5.1.** *Sean ``A`` sdp y ``b\in\mathbb{R}^n``. Sea ``\phi:\mathbb{R}^n\to\mathbb{R}`` dada por ``\phi(x)=\frac12x^TAx-x^Tb``. Entonces, el mínimo de ``\phi`` en ``\mathbb{R}^n`` se alcanza solo en la solución ``x_*`` de ``Ax=b`` y el valor mínimo de ``\phi`` es ``-\frac{b^Tx_*}{2}``.*

**Demostración.** Sea ``x_*:=A^{-1}b`` y ``y\in\mathbb{R}^n`` arbitrario, ``y\neq x_*``. Vamos a probar que ``\phi(x_*)<\phi(y)``. Llamando ``z:=y-x_*`` tenemos que

"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000006
md"""
```math
\begin{align*}
\phi(y) &= \phi(x_*+z) = \frac12(x_*+z)^TA(x_*+z)-(x_*+z)^Tb\\
&= \frac12x_*^TAx_*+x_*^TAz+\frac12z^TAz-x_*^Tb-z^Tb\\
&= \phi(x_*)+\frac12z^TAz.
\end{align*}
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000007
md"""
Puesto que ``A`` es sdp, tenemos que ``\frac12z^TAz>0``, y así,

```math
\phi(x_*)<\phi(y), \qquad \forall y\in\mathbb{R}^n,\ y\neq x_*.
```

Por lo tanto, el mínimo de ``\phi`` se alcanza solo en ``x_*`` y

```math
\phi(x_*) = \frac12x_*^TAx_*-x_*^Tb = -\frac{x_*^Tb}{2}.
```

"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000008
md"""
*Verificación.* Escribiendo ``\phi`` en coordenadas,
```math
\phi(x) = \frac12\sum_{i=1}^n\sum_{l=1}^na_{il}x_ix_l-\sum_{i=1}^nb_ix_i,
```
y derivando respecto de ``x_j`` (los únicos términos que dependen de ``x_j`` son los que tienen ``i=j`` o ``l=j``),
```math
\frac{\partial\phi}{\partial x_j}(x) = \frac12\left(\sum_{l=1}^na_{jl}x_l+\sum_{i=1}^na_{ij}x_i\right)-b_j = \sum_{l=1}^na_{jl}x_l-b_j = (Ax-b)_j,
```
donde en la antepenúltima igualdad usamos que ``A`` es simétrica (``a_{ij}=a_{ji}``).

Alternativamente, y sin derivar: para ``x,z\in\mathbb{R}^n`` y ``t\in\mathbb{R}``, la cuenta hecha en la demostración del Teorema anterior (con ``x`` en lugar de ``x_*``) da
```math
\phi(x+tz) = \phi(x)+t\,z^T(Ax-b)+\frac{t^2}{2}z^TAz, \tag{5.1}
```
de donde ``\lim_{t\to0}\frac{\phi(x+tz)-\phi(x)}{t}=z^T(Ax-b)`` para todo ``z``, es decir, ``\nabla\phi(x)=Ax-b``.

Notemos que la simetría de ``A`` es esencial: si ``A`` no fuese simétrica se obtendría ``\nabla\phi(x)=\frac12(A+A^T)x-b``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000009
md"""
Del teorema anterior se sigue que si ``A`` es sdp,

Hallar ``x_*`` tal que ``Ax_*=b``

es equivalente a

Hallar ``x_*`` tal que ``x_*`` minimiza ``\phi`` en ``\mathbb{R}^n``.

Por otro lado, el siguiente lema será útil al momento de diseñar y entender los métodos para sistemas asociados a matrices sdp.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000000a
md"""
**Lema 5.3.** *Sean ``A`` sdp, ``b\in\mathbb{R}^n``, y ``x_*`` la solución de ``Ax=b``. Sea ``S\subset\mathbb{R}^n``. Entonces, ``y\in S`` minimiza ``\phi`` sobre ``S`` si y solo si ``y\in S`` minimiza ``\|\cdot-x_*\|_A`` sobre ``S``.*

**Demostración.** Se puede demostrar que (ejercicio)
```math
\phi(y) = \phi(x_*)+\frac12\|y-x_*\|_A^2, \qquad \forall y\in\mathbb{R}^n.
```
De esta última afirmación se sigue lo que afirma este lema.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000000b
md"""
## 5.2 Método del descenso más pronunciado (Steepest descent)

La idea es diseñar un método iterativo para hallar la solución ``x_*`` del sistema ``Ax=b``, "moviéndonos" en la dirección del máximo descenso de ``\phi``. Más precisamente, dado ``x_k\in\mathbb{R}^n``, nos "movemos" en la dirección del descenso más pronunciado ``r_k:=-\nabla\phi(x_k)=b-Ax_k`` para obtener ``x_{k+1}``; esto es, definimos
```math
x_{k+1} = x_k+\alpha r_k,
```
eligiendo ``\alpha\in\mathbb{R}`` de manera que ``\phi(x_k+\alpha r_k)`` sea lo más chico posible.

Si consideramos ``h(\alpha):=\phi(x_k+\alpha r_k)``, para ``\alpha\in\mathbb{R}``, tenemos que el mínimo de ``h`` se alcanza en ``\alpha=\frac{r_k^Tr_k}{r_k^TAr_k}``.

En efecto, usando (5.1) con ``x=x_k``, ``z=r_k`` y ``t=\alpha``, y recordando que ``Ax_k-b=-r_k``,
```math
h(\alpha) = \phi(x_k)+\alpha r_k^T(Ax_k-b)+\frac{\alpha^2}{2}r_k^TAr_k = \phi(x_k)-\alpha r_k^Tr_k+\frac{\alpha^2}{2}r_k^TAr_k.
```
Es decir, ``h`` es un polinomio cuadrático en ``\alpha`` con coeficiente principal ``\frac12r_k^TAr_k>0`` (pues ``A`` es sdp y ``r_k\neq0``); luego ``h`` tiene un único mínimo, que se alcanza donde ``h'(\alpha)=0``:
```math
h'(\alpha) = -r_k^Tr_k+\alpha r_k^TAr_k = 0 \iff \alpha = \frac{r_k^Tr_k}{r_k^TAr_k}.
```

Así, dado ``x_0\in\mathbb{R}^n``, la iteración del *método del descenso más pronunciado* está dada por
```math
x_{k+1} = x_k+\alpha r_k, \qquad k\geq0,
```
donde ``r_k=b-Ax_k`` y ``\alpha=\dfrac{r_k^Tr_k}{r_k^TAr_k}``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000000c
md"""
**Observación 5.4.** Si quisiéramos minimizar ``\phi(x_k+\alpha d_k)`` como función de ``\alpha``, donde ``d_k\in\mathbb{R}^n`` es alguna dirección dada, entonces dicho valor de ``\alpha`` vendría dado por ``\alpha=\dfrac{d_k^Tr_k}{d_k^TAd_k}``, donde ``r_k=b-Ax_k``.

**Algoritmo Steepest descent (teórico)**

```octave
Sean A sdp y b en R^n. Sea x0 en R^n.
r0 = b - A*x0;
k = 0;
mientras rk != 0          % ||rk|| > tol
    alpha = (rk' * rk) / (rk' * A * rk);
    x{k+1} = xk + alpha * rk;
    r{k+1} = b - A * x{k+1};
    k = k + 1;
fin
```

**Observación 5.5.** Recordemos que para la iteración de Richardson
```math
x_{k+1} = x_k-\alpha(Ax_k-b) = x_k+\alpha r_k,
```
donde ``r_k=b-Ax_k``. Así, el método del descenso más pronunciado es como la iteración de Richardson, pero eligiendo un valor de ``\alpha`` en cada paso.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000000d
md"""
**Observación 5.6.** En principio en cada iteración del método de Steepest descent se calculan dos productos mat-vec por iteración: ``Ar_k`` y ``Ax_{k+1}``.

Notemos que
```math
r_{k+1} = b-Ax_{k+1} = b-A(x_k+\alpha r_k) = b-Ax_k-\alpha Ar_k = r_k-\alpha Ar_k.
```
Reemplazando en el algoritmo anterior la línea ``r_{k+1}=b-Ax_{k+1}``; por ``r_{k+1}=r_k-\alpha Ar_k``; podemos reducir el costo a un solo producto mat-vec por iteración.

**Observación 5.7.** Notemos que en el algoritmo anterior no hace falta conocer la matriz ``A``, es suficiente disponer de una rutina matriz-vector, esto es, saber cómo actúa ``A`` sobre cada vector ``x``, o sea, conocer la acción ``Ax``, para cada ``x``. Esto puede resultar mucho más eficiente en algunos casos.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000000e
md"""
**Ejemplo 5.8.** Si ``A`` es de ``2\times2``, las curvas de nivel de ``\phi`` son elipses. Cuanto mayor es el cociente ``\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}`` más alargadas son. Pensemos por ejemplo en ``A=\begin{pmatrix}\lambda_{\max}(A)&0\\0&\lambda_{\min}(A)\end{pmatrix}``. Cuanto más alargadas, más lenta será la convergencia.

Expliquemos esto con más detalle, porque es la clave para entender la lentitud del método.

**Las curvas de nivel son elipses centradas en ``x_*``.** Por la identidad usada en la demostración del Lema 5.3,
```math
\phi(y) = \phi(x_*)+\frac12\|y-x_*\|_A^2,
```
de modo que el conjunto de nivel ``\{y:\phi(y)=c\}`` es ``\{y:\|y-x_*\|_A^2=2(c-\phi(x_*))\}``: un elipsoide centrado en ``x_*`` (y vacío si ``c<\phi(x_*)``). Si ``A=V\Lambda V^T`` con ``V`` ortogonal y ``\Lambda=\text{diag}(\lambda_1,\dots,\lambda_n)``, escribiendo ``w=V^T(y-x_*)`` el conjunto de nivel se lee
```math
\sum_{i=1}^n\lambda_iw_i^2 = 2(c-\phi(x_*)),
```
es decir, un elipsoide con ejes en las direcciones de los autovectores de ``A`` y semieje ``i``-ésimo proporcional a ``1/\sqrt{\lambda_i}``. Por lo tanto,
```math
\frac{\text{semieje mayor}}{\text{semieje menor}} = \sqrt{\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}} = \sqrt{\kappa_2(A)}.
```
Cuando ``A`` es un múltiplo de la identidad las curvas de nivel son esferas; cuanto mayor es ``\kappa_2(A)``, más alargadas son.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000000f
md"""
**Por qué eso hace lento al método.** El residuo ``r_k=-\nabla\phi(x_k)`` es ortogonal (en el sentido euclídeo) a la curva de nivel que pasa por ``x_k``, pero *no* apunta hacia el centro ``x_*``, salvo que la curva de nivel sea una esfera o que ``x_k-x_*`` sea un autovector de ``A``. Cuanto más alargadas son las elipses, mayor puede ser ese desvío.

Además, la elección óptima de ``\alpha`` hace que dos residuos consecutivos sean ortogonales:
```math
r_{k+1}^Tr_k = 0, \qquad k\geq0. \tag{5.2}
```
En efecto, si ``h(\alpha)=\phi(x_k+\alpha r_k)``, la regla de la cadena y la Observación 5.2 dan
```math
h'(\alpha) = \nabla\phi(x_k+\alpha r_k)^Tr_k = \big(A(x_k+\alpha r_k)-b\big)^Tr_k = -r(\alpha)^Tr_k,
```
donde ``r(\alpha):=b-A(x_k+\alpha r_k)``; evaluando en el ``\alpha`` que anula ``h'``, que es el que define ``x_{k+1}``, resulta ``r(\alpha)=r_{k+1}`` y se obtiene (5.2). Geométricamente: el método avanza en línea recta hasta ser tangente a una curva de nivel, y ahí dobla en ángulo recto. Sobre elipses muy alargadas esto produce el característico *zigzag* de pasos cortos, en el que cada iteración deshace parte del avance de la anterior (Figura 5.1).
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000010
md"""
**Figura 5.1** (descripción, no reproducida gráficamente). Descenso más pronunciado para ``A=\text{diag}(\lambda_{\min},\lambda_{\max})`` con ``\kappa_2(A)=9`` y ``b=0`` (o sea ``x_*=0``), partiendo de ``x_0-x_*=(1/\lambda_{\min},1/\lambda_{\max})``. Los residuos consecutivos son ortogonales y el método zigzaguea; en cada paso el error en norma ``A`` se reduce exactamente por el factor ``\frac{\kappa_2(A)-1}{\kappa_2(A)+1}=0.8``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000011
md"""
**El peor caso se alcanza.** El ejemplo de la Figura 5.1 muestra que la cota de convergencia que probaremos a continuación no se puede mejorar. Tomemos ``n=2``, ``A=\text{diag}(\lambda_1,\lambda_2)`` con ``0<\lambda_1<\lambda_2``, y ``b=0``, de modo que ``x_*=0`` y ``r_k=-Ax_k``. Si ``x_k=\left(\frac{1}{\lambda_1},\frac{1}{\lambda_2}\right)``, entonces ``r_k=-(1,1)`` y
```math
\alpha = \frac{r_k^Tr_k}{r_k^TAr_k} = \frac{2}{\lambda_1+\lambda_2},
```
con lo cual
```math
x_{k+1} = \left(\frac1{\lambda_1}(1-\alpha\lambda_1),\ \frac1{\lambda_2}(1-\alpha\lambda_2)\right) = \frac{\lambda_2-\lambda_1}{\lambda_2+\lambda_1}\left(\frac1{\lambda_1},-\frac1{\lambda_2}\right).
```
El nuevo error tiene la misma forma que el anterior (salvo un cambio de signo en la segunda coordenada, que se repite alternadamente), y por lo tanto
```math
\|x_{k+1}-x_*\|_A = \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A \qquad\text{para todo } k\geq0,
```
o sea, la desigualdad del teorema siguiente vale con *igualdad* en cada paso. Notemos que, en cambio, si ``x_0-x_*`` es un autovector de ``A``, el método converge en una sola iteración: el problema no es la matriz mal condicionada *en sí*, sino la combinación de un espectro disperso con un error inicial que mezcla las direcciones extremas.

Finalmente, establecemos un resultado de convergencia para el método que estamos estudiando.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000012
md"""
**Teorema 5.9.** *Sea ``A`` sdp, ``b\in\mathbb{R}^n``, y ``x_*`` la solución de ``Ax=b``. Sea ``x_0\in\mathbb{R}^n`` y sea ``\{x_k\}_{k\in\mathbb{N}}`` la sucesión generada por el método del descenso más pronunciado. Entonces,*
```math
\|x_{k+1}-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A, \qquad k\geq0. \tag{5.3}
```

Notemos que es la misma estimación de error de Richardson cuando se elige el ``\alpha`` óptimo (``\alpha_{\text{opt}}``). ¿Cuál le parece que será entonces la "ventaja" de este método respecto del anterior?

**Demostración.** Sea ``k\geq0`` arbitrario y fijo, y supongamos que ``x_k`` está dado. La iteración de Richardson con ``\alpha_{\text{opt}}`` da
```math
x_{k+1}^R = x_k-\alpha_{\text{opt}}(Ax_k-b) = x_k+\alpha_{\text{opt}}r_k,
```
con ``r_k=b-Ax_k``. Por otro lado, para el descenso más pronunciado,
```math
x_{k+1} = x_k+\alpha r_k,
```
donde ``\alpha`` se elige de modo de minimizar ``h(\alpha):=\phi(x_k+\alpha r_k)``. Por el Lema 5.3, esto es equivalente a minimizar ``\|x_k+\alpha r_k-x_*\|_A``. Así,
```math
\|\underbrace{x_k+\alpha r_k}_{x_{k+1}}-x_*\|_A \leq \|\underbrace{x_k+\alpha_{\text{opt}}r_k}_{x_{k+1}^R}-x_*\|_A.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000013
md"""
Del teorema de convergencia para la iteración de Richardson se tiene que
```math
\|x_{k+1}^R-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A,
```
y por lo tanto,
```math
\|x_{k+1}-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A,
```
como queríamos demostrar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000014
md"""
## 5.3 Método del gradiente conjugado

Seguimos considerando la resolución de un sistema ``Ax=b``, donde ``A`` es sdp.

Pregunta: ¿Se puede "mejorar" el método del descenso más pronunciado eligiendo "mejor" las direcciones de búsqueda?

El método del gradiente conjugado que estudiaremos en esta sección surgió en 1950 como un método directo. Posteriormente se empezó a usar como método iterativo para grandes sistemas de ecuaciones lineales y ha sustituido el uso de los métodos iterativos estacionarios como Jacobi, Gauss-Seidel, etc.

Recordemos el algoritmo del método del descenso más pronunciado:
```octave
Sea x0 en R^n.
r0 = b - A*x0;
k = 0;
mientras rk != 0
    x{k+1} = xk + alpha_k * rk, donde elegimos adecuadamente el valor alpha_k;
    r{k+1} = b - A * x{k+1};
    k = k + 1;
fin
```

Recordemos que ``\alpha_k`` se elige de manera que
```math
x_{k+1} \text{ minimiza } \phi \text{ sobre } \{x_k+\alpha r_k \mid \alpha\in\mathbb{R}\}.
```
O sea,
```math
\begin{aligned}
x_1 &\text{ minimiza } \phi \text{ sobre } \{x_0+\alpha r_0 \mid \alpha\in\mathbb{R}\}\\
x_2 &\text{ minimiza } \phi \text{ sobre } \{x_0+\alpha_0 r_0+\alpha r_1 \mid \alpha\in\mathbb{R}\}\\
x_3 &\text{ minimiza } \phi \text{ sobre } \{x_0+\alpha_0 r_0+\alpha_1 r_1+\alpha r_2 \mid \alpha\in\mathbb{R}\}\\
&\ \ \vdots
\end{aligned}
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000015
md"""
Dado ``x_0\in\mathbb{R}^n``, definimos una nueva sucesión ``\{x_k\}_{k\in\mathbb{N}}`` de la siguiente manera:
```math
x_{k+1} \text{ minimiza } \phi \text{ sobre } \underbrace{x_0+\text{span}\{r_0,r_1,\dots,r_k\}}_{:=\{x_0+\alpha_0r_0+\alpha_1r_1+\cdots+\alpha_kr_k\,\mid\,\alpha_0,\dots,\alpha_k\in\mathbb{R}\}}, \qquad k\geq0,
```
donde ``r_i=b-Ax_i``, para ``i=0,1,\dots,k``.

Recordando que minimizar ``\phi`` es equivalente a minimizar la ``A``-distancia a la solución ``x_*``, se tiene que

``x_{k+1}`` es el punto de ``x_0+\text{span}\{r_0,r_1,\dots,r_k\}`` más cercano (en norma ``A``) a la solución ``x_*``.

Verificar que tal punto siempre existe, ¡para asegurarse de que ``x_{k+1}`` está bien definido!

**Método del gradiente conjugado.** (Versión preliminar 1)

```octave
Sean A sdp y b en R^n. Sea x0 en R^n.
r0 = b - A*x0;
k = 0;
mientras rk != 0
    x{k+1} = argmin_{x en x0+span{r0,r1,...,rk}} phi(x);
    r{k+1} = b - A * x{k+1};
    k = k + 1;
fin
```

A continuación, nos abocaremos al estudio de ciertas propiedades interesantes del método del gradiente conjugado y de su convergencia. Posteriormente, analizaremos si hay alguna manera computacional relativamente económica para calcular efectivamente los ``x_k``'s.

En el siguiente lema demostramos que los residuos ``r_k`` son ortogonales.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000016
md"""
**Lema 5.10.** *Sea ``\{r_k\}_{k\in\mathbb{N}_0}`` la sucesión de residuos del método del gradiente conjugado. Entonces,*
```math
r_k\perp\{r_0,r_1,\dots,r_{k-1}\}, \qquad k\geq1.
```

**Demostración.** Sea ``k\geq1`` arbitrario y fijo. Puesto que ``x_k=\text{argmin}_{x\in x_0+\text{span}\{r_0,\dots,r_{k-1}\}}\phi(x)``, se tiene que
```math
\phi(x_k)\leq\phi(x), \qquad \forall x\in x_0+\text{span}\{r_0,\dots,r_{k-1}\}.
```
Sea ``\xi\in\text{span}\{r_0,\dots,r_{k-1}\}`` arbitrario. Entonces ``\phi(x_k)\leq\phi(x_k+t\xi)``, para todo ``t\in\mathbb{R}``. Así, la función ``f_\xi(t):=\phi(x_k+t\xi)`` alcanza su mínimo en ``t=0``, y por lo tanto, ``f_\xi'(0)=0``.

Puesto que
```math
f_\xi'(t) = \xi^T\nabla\phi(x_k+t\xi), \qquad \forall t\in\mathbb{R},
```
se tiene que
```math
f_\xi'(0) = \xi^T\nabla\phi(x_k) = \xi^T\underbrace{(Ax_k-b)}_{=-r_k} = 0,
```
y así, ``r_k\perp\xi``. Puesto que ``\xi\in\text{span}\{r_0,\dots,r_{k-1}\}`` es arbitrario,
```math
r_k \perp \{r_0,r_1,\dots,r_{k-1}\}.
```

**Observación 5.11.** Si ``r_k\neq0`` entonces ``\{r_0,r_1,\dots,r_k\}`` es linealmente independiente.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000017
md"""
**Lema 5.12.** *Si ``r_k\neq0`` entonces*
```math
\text{span}\{r_0,r_1,\dots,r_k\} = \text{span}\{r_0,Ar_0,\dots,A^kr_0\}, \qquad \text{para todo } k\geq0.
```
El espacio ``\mathcal{K}_{k+1}:=\text{span}\{r_0,Ar_0,\dots,A^kr_0\}`` se llama *espacio de Krylov*.

**Demostración.** Demostraremos el resultado haciendo inducción sobre ``k``.

Notemos que si ``k=0``, ``\text{span}\{r_0\}=\mathcal{K}_1``. Supongamos ahora que
```math
\text{span}\{r_0,\dots,r_k\} = \mathcal{K}_{k+1}, \qquad\text{(Hip. ind.)}
```
y veamos que ``\text{span}\{r_0,\dots,r_k,r_{k+1}\}=\mathcal{K}_{k+2}``, siempre que ``r_{k+1}\neq0``.

Demostramos primero que ``\text{span}\{r_0,\dots,r_k,r_{k+1}\}\subset\mathcal{K}_{k+2}``. Para esto, notemos que por hipótesis inductiva,
```math
r_i\in\mathcal{K}_{k+1}\subset\mathcal{K}_{k+2}, \qquad\text{para } i=0,\dots,k.
```
Solo resta ver entonces que ``r_{k+1}\in\mathcal{K}_{k+2}``. Puesto que ``r_{k+1}=b-Ax_{k+1}``, donde ``x_{k+1}=x_0+\xi``, con ``\xi\in\text{span}\{r_0,\dots,r_k\}\underset{\text{por hip. ind.}}{=}\mathcal{K}_{k+1}``, se tiene que
```math
r_{k+1} = b-Ax_{k+1} = b-Ax_0-A\xi = \underbrace{r_0}_{\in\mathcal{K}_{k+2}}-\underbrace{A\xi}_{\in\mathcal{K}_{k+2},\ \text{pues}\ \xi\in\mathcal{K}_{k+1}} \in\mathcal{K}_{k+2}.
```

Por lo tanto, ``\text{span}\{r_0,\dots,r_k,r_{k+1}\}\subset\mathcal{K}_{k+2}``.

Ahora bien, puesto que ``r_{k+1}\neq0``, tenemos que ``\text{span}\{r_0,\dots,r_k,r_{k+1}\}`` es un subespacio de dimensión ``k+2`` del espacio de Krylov ``\mathcal{K}_{k+2}`` que tiene dimensión a lo sumo ``k+2``, y así, necesariamente, ``\text{span}\{r_0,\dots,r_k,r_{k+1}\}=\mathcal{K}_{k+2}``, como queríamos demostrar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000018
md"""
## 5.4 Cotas de error para el método del gradiente conjugado

Recordemos que de la definición de ``x_k`` se sigue que ``x_k`` minimiza ``\phi`` en ``x_0+\mathcal{K}_k``, y esto ocurre, si y solo si
```math
\|x_k-x_*\|_A \leq \|w-x_*\|_A, \qquad \forall w\in x_0+\mathcal{K}_k.
```
Si ``w\in x_0+\mathcal{K}_k``, entonces ``w=x_0+\sum_{i=0}^{k-1}\gamma_iA^ir_0``, con ``\gamma_0,\dots,\gamma_{k-1}\in\mathbb{R}``. Entonces,
```math
\begin{align*}
w-x_* &= x_0-x_*+\sum_{i=0}^{k-1}\gamma_iA^ir_0\\
&= x_0-x_*+\sum_{i=0}^{k-1}\gamma_iA^i(b-Ax_0)\\
&= x_0-x_*+\sum_{i=0}^{k-1}\gamma_iA^i(Ax_*-Ax_0)\\
&= x_0-x_*+\sum_{i=0}^{k-1}\gamma_iA^{i+1}(x_*-x_0)\\
&= \left(I-\sum_{i=0}^{k-1}\gamma_iA^{i+1}\right)(x_0-x_*)\\
&= p_k(A)(x_0-x_*),
\end{align*}
```
donde ``p_k(z)=1-\sum_{i=0}^{k-1}\gamma_iz^{i+1}``. Notemos que ``p_k\in\mathcal{P}_k`` y que ``p_k(0)=1``.

Definimos el conjunto ``\overline{\mathcal{P}}_k`` de polinomios residuales de grado ``\leq k`` de la siguiente manera:
```math
\overline{\mathcal{P}}_k := \{p\in\mathcal{P}_k \mid p(0)=1\}.
```
Así, hemos demostrado el siguiente teorema.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000019
md"""
**Teorema 5.13.** *Sean ``A`` sdp y ``b\in\mathbb{R}^n``. Sea ``x_*`` la solución del sistema ``Ax=b``. Si ``\{x_k\}_{k\in\mathbb{N}}`` denota la sucesión generada por el método del gradiente conjugado partiendo de un vector inicial ``x_0\in\mathbb{R}^n``, entonces,*
```math
\|x_k-x_*\|_A = \min_{p\in\overline{\mathcal{P}}_k}\|p(A)(x_0-x_*)\|_A, \qquad k\geq0.
```

**Corolario 5.14.** *Sean ``A`` sdp y ``b\in\mathbb{R}^n``. Sea ``x_*`` la solución del sistema ``Ax=b``. Si ``\{x_k\}_{k\in\mathbb{N}}`` denota la sucesión generada por el método del gradiente conjugado partiendo de un vector inicial ``x_0\in\mathbb{R}^n``, entonces,*
```math
\|x_k-x_*\|_A \leq \max_{\lambda\in\sigma(A)}|p(\lambda)|\,\|x_0-x_*\|_A, \qquad \forall p\in\overline{\mathcal{P}}_k.
```

**Demostración.** Del teorema anterior se sigue que
```math
\|x_k-x_*\|_A \leq \|p(A)\|_A\|x_0-x_*\|_A,
```
para todo ``p\in\overline{\mathcal{P}}_k``. Por lo tanto es suficiente demostrar que
```math
\|p(A)\|_A \leq \max_{\lambda\in\sigma(A)}|p(\lambda)|, \qquad \forall p\in\overline{\mathcal{P}}_k. \tag{5.4}
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000001a
md"""
Sea ``p\in\overline{\mathcal{P}}_k``. Entonces, para todo ``x\in\mathbb{R}^n`` se tiene que
```math
\begin{align*}
\|p(A)x\|_A^2 &= (p(A)x)^TA(p(A)x) = (A^{1/2}p(A)x)^T(A^{1/2}p(A)x) = \|A^{1/2}p(A)x\|_2^2\\
&= \|p(A)A^{1/2}x\|_2^2 \leq \|p(A)\|_2^2\|A^{1/2}x\|_2^2.
\end{align*}
```
Puesto que ``\|A^{1/2}x\|_2^2=(A^{1/2}x)^T(A^{1/2}x)=x^TA^{1/2}A^{1/2}x=x^TAx=\|x\|_A^2``, tenemos que
```math
\|p(A)x\|_A \leq \|p(A)\|_2\|x\|_A, \qquad\text{para todo } x\in\mathbb{R}^n,
```
y así,
```math
\|p(A)\|_A \leq \|p(A)\|_2.
```
En consecuencia, para demostrar (5.4) es suficiente demostrar que
```math
\|p(A)\|_2 \leq \max_{\lambda\in\sigma(A)}|p(\lambda)|, \qquad \forall p\in\overline{\mathcal{P}}_k. \tag{5.5}
```

Puesto que ``A`` es simétrica, tenemos que ``A`` es diagonalizable y que existe una base ortonormal de vectores propios. Si ``\lambda_1,\dots,\lambda_n`` son los autovalores de ``A``, entonces,
```math
A = V\Lambda V^T, \qquad\text{donde } \Lambda = \begin{pmatrix}\lambda_1&0&\dots&0\\0&\lambda_2&\dots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\dots&\lambda_n\end{pmatrix},
```
y ``V`` es una matriz ortogonal (``VV^T=I``) cuyas columnas son correspondientes vectores propios.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000001b
md"""
Sea ``p\in\overline{\mathcal{P}}_k``. Entonces, ``p(z)=1-\sum_{i=0}^{k-1}\gamma_iz^{i+1}`` para algunos ``\gamma_0,\gamma_1,\dots,\gamma_{k-1}\in\mathbb{R}``. Así, puesto que ``VV^T=I``,
```math
\begin{align*}
p(A) &= I-\sum_{i=0}^{k-1}\gamma_iA^{i+1} = I-\sum_{i=0}^{k-1}\gamma_i(V\Lambda V^T)^{i+1}\\
&= I-\sum_{i=0}^{k-1}\gamma_iV\Lambda^{i+1}V^T = V\left(I-\sum_{i=0}^{k-1}\gamma_i\Lambda^{i+1}\right)V^T = Vp(\Lambda)V^T.
\end{align*}
```
Es fácil verificar que
```math
p(\Lambda) = \begin{pmatrix}p(\lambda_1)&0&\dots&0\\0&p(\lambda_2)&\dots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\dots&p(\lambda_n)\end{pmatrix},
```
y en consecuencia, ``p(A)`` es diagonalizable y
```math
\sigma(p(A)) = \{p(\lambda) \mid \lambda\in\sigma(A)\}.
```
Finalmente, puesto que ``p(A)`` es simétrica, ``\|p(A)\|_2=\rho(p(A))`` y así, (5.5) queda demostrada.

La utilidad del resultado anterior radica en que si se dispone de información sobre el espectro de ``A``, se pueden construir un polinomio residual conveniente y obtener una cota para el error de las iteraciones del método del gradiente conjugado.

Veamos algunos ejemplos de cómo hacer esto. En todos ellos la estrategia es la misma: elegir un ``p\in\overline{\mathcal{P}}_k`` que sea *chico sobre el espectro de* ``A``, sin preocuparse por lo que ``p`` hace fuera de él, y aplicar el Corolario 5.14.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000001c
md"""
**Ejemplo 5.15** (Espectro agrupado alrededor de un punto). Supongamos que ``\sigma(A)\subset[1-\varepsilon,1+\varepsilon]`` con ``0<\varepsilon<1``. Tomando
```math
p(z) = (1-z)^k \in \overline{\mathcal{P}}_k,
```
se tiene ``|p(\lambda)|=|1-\lambda|^k\leq\varepsilon^k`` para todo ``\lambda\in\sigma(A)``, y por lo tanto
```math
\|x_k-x_*\|_A \leq \varepsilon^k\|x_0-x_*\|_A.
```
Si ``\varepsilon`` es chico la convergencia es muy rápida, *independientemente de* ``n``. Este es exactamente el efecto que busca un buen precondicionador (Sección 5.1.6): no achicar los autovalores, sino *agruparlos*.

**Ejemplo 5.16** (Pocos autovalores distintos). Si ``A`` tiene solamente ``m`` autovalores distintos ``\mu_1,\dots,\mu_m``, tomamos
```math
p(z) = \prod_{j=1}^m\left(1-\frac{z}{\mu_j}\right) \in \overline{\mathcal{P}}_m.
```
Entonces ``p(\lambda)=0`` para todo ``\lambda\in\sigma(A)`` y el corolario da ``\|x_m-x_*\|_A\leq0``, es decir, ``x_m=x_*``: el método termina en a lo sumo ``m`` pasos. Volveremos sobre esto en la Sección 5.1.4.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000001d
md"""
**Ejemplo 5.17** (Autovalores aislados: el espectro "casi" agrupado). El caso realmente interesante es el intermedio: la mayoría de los autovalores agrupados y unos pocos lejos. Supongamos que
```math
\sigma(A) \subset [a,b]\cup\{\mu_1,\dots,\mu_m\},
```
con ``0<a\leq b`` y ``m`` autovalores aislados ``\mu_j>b``. La idea es *gastar* ``m`` factores del polinomio en anular los autovalores aislados, y usar los ``k-m`` grados de libertad restantes sobre el intervalo. Para ``k\geq m`` tomamos
```math
p(z) = \underbrace{\prod_{j=1}^m\left(1-\frac{z}{\mu_j}\right)}_{=:q(z)} \cdot p_{k-m}(z),
```
donde ``p_{k-m}\in\overline{\mathcal{P}}_{k-m}`` es el polinomio de Tchebychev normalizado sobre ``[a,b]`` que aparece en la demostración del Teorema siguiente. Como ``p(0)=q(0)p_{k-m}(0)=1``, resulta ``p\in\overline{\mathcal{P}}_k``. Además ``p`` se anula en cada ``\mu_j``, y sobre ``[a,b]``,
```math
\max_{\lambda\in[a,b]}|p(\lambda)| \leq C\max_{\lambda\in[a,b]}|p_{k-m}(\lambda)| \leq 2C\left(\frac{\sqrt{b/a}-1}{\sqrt{b/a}+1}\right)^{k-m}, \qquad C:=\max_{\lambda\in[a,b]}|q(\lambda)|.
```
En consecuencia,
```math
\|x_k-x_*\|_A \leq 2C\left(\frac{\sqrt{b/a}-1}{\sqrt{b/a}+1}\right)^{k-m}\|x_0-x_*\|_A, \qquad k\geq m.
```
Lo notable es que, salvo las primeras ``m`` iteraciones y la constante ``C``, la velocidad de convergencia está gobernada por ``b/a`` y *no* por ``\kappa_2(A)=\mu_{\max}/a``. Por ejemplo, si ``\sigma(A)\subset[1,2]\cup\{1000\}``, entonces ``\kappa_2(A)=1000`` y la cota general del Teorema siguiente da un factor
```math
\frac{\sqrt{1000}-1}{\sqrt{1000}+1} \approx 0.939
```
por iteración, mientras que con ``m=1``, ``a=1``, ``b=2`` y ``C=\max_{\lambda\in[1,2]}|1-\lambda/1000|\leq1`` obtenemos
```math
\|x_k-x_*\|_A \leq 2\left(\frac{\sqrt2-1}{\sqrt2+1}\right)^{k-1}\|x_0-x_*\|_A \approx 2(0.172)^{k-1}\|x_0-x_*\|_A.
```
Para reducir el error en un factor ``10^{-6}``, la primera cota pide unas 230 iteraciones y la segunda, unas 9. Esto explica una observación experimental frecuente: *el gradiente conjugado suele andar mucho mejor que lo que predice la cota basada solo en* ``\kappa_2(A)``, *sobre todo cuando el espectro está agrupado con pocos valores extremos*.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000001e
md"""
**Observación 5.18.** Vale la pena insistir en la lógica de estos argumentos: el Corolario 5.14 vale para *todo* ``p\in\overline{\mathcal{P}}_k``, y el gradiente conjugado no "sabe" cuál elegimos ni lo calcula. Nosotros elegimos un ``p`` conveniente solo para *estimar* el error; el método siempre hace algo al menos tan bueno como el mejor de esos ``p``. Por eso todas las cotas obtenidas de esta manera son cotas superiores, y cualquier información adicional sobre ``\sigma(A)`` puede dar una cota mejor.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000001f
md"""
**Teorema 5.19.** *Sea ``A`` sdp, ``b\in\mathbb{R}^n``, y ``x_*`` la solución de ``Ax=b``. Sea ``x_0\in\mathbb{R}^n`` y sea ``\{x_k\}_{k\in\mathbb{N}}`` la sucesión generada por el método del gradiente conjugado. Entonces,*
```math
\|x_k-x_*\|_A \leq 2\left(\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\right)^k\|x_0-x_*\|_A, \qquad k\geq0. \tag{5.6}
```

**Demostración.** Sean ``0<\lambda_1\leq\lambda_2\leq\dots\leq\lambda_n`` los autovalores de ``A``, y sean ``\lambda_{\min}=\lambda_1`` y ``\lambda_{\max}=\lambda_n`` el menor y el mayor de ellos.

Recordemos los polinomios de Tchebychev definidos por:
```math
\begin{cases}
T_0(x)=1,\\
T_1(x)=x,\\
T_k(x)=2xT_{k-1}(x)-T_{k-2}(x), \qquad k=2,3,\dots
\end{cases}
```
Sabemos que ``T_k\in\mathcal{P}_k`` y que
```math
T_k(x) = \cos(k\cos^{-1}(x)), \qquad x\in[-1,1],
```
donde la función ``\cos:[0,\pi]\to[-1,1]``. Teniendo en cuenta la biyección ``x:[\lambda_{\min},\lambda_{\max}]\to[-1,1]`` dada por ``x(\lambda)=\dfrac{\lambda-\frac{\lambda_{\max}+\lambda_{\min}}{2}}{\frac{\lambda_{\max}-\lambda_{\min}}{2}}``, definimos
```math
p_k(\lambda) := \frac{T_k(x(\lambda))}{T_k(x(0))} \in \overline{\mathcal{P}}_k.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000020
md"""
Para aplicar el Corolario 5.14, tenemos que acotar ``\max_{\lambda\in\sigma(A)}|p_k(\lambda)|``. Puesto que ``\max_{x\in[-1,1]}|T_k(x)|\leq1``, se tiene que
```math
\max_{\lambda\in\sigma(A)}|p_k(\lambda)| \leq \max_{\lambda_{\min}\leq\lambda\leq\lambda_{\max}}|p_k(\lambda)| = \max_{\lambda_{\min}\leq\lambda\leq\lambda_{\max}}\frac{|T_k(x(\lambda))|}{|T_k(x(0))|} \leq \frac{1}{|T_k(x(0))|}, \tag{5.7}
```
y en consecuencia, es suficiente acotar por debajo ``|T_k(x(0))|``. Por un lado, llamando ``c=\kappa_2(A)=\dfrac{\lambda_{\max}}{\lambda_{\min}}``, tenemos
```math
x(0) = -\frac{\lambda_{\max}+\lambda_{\min}}{\lambda_{\max}-\lambda_{\min}} = -\frac{c+1}{c-1} = \frac12\left(\underbrace{\left(-\frac{\sqrt c-1}{\sqrt c+1}\right)}_{=:z}+\underbrace{\left(-\frac{\sqrt c+1}{\sqrt c-1}\right)}_{=\frac1z}\right).
```
Por otro lado, es fácil ver que (ejercicio!), si ``z\neq0``,
```math
x(0) = \frac{z+z^{-1}}{2} \implies T_k(x(0)) = \frac{z^k+z^{-k}}{2}.
```
Así,
```math
T_k(x(0)) = \frac12\left(\left(-\frac{\sqrt c-1}{\sqrt c+1}\right)^k+\left(-\frac{\sqrt c+1}{\sqrt c-1}\right)^k\right) = (-1)^k\frac12\left(\left(\frac{\sqrt c-1}{\sqrt c+1}\right)^k+\left(\frac{\sqrt c+1}{\sqrt c-1}\right)^k\right),
```
y luego,
```math
|T_k(x(0))| \geq \frac12\left(\frac{\sqrt c+1}{\sqrt c-1}\right)^k. \tag{5.8}
```
De (5.7) y (5.8), teniendo en cuenta el Corolario 5.14, concluye la demostración.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000021
md"""
**Observación 5.20** (Cuántas iteraciones hacen falta). De la cota (5.6) se deduce cuántas iteraciones alcanzan para reducir el error en un factor dado, y es de aquí que sale la afirmación de que el gradiente conjugado necesita ``O(\sqrt{\kappa_2(A)})`` iteraciones.

Fijemos una tolerancia ``\varepsilon>0`` y preguntémonos cuántas iteraciones garantizan que
```math
\|x_k-x_*\|_A \leq \varepsilon\|x_0-x_*\|_A.
```
Escribiendo ``\kappa=\kappa_2(A)`` y ``r=\frac{\sqrt\kappa-1}{\sqrt\kappa+1}``, por (5.6) basta con que ``2r^k\leq\varepsilon``, es decir, con que
```math
k\ln\frac1r \geq \ln\frac2\varepsilon \iff k \geq \frac{\ln(2/\varepsilon)}{\ln(1/r)}. \tag{5.9}
```
Toda la dependencia en ``\kappa`` está en ``\ln\frac1r``, que podemos pensar como la *tasa de convergencia por iteración*. Para estimarla, escribimos ``r`` como "uno menos algo chico":
```math
r = \frac{\sqrt\kappa-1}{\sqrt\kappa+1} = 1-t, \qquad t:=\frac{2}{\sqrt\kappa+1}\in(0,1),
```
y usamos que ``\ln(1-t)\leq-t``, o sea ``\ln\frac1{1-t}\geq t``, de donde
```math
\ln\frac1r \geq \frac{2}{\sqrt\kappa+1}.
```
Reemplazando en (5.9), vemos que alcanza con tomar
```math
k \geq \frac{\sqrt{\kappa_2(A)}+1}{2}\ln\frac2\varepsilon, \tag{5.10}
```
ya que ese ``k`` verifica ``k\ln\frac1r\geq\ln\frac2\varepsilon``. Con ``\varepsilon`` fijo, el segundo factor es una constante, y por lo tanto
```math
k = O\big(\sqrt{\kappa_2(A)}\big).
```

La estimación de ``\ln\frac1r`` no es floja: como ``\ln\frac1{1-t}=t+\frac{t^2}2+\cdots``, para ``\kappa`` grande vale ``\ln\frac1r=\frac2{\sqrt\kappa}+O(\kappa^{-3/2})``, de modo que el orden ``\sqrt\kappa`` es el correcto.

El mismo razonamiento aplicado al método del descenso más pronunciado, cuyo factor por iteración es ``\frac{\kappa-1}{\kappa+1}=1-\frac2{\kappa+1}``, da ``k=O(\kappa_2(A))``: reaparece aquí la diferencia entre ``\kappa`` y ``\sqrt\kappa`` señalada más arriba.

Dos advertencias. Primero, (5.10) es una cota *superior*: dice que esa cantidad de iteraciones *alcanza*, no que sea necesaria; como ya observamos, el gradiente conjugado suele comportarse bastante mejor que lo que predice (5.6), sobre todo si el espectro está agrupado. Segundo, la dependencia en la tolerancia es apenas logarítmica: pedir una cifra más de precisión cuesta un número *fijo* de iteraciones adicionales, mientras que duplicar ``\sqrt\kappa`` duplica el total. Por eso tiene sentido hablar de "las iteraciones necesarias para una tolerancia fija" sin especificar cuál es.
"""

# ╔═╡ e5f6a7b8-dddd-4a1b-8c2d-000000000001
md"""
**Nota** (de las diapositivas de la Clase 9). Conviene tener a la vista, numéricamente, la diferencia entre la tasa del descenso más pronunciado (equivalentemente, la de Richardson con parámetro óptimo) y la del gradiente conjugado, ambas en función de ``\kappa=\kappa_2(A)``:

| Método | Cota del factor de reducción por iteración | ``\kappa=100`` |
|---|---|---|
| Richardson óptimo / SD | ``\left(\dfrac{\kappa-1}{\kappa+1}\right)^k`` | ``(0.980)^k`` |
| GC (Tchebychev) | ``2\left(\dfrac{\sqrt\kappa-1}{\sqrt\kappa+1}\right)^k`` | ``2\,(0.818)^k`` |

Para reducir el error en un factor ``10^{-2}`` con ``\kappa=100``, el descenso más pronunciado necesita del orden de 230 iteraciones, mientras que el gradiente conjugado necesita del orden de 30: la diferencia entre ``\kappa`` y ``\sqrt\kappa`` en el exponente crece a medida que empeora el condicionamiento, y la constante 2 de la cota del GC resulta asintóticamente irrelevante frente a esa diferencia. Recordar también que esta es solo una cota (usa un polinomio ``p`` particular); en la práctica el gradiente conjugado suele comportarse todavía mejor, sobre todo si el espectro está agrupado.

Vale además la siguiente observación general: en cada paso, el gradiente conjugado es al menos tan bueno como el descenso más pronunciado, porque el afín sobre el que minimiza SD (la recta ``x_k+\alpha r_k``) está contenido en el afín sobre el que minimiza GC (``x_0+\mathcal{K}_{k+1}``).
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000022
md"""
### 5.4.1 El método del gradiente conjugado como método directo

Es claro que el método del gradiente conjugado converge en a lo sumo ``n`` pasos, ya que si ``r_n=b-Ax_n`` fuera no nulo, el conjunto ``\{r_0,r_1,\dots,r_n\}`` sería linealmente independiente en ``\mathbb{R}^n``, lo cual no es posible. Más aún, tenemos el siguiente resultado.

**Teorema 5.21.** *Sea ``A`` sdp y ``b\in\mathbb{R}^n``. Si ``A`` tiene ``k`` autovalores distintos entonces el método del gradiente conjugado para resolver el sistema ``Ax=b`` converge a la solución en a lo sumo ``k`` pasos, partiendo de cualquier valor inicial.*

**Demostración.** Sean ``\lambda_1,\lambda_2,\dots,\lambda_k`` los ``k\leq n`` autovalores distintos de ``A``. Definamos
```math
p_k(z) := \frac{(\lambda_1-z)(\lambda_2-z)\dots(\lambda_k-z)}{\lambda_1\lambda_2\dots\lambda_k} \in \overline{\mathcal{P}}_k.
```
Puesto que ``\max_{\lambda\in\sigma(A)}|p_k(\lambda)|=0``, del Corolario 5.14 se sigue que
```math
\|x_k-x_*\|_A \leq \underbrace{\max_{\lambda\in\sigma(A)}|p_k(\lambda)|}_{=0}\|x_0-x_*\|_A = 0,
```
y así, ``x_k=x_*``, como queríamos demostrar.

**Teorema 5.22.** *Sea ``A`` sdp y ``b\in\mathbb{R}^n``. Si ``b`` es una combinación lineal de ``k`` autovectores de ``A`` entonces el método del gradiente conjugado para resolver el sistema ``Ax=b``, partiendo de ``x_0=0``, converge a la solución en a lo sumo ``k`` pasos.*

**Demostración.** La demostración se deja como ejercicio. Sugerencia: Utilizar el Teorema 5.13.

En la práctica, no es cierto en general que el método del gradiente conjugado converja en un número finito de pasos, debido a los errores de redondeo en los cálculos computacionales. En la actualidad, este método se utiliza como método iterativo para sistemas con matrices de gran tamaño (por ejemplo, ``n=100000``). En este caso, no se pretende realizar las ``n`` iteraciones para *alcanzar* la raíz. En vez de esto, se utiliza algún criterio de parada en las iteraciones, dependiendo de una tolerancia admisible para el error.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000023
md"""
## 5.5 Implementación del método del gradiente conjugado

En esta sección suponemos que ``A`` es sdp, ``b\in\mathbb{R}^n``, y ``x_*`` denota la solución de ``Ax=b``. Dado ``x_0\in\mathbb{R}^n``, consideramos ``\{x_k\}_{k\in\mathbb{N}}`` la sucesión generada por el método del gradiente conjugado.

Recordemos que
```math
x_{k+1} = \operatorname*{argmin}_{x\in x_0+\mathcal{K}_{k+1}}\phi(x), \qquad k\geq0,
```
donde ``\mathcal{K}_{k+1}=\text{span}\{r_0,r_1,\dots,r_k\}=\text{span}\{r_0,Ar_0,\dots,A^kr_0\}``.

Veremos que una vez que ``x_k`` está determinado, o bien, ``x_k=x_*``, o bien, existe una dirección de búsqueda ``d_{k+1}`` (fácil de obtener) tal que
```math
x_{k+1} = x_k+\alpha_{k+1}d_{k+1}, \qquad\text{para algún } \alpha_{k+1}\in\mathbb{R}.
```

**Lema 5.23.** *Para todo ``k\in\mathbb{N}``,*
```math
x_{k+1}-x_k\in\mathcal{K}_{k+1}, \qquad\text{y}\qquad (x_{k+1}-x_k)^TA\xi=0, \quad\forall\xi\in\mathcal{K}_k.
```

**Demostración.** Sea ``k\in\mathbb{N}``. Sabemos que ``x_k=x_0+w_k``, con ``w_k\in\mathcal{K}_k``, y que ``x_{k+1}=x_0+w_{k+1}``, con ``w_{k+1}\in\mathcal{K}_{k+1}``. Así, puesto que ``\mathcal{K}_k\subset\mathcal{K}_{k+1}``,
```math
x_{k+1}-x_k = \underbrace{w_{k+1}}_{\in\mathcal{K}_{k+1}}-\underbrace{w_k}_{\in\mathcal{K}_k} \in \mathcal{K}_{k+1}.
```
Además, si ``\xi\in\mathcal{K}_k``,
```math
(x_{k+1}-x_k)^TA\xi = \xi^TA(x_{k+1}-x_k) = \xi^T(\underbrace{Ax_{k+1}-b}_{=-r_{k+1}}+\underbrace{b-Ax_k}_{=r_k}) = \xi^T(r_k-r_{k+1}) = 0,
```
ya que tanto ``r_k`` como ``r_{k+1}`` son ortogonales a ``\mathcal{K}_k=\text{span}\{r_0,r_1,\dots,r_{k-1}\}``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000024
md"""
Para cada ``k\in\mathbb{N}``, si ``r_k\neq0``, existe una única dirección en ``\mathcal{K}_{k+1}=\text{span}\{r_0,r_1,\dots,r_k\}`` que es ``A``-ortogonal a ``\mathcal{K}_k=\text{span}\{r_0,r_1,\dots,r_{k-1}\}``, puesto que ``\mathcal{K}_k\subset\mathcal{K}_{k+1}`` y que ``\dim(\mathcal{K}_k)=k`` y ``\dim(\mathcal{K}_{k+1})=k+1``. En particular, existe un único ``w_k\in\mathcal{K}_k`` tal que ``d_{k+1}:=r_k+w_k`` es ``A``-ortogonal a ``\mathcal{K}_k``.

Así, tenemos la siguiente consecuencia del lema anterior.

**Corolario 5.24.** *Sea ``d_1=r_0`` y sea ``\{d_{k+1}\}_{k\in\mathbb{N}}`` la sucesión definida arriba. Entonces, si ``k\in\mathbb{N}_0`` y ``r_k\neq0``, se tiene que*
```math
x_{k+1} = x_k+\alpha_{k+1}d_{k+1}, \qquad\text{donde } \alpha_{k+1} = \frac{d_{k+1}^Tr_k}{d_{k+1}^TAd_{k+1}}.
```

**Demostración.** Sea ``k\in\mathbb{N}``. Por el lema anterior tenemos que ``x_{k+1}-x_k\in\mathcal{K}_{k+1}``, y ``x_{k+1}-x_k\perp_A\mathcal{K}_k``, es decir, ``x_{k+1}-x_k`` es una dirección en ``\mathcal{K}_{k+1}`` que es ``A``-ortogonal a ``\mathcal{K}_k``. Por lo tanto, existe ``\alpha_{k+1}\in\mathbb{R}`` tal que
```math
x_{k+1}-x_k = \alpha_{k+1}d_{k+1}, \qquad\text{o sea,}\qquad x_{k+1} = x_k+\alpha_{k+1}d_{k+1}.
```
Puesto que ``x_{k+1}=\operatorname*{argmin}_{x\in x_0+\mathcal{K}_{k+1}}\phi(x)``, tenemos que ``\alpha_{k+1}`` debe satisfacer que
```math
\alpha_{k+1} = \operatorname*{argmin}_{\alpha\in\mathbb{R}}\phi(x_k+\alpha d_{k+1}).
```
Teniendo en cuenta la Observación 5.4, concluimos que
```math
\alpha_{k+1} = \frac{d_{k+1}^Tr_k}{d_{k+1}^TAd_{k+1}}.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000025
md"""
Ahora estamos en condiciones de escribir la siguiente versión mejorada del algoritmo anterior:

**Método del gradiente conjugado.** (Versión preliminar 2)

```octave
Sean A sdp y b en R^n. Sea x0 en R^n.
r0 = b - A*x0;
k = 0;
mientras rk != 0
    Calcular la direccion de busqueda d{k+1};
    alpha{k+1} = (d{k+1}' * rk) / (d{k+1}' * A * d{k+1});
    x{k+1} = xk + alpha{k+1} * d{k+1};
    r{k+1} = b - A * x{k+1};
    k = k + 1;
fin
```

Antes de determinar una fórmula para calcular las direcciones de búsqueda ``d_k``, demostramos que éstas también generan el espacio de Krylov. Más precisamente, tenemos el siguiente resultado.

**Lema 5.25.** *Sea ``k\in\mathbb{N}``. Si ``r_{k-1}\neq0``, entonces ``\{d_1,\dots,d_k\}`` es ``A``-ortogonal y*
```math
\mathcal{K}_k = \text{span}\{d_1,\dots,d_k\}.
```

**Demostración.** Puesto que ``d_1=r_0``, tenemos que ``\mathcal{K}_1=\text{span}\{d_1\}``. Supongamos que
```math
\mathcal{K}_k = \text{span}\underbrace{\{d_1,\dots,d_k\}}_{A\text{-ortogonal}}, \qquad\text{(Hip. ind.)}
```
y veamos que si ``r_k\neq0``,
```math
\mathcal{K}_{k+1} = \text{span}\underbrace{\{d_1,\dots,d_{k+1}\}}_{A\text{-ortogonal}}.
```
Por hipótesis inductiva se tiene que ``\{d_1,\dots,d_k\}`` es ``A``-ortogonal y que ``\{d_1,\dots,d_k\}\subset\mathcal{K}_k\subset\mathcal{K}_{k+1}``. Además, de la definición de ``d_{k+1}`` se sigue que ``d_{k+1}\in\mathcal{K}_{k+1}``, y que ``d_{k+1}`` es ``A``-ortogonal a ``\mathcal{K}_k=\text{span}\{d_1,\dots,d_k\}``. En consecuencia,
```math
\text{span}\underbrace{\{d_1,\dots,d_{k+1}\}}_{A\text{-ortogonal}} \subset \mathcal{K}_{k+1}.
```
Puesto que ``r_k\neq0``, se tiene que ``d_1,\dots,d_{k+1}`` son no nulos, y en consecuencia, generan un subespacio de dimensión ``k+1`` del espacio de Krylov ``\mathcal{K}_{k+1}`` que tiene dimensión ``k+1``. Así,
```math
\text{span}\underbrace{\{d_1,\dots,d_{k+1}\}}_{A\text{-ortogonal}} = \mathcal{K}_{k+1},
```
como queríamos demostrar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000026
md"""
**Teorema 5.26.** *Sea ``k\in\mathbb{N}``. Si ``r_{k-1}\neq0``, entonces*
```math
d_{k+1} = r_k+\beta_{k+1}d_k,
```
*donde* ``\beta_{k+1} = -\dfrac{r_k^TAd_k}{d_k^TAd_k}``.

**Demostración.** Sabemos que existe un único ``w_k\in\mathcal{K}_k`` tal que ``d_{k+1}=r_k+w_k``. Por el lema anterior, tenemos que ``\mathcal{K}_k=\text{span}\{d_1,\dots,d_k\}`` y por lo tanto, ``w_k=\sum_{i=1}^k\gamma_id_i``, para algunos ``\gamma_1,\dots,\gamma_k\in\mathbb{R}``, y luego,
```math
d_{k+1} = r_k+\sum_{i=1}^k\gamma_id_i.
```
De la definición de ``d_{k+1}`` se sigue que ``d_{k+1}\perp_A\mathcal{K}_k=\text{span}\{d_1,\dots,d_k\}``, y por lo tanto,
```math
d_{k+1}^TAd_j = 0, \qquad j=1,\dots,k.
```
Así,
```math
d_{k+1}^TAd_j = \left(r_k^T+\sum_{i=1}^k\gamma_id_i^T\right)Ad_j = r_k^TAd_j+\sum_{i=1}^k\gamma_i\underbrace{d_i^TAd_j}_{=0,\text{ si } i\neq j} = r_k^TAd_j+\gamma_jd_j^TAd_j = 0,
```
y de aquí,
```math
\gamma_j = -\frac{r_k^TAd_j}{d_j^TAd_j}, \qquad j=1,\dots,k.
```
En consecuencia, para completar la demostración, es suficiente probar que ``\gamma_j=0``, para ``j=1,\dots,k-1``, y esto último es equivalente a probar que
```math
r_k^TAd_j = 0, \qquad j=1,\dots,k-1.
```
Finalmente, esta última afirmación se cumple ya que ``Ad_j\in\mathcal{K}_{j+1}\subset\mathcal{K}_k``, para ``j=1,\dots,k-1``, teniendo en cuenta que ``r_k\perp\{r_0,r_1,\dots,r_{k-1}\}=\mathcal{K}_k``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000027
md"""
Ahora que tenemos una fórmula para calcular las direcciones de búsqueda ``d_k`` el algoritmo queda de la siguiente manera:

**Método del gradiente conjugado.** (Versión preliminar 3)

```octave
Sean A sdp y b en R^n. Sea x0 en R^n.
r0 = b - A*x0;
k = 0;
mientras rk != 0
    si k == 0
        d1 = r0;
    si no
        beta{k+1} = -(rk' * A * dk) / (dk' * A * dk);
        d{k+1} = rk + beta{k+1} * dk;
    fin si
    alpha{k+1} = (d{k+1}' * rk) / (d{k+1}' * A * d{k+1});
    x{k+1} = xk + alpha{k+1} * d{k+1};
    r{k+1} = b - A * x{k+1};
    k = k + 1;
fin
```

En el siguiente lema damos expresiones alternativas para el cálculo de los coeficientes ``\alpha_k`` y ``\beta_k``, de manera de hacer más eficiente nuestro algoritmo.

**Lema 5.27.** *Si ``r_k\neq0``, entonces,*
```math
\alpha_{k+1} = \frac{\|r_k\|_2^2}{\|d_{k+1}\|_A^2} = \frac{r_k^Tr_k}{d_{k+1}^TAd_{k+1}}, \qquad\text{y}\qquad \beta_{k+1} = \frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2} = \frac{r_k^Tr_k}{r_{k-1}^Tr_{k-1}}.
```

**Demostración.** Notemos que
```math
d_{k+1}^Tr_k = (r_k+\beta_{k+1}d_k)^Tr_k = r_k^Tr_k+\beta_{k+1}\underbrace{d_k^Tr_k}_{=0} = r_k^Tr_k,
```
ya que ``d_k\in\mathcal{K}_k`` y ``r_k\perp\mathcal{K}_k``. En consecuencia, ``\alpha_{k+1}=\dfrac{d_{k+1}^Tr_k}{d_{k+1}^TAd_{k+1}}=\dfrac{r_k^Tr_k}{d_{k+1}^TAd_{k+1}}``.

Por otro lado, puesto que ``r_k=b-Ax_k=b-Ax_{k-1}+A(x_{k-1}-x_k)=r_{k-1}-\alpha_kAd_k`` se tiene que
```math
Ad_k = \frac{1}{\alpha_k}(r_{k-1}-r_k).
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000028
md"""
Así, usando que ``\alpha_k=\dfrac{r_{k-1}^Tr_{k-1}}{d_k^TAd_k}`` y que ``r_k\perp r_{k-1}``,
```math
\beta_{k+1} = -\frac{r_k^TAd_k}{d_k^TAd_k} = -\frac{r_k^T(r_{k-1}-r_k)}{\alpha_kd_k^TAd_k} = \frac{r_k^Tr_k}{r_{k-1}^Tr_{k-1}},
```
como queríamos demostrar.

Reemplazando las fórmulas que indica el teorema, obtenemos el siguiente algoritmo del Gradiente Conjugado, que aún no es implementable, porque los vectores y cantidades que aparecen tienen subíndices.

**Método del gradiente conjugado.** (Versión final, con subíndices)

```octave
Sean A sdp y b en R^n.
Sea x0 en R^n.
r0 = b - A*x0;
rho0 = ||r0||_2^2;
k = 0;
mientras rho_k != 0
    si k == 0
        d1 = r0;
    si no
        beta{k+1} = rho_k / rho_{k-1};      % = (rk' rk) / (r{k-1}' r{k-1})
        d{k+1} = rk + beta{k+1} * dk;
    fin si
    w{k+1} = A * d{k+1};
    alpha{k+1} = rho_k / (d{k+1}' * w{k+1});    % = (rk' rk) / (d{k+1}' A d{k+1})
    x{k+1} = xk + alpha{k+1} * d{k+1};
    r{k+1} = b - A * x{k+1};
    rho_{k+1} = ||r{k+1}||_2^2;
    k = k + 1;
fin
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000029
md"""
Antes de escribir cómo queda el algoritmo del método del gradiente conjugado hacemos una observación respecto del criterio de parada. Sabemos que
```math
\frac{\|x_k-x_*\|_2}{\|x_*\|_2} \leq \kappa_2(A)\frac{\|\overbrace{b-Ax_k}^{=r_k}\|_2}{\|b\|_2}.
```
Dada una tolerancia ``\text{tol}``, lo usual es pedir que
```math
\|r_k\|_2 < \text{tol}\,\|b\|_2,
```
y así, el error relativo a la salida del algoritmo será
```math
\frac{\|x_k-x_*\|_2}{\|x_*\|_2} \leq \kappa_2(A)\,\text{tol}.
```

Ahora sí, presentamos el algoritmo para el método del gradiente conjugado.

**Método del gradiente conjugado.**

```octave
Sean A sdp y b en R^n y sean
    x en R^n (iterando inicial),
    tol en R_+ (tolerancia para el error relativo),
    maxiter en N (tope de iteraciones).
r = b - A*x;      % residuo inicial
rho = r' * r;      % rho = ||r||_2^2
TOL = tol^2 * b' * b;      % TOL = tol^2 ||b||_2^2
k = 0;
mientras (rho > TOL y k < maxiter)      % ||r||_2/||b||_2 > tol
    si k == 0
        d = r;
    sino
        beta = rho / rho_old;
        d = r + beta * d;
    fin si
    w = A * d;
    alpha = rho / (d' * w);
    x = x + alpha * d;
    r = r - alpha * w;      % r{k+1} = b - A x{k+1} = b - A xk + A(xk - x{k+1}) = rk - alpha{k+1} A d{k+1}
    rho_old = rho;
    rho = r' * r;
    k = k + 1;
fin
```

**Observación 5.28.** (1) Notemos que, al igual que ocurría con el método del descenso más pronunciado, no es necesario almacenar la matriz ``A``, sino que es suficiente con disponer de una rutina para calcular el producto ``Ax``, para cada vector ``x``.

(2) *Almacenamiento.* Se necesitan guardar cuatro vectores: ``x``, ``w``, ``d`` y ``r``.

(3) *Costo por iteración.*
- 1 producto mat-vec (``Ad``).
- 2 productos escalares (``d^Tw`` y ``r^Tr``).
- 3 operaciones axpy, es decir, del tipo ``ax+y``, con ``a\in\mathbb{R}``, ``x,y\in\mathbb{R}^n`` (``d=r+\beta d``, ``x=x+\alpha d`` y ``r=r-\alpha w``).

(4) Notar que no es necesario almacenar una base para el espacio de Krylov (``d_1,\dots,d_k,\dots``). Esta es una ventaja que se tiene por el hecho de que ``A`` es sdp. Existen métodos de espacios de Krylov para matrices que no son sdp (por ejemplo, GMRes), pero en este caso es necesario ir guardando las direcciones de búsqueda.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000002a
md"""
## 5.6 Precondicionamiento de un sistema y método del GC precondicionado

En general, cuando queremos resolver un sistema ``Ax=b`` con un método iterativo, la convergencia depende del número de condición de ``A``, ``\kappa(A)``. Si pudiéramos determinar una matriz ``M`` invertible tal que ``M`` sea como ``A^{-1}`` en algún sentido, entonces el sistema ``Ax=b``, sería equivalente al sistema ``MAx=Mb``.

Si ``M\approx A^{-1}`` entonces ``MA\approx I``, y por lo tanto, el sistema ``MAx=Mb`` está mejor condicionado.

### 5.6.1 Gradiente conjugado precondicionado (GCP)

Vimos que si ``\kappa_2(A)`` es grande, la convergencia del método del gradiente conjugado puede ser muy lenta. Supongamos que disponemos de una matriz ``M`` sdp tal que ``M\approx A^{-1}`` en el siguiente sentido: existen constantes ``c_1,c_2`` tales que
```math
c_1 \leq \frac{x^TMx}{x^TA^{-1}x} \leq c_2, \qquad \forall x\neq0.
```
Si ``c_2/c_1`` es pequeño (cercano a 1) diremos que ``M`` es un buen precondicionador de ``A``.

El sistema ``Ax=b`` es equivalente al sistema ``MAx=Mb``. Aunque ``A`` y ``M`` son sdp, no podemos garantizar que ``MA`` sea sdp y por lo tanto, no podemos utilizar el método del gradiente conjugado para resolver el sistema ``MAx=Mb``.

Si ``S`` denota la raíz cuadrada de ``M``, es decir, si ``S`` es sdp y tal que ``S^2=M``, se tiene que
```math
Ax=b \iff MAx=Mb \iff S^2Ax=S^2b \iff SAx=Sb.
```
Así, resolver el sistema ``Ax=b`` es equivalente a resolver
```math
\begin{cases}SASy=Sb\\x=Sy.\end{cases}
```

Se deja como ejercicio verificar que ``SAS`` es sdp, y por lo tanto, podemos aplicar el método del gradiente conjugado para resolver el sistema ``SASy=Sb``. Más aún, se puede demostrar (ejercicio!) que
```math
\kappa_2(SAS) \leq \frac{c_2}{c_1}.
```
Si ``\frac{c_2}{c_1}`` es pequeño, la convergencia del método del gradiente conjugado aplicado al sistema ``SASy=Sb`` será más rápida. Pero no conocemos y no queremos calcular ``S=M^{1/2}``!
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000002b
md"""
Supongamos por el momento que disponemos de la matriz ``S`` y apliquemos el método del gradiente conjugado para resolver el sistema ``SASy=Sb``. Llamemos ``\{y_k\}_{k\in\mathbb{N}}`` a la sucesión generada por dicho algoritmo, y denotemos por ``\tilde r_k:=Sb-SASy_k`` al residuo en el paso ``k``, y por ``\tilde d_{k+1}`` a la dirección de búsqueda para obtener ``y_{k+1}`` a partir de ``y_k``.

**Método del gradiente conjugado aplicado a** ``SASy=Sb``.

```octave
Sea y0 en R^n.
r0_tilde = S*b - S*A*S*y0;
rho0_tilde = ||r0_tilde||_2^2;
k = 0;
mientras rho_k_tilde != 0
    si k == 0
        d1_tilde = r0_tilde;
    sino
        beta{k+1}_tilde = rho_k_tilde / rho_{k-1}_tilde;
        d{k+1}_tilde = r_k_tilde + beta{k+1}_tilde * dk_tilde;
    fin si
    w{k+1}_tilde = S*A*S*d{k+1}_tilde;
    alpha{k+1}_tilde = rho_k_tilde / (d{k+1}_tilde' * w{k+1}_tilde);
    y{k+1} = yk + alpha{k+1}_tilde * d{k+1}_tilde;
    r{k+1}_tilde = r_k_tilde - alpha{k+1}_tilde * w{k+1}_tilde;
    rho_{k+1}_tilde = ||r{k+1}_tilde||_2^2;
    k = k + 1;
fin
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000002c
md"""
Definamos ahora ``x_k:=Sy_k``. Si ``y_k`` converge a la solución ``y_*`` de ``SASy=Sb``, entonces ``x_k`` converge a la solución ``x_*`` de ``Ax=b``. Más concretamente, puesto que ``\kappa_2(SAS)\leq\frac{c_2}{c_1}``, tenemos que
```math
\|y_k-y_*\|_{SAS} \leq 2\left(\frac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\right)^k\|y_0-y_*\|_{SAS}, \qquad k\geq0. \tag{5.11}
```
Además, ya que ``\|z\|_{SAS}^2=z^TSASz=(Sz)^TA(Sz)=\|Sz\|_A^2``, para todo ``z\in\mathbb{R}^n``, la desigualdad anterior nos da
```math
\|x_k-x_*\|_A \leq 2\left(\frac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\right)^k\|x_0-x_*\|_A, \qquad k\geq0. \tag{5.12}
```
Por lo tanto, si ``M=S^2`` es un buen precondicionador de ``A``, esta nueva sucesión ``x_k`` converge más rápidamente a la solución ``x_*`` del sistema ``Ax=b`` que la sucesión generada por el método del gradiente conjugado aplicado directamente al sistema ``Ax=b``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000002d
md"""
A continuación, reescribiremos el algoritmo anterior pensando en calcular la sucesión ``x_k:=Sy_k`` en vez de ``y_k``. Para esto, definimos el residuo ``r_k:=b-Ax_k`` y notamos que
```math
r_k = b-Ax_k = b-ASy_k = S^{-1}(Sb-SASy_k) = S^{-1}\tilde r_k.
```
Por otro lado, puesto que ``y_{k+1}=y_k+\tilde\alpha_{k+1}\tilde d_{k+1}`` se tiene que
```math
Sy_{k+1} = Sy_k+\tilde\alpha_{k+1}S\tilde d_{k+1},
```
o equivalentemente, que
```math
x_{k+1} = x_k+\tilde\alpha_{k+1}d_{k+1},
```
donde la dirección de búsqueda ``d_{k+1}`` está dada por ``d_{k+1}:=S\tilde d_{k+1}``. Así,
```math
\begin{cases}
d_1 = S\tilde d_1 = S\tilde r_0 = S^2r_0 = Mr_0,\\
d_{k+1} = S\tilde d_{k+1} = S(\tilde r_k+\tilde\beta_{k+1}\tilde d_k) = S\tilde r_k+\tilde\beta_{k+1}d_k = S^2r_k+\tilde\beta_{k+1}d_k = Mr_k+\tilde\beta_{k+1}d_k, \quad k\geq1.
\end{cases}
```

Además, para los coeficientes ``\tilde\alpha_{k+1}`` y ``\tilde\beta_{k+1}`` se tiene que
```math
\tilde\alpha_{k+1} = \frac{\tilde r_k^T\tilde r_k}{\tilde d_{k+1}^TSAS\tilde d_{k+1}} = \frac{(Sr_k)^T(Sr_k)}{(S\tilde d_{k+1})^TA(S\tilde d_{k+1})} = \frac{r_k^TSSr_k}{d_{k+1}^TAd_{k+1}} = \frac{r_k^TMr_k}{d_{k+1}^TAd_{k+1}},
```
y
```math
\tilde\beta_{k+1} = \frac{\tilde r_k^T\tilde r_k}{\tilde r_{k-1}^T\tilde r_{k-1}} = \frac{r_k^TMr_k}{r_{k-1}^TMr_{k-1}}.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000002e
md"""
Antes de escribir cómo queda el algoritmo, hacemos una observación sobre el criterio de parada. Puesto que
```math
\frac{\|y_k-y_*\|_2}{\|y_*\|_2} \leq \kappa_2(SAS)\frac{\|\tilde r_k\|_2}{\|Sb\|_2},
```
se tiene que
```math
\frac{\|x_k-x_*\|_{M^{-1}}}{\|x_*\|_{M^{-1}}} \leq \kappa_2(SAS)\frac{\|r_k\|_M}{\|b\|_M}.
```
En general, se suele parar el algoritmo cuando
```math
\frac{\|r_k\|_M}{\|b\|_M} < \text{tol},
```
donde ``\text{tol}`` es una tolerancia dada.

Ahora bien, el criterio de parada controla el error en la norma ``\|\cdot\|_{M^{-1}}``, mientras que toda la teoría de convergencia del método está escrita en la norma ``\|\cdot\|_A``. El siguiente lema muestra que, si ``M`` es un buen precondicionador, ambas normas son *equivalentes con constantes cercanas entre sí*, de modo que controlar una es esencialmente lo mismo que controlar la otra.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000002f
md"""
**Lema 5.29.** *Sean ``A`` y ``M`` sdp, y supongamos que existen ``0<c_1\leq c_2`` tales que*
```math
c_1 \leq \frac{x^TMx}{x^TA^{-1}x} \leq c_2, \qquad \forall x\in\mathbb{R}^n\setminus\{0\}.
```
*Entonces, para todo* ``v\in\mathbb{R}^n``,
```math
\frac{1}{c_2}\|v\|_A^2 \leq \|v\|_{M^{-1}}^2 \leq \frac{1}{c_1}\|v\|_A^2,
```
*o equivalentemente,*
```math
\sqrt{c_1}\|v\|_{M^{-1}} \leq \|v\|_A \leq \sqrt{c_2}\|v\|_{M^{-1}}. \tag{5.13}
```

**Demostración.** Observemos primero que si ``B`` es sdp, entonces para todo ``v\in\mathbb{R}^n``
```math
v^TB^{-1}v = \max_{x\neq0}\frac{(x^Tv)^2}{x^TBx}. \tag{5.14}
```
En efecto, escribiendo ``x^Tv=(B^{1/2}x)^T(B^{-1/2}v)``, la desigualdad de Cauchy-Schwarz da ``(x^Tv)^2\leq(x^TBx)(v^TB^{-1}v)`` para todo ``x``, y la igualdad se alcanza tomando ``x=B^{-1}v``.

Sea ahora ``v\in\mathbb{R}^n``, ``v\neq0``. Aplicando (5.14) con ``B=M`` y con ``B=A^{-1}`` obtenemos
```math
\|v\|_{M^{-1}}^2 = v^TM^{-1}v = \max_{x\neq0}\frac{(x^Tv)^2}{x^TMx}, \qquad \|v\|_A^2 = v^TAv = \max_{x\neq0}\frac{(x^Tv)^2}{x^TA^{-1}x}.
```
Por hipótesis, ``c_1\,x^TA^{-1}x\leq x^TMx\leq c_2\,x^TA^{-1}x`` para todo ``x\neq0``, y por lo tanto
```math
\frac{1}{c_2}\cdot\frac{(x^Tv)^2}{x^TA^{-1}x} \leq \frac{(x^Tv)^2}{x^TMx} \leq \frac{1}{c_1}\cdot\frac{(x^Tv)^2}{x^TA^{-1}x}, \qquad \forall x\neq0.
```
Tomando máximo sobre ``x\neq0`` en los tres miembros resulta
```math
\frac{1}{c_2}\|v\|_A^2 \leq \|v\|_{M^{-1}}^2 \leq \frac{1}{c_1}\|v\|_A^2,
```
que es lo que queríamos probar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000030
md"""
**Observación 5.30.** La equivalencia (5.13) tiene "ancho" ``\sqrt{c_2/c_1}``, que es precisamente la cantidad que mide la calidad del precondicionador. En el caso extremo ``M=A^{-1}`` se tiene ``c_1=c_2=1`` y ``M^{-1}=A``: las dos normas coinciden. Si ``\frac{c_2}{c_1}\approx1``, las normas ``\|\cdot\|_{M^{-1}}`` y ``\|\cdot\|_A`` son casi indistinguibles.

Combinando (5.13) (aplicada a ``v=x_k-x_*`` y a ``v=x_*``) con la estimación anterior y con ``\kappa_2(SAS)\leq\frac{c_2}{c_1}``, obtenemos una cota para el error relativo *en la norma ``A``* al salir del algoritmo:
```math
\frac{\|x_k-x_*\|_A}{\|x_*\|_A} \leq \sqrt{\frac{c_2}{c_1}}\frac{\|x_k-x_*\|_{M^{-1}}}{\|x_*\|_{M^{-1}}} \leq \sqrt{\frac{c_2}{c_1}}\,\kappa_2(SAS)\frac{\|r_k\|_M}{\|b\|_M} \leq \left(\frac{c_2}{c_1}\right)^{3/2}\text{tol}.
```
Es decir: el criterio de parada ``\frac{\|r_k\|_M}{\|b\|_M}<\text{tol}`` garantiza un error relativo en norma ``A`` del orden de ``\text{tol}``, con una constante que *no* depende de ``\kappa_2(A)`` sino solamente de la calidad del precondicionador. Comparar con lo que ocurre en el gradiente conjugado sin precondicionar, donde la constante es ``\kappa_2(A)``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000031
md"""
**Método del gradiente conjugado precondicionado.**

```octave
Sean A sdp y b en R^n y sean
    M sdp (precondicionador de A, i.e., M aproximadamente A^{-1}),
    x en R^n (iterando inicial),
    tol en R_+ (tolerancia para el error relativo en norma M),
    maxiter en N (tope de iteraciones).
r = b - A*x;      % residuo inicial
z = M*r,
rho = r' * z;      % rho = ||r||_M^2
TOL = tol^2 * b' * M * b;      % TOL = tol^2 ||b||_M^2
k = 0;
mientras (rho > TOL y k < maxiter)      % ||r||_M/||b||_M > tol
    si k == 0
        d = z;      % d = Mr
    sino
        beta = rho / rho_old;
        d = z + beta * d;
    fin si
    w = A * d;
    alpha = rho / (d' * w);
    x = x + alpha * d;
    r = r - alpha * w;
    z = M * r;
    rho_old = rho;
    rho = r' * z;
    k = k + 1;
fin
```

**Observación 5.31.** (1) Notemos que finalmente para implementar el GCP no se necesita conocer (calcular) la raíz cuadrada ``S`` del precondicionador ``M``. Esto fue una herramienta teórica en la deducción del método.

(2) El costo de este algoritmo es solo un producto mat-vec más por iteración que el gradiente conjugado sin precondicionar.

(3) No es necesario conocer ``M``, basta con disponer de una rutina ``\xi\to M\xi``.

En general, un *buen precondicionador* para una matriz ``A`` sdp es una matriz ``M`` sdp tal que

- La acción ``\xi\to M\xi`` es "fácil" de calcular.
- ``M\approx A^{-1}`` en el siguiente sentido: existen constantes ``c_1,c_2`` tales que
```math
c_1 \leq \frac{x^TMx}{x^TA^{-1}x} \leq c_2, \qquad \forall x\neq0, \qquad\text{con } \frac{c_2}{c_1}\approx1.
```

Así, la sucesión ``\{x_k\}_{k\in\mathbb{N}_0}`` generada por el gradiente conjugado precondicionado satisface
```math
\|x_k-x_*\|_A \leq 2\left(\frac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\right)^k\|x_0-x_*\|_A, \qquad k\geq0.
```

Si tomáramos ``M=I``, la acción ``\xi\to M\xi`` es trivial de calcular, pero ``\frac{c_2}{c_1}`` puede ser muy grande. Por otro lado, la "mejor elección" para satisfacer la segunda condición sería tomar ``M=A^{-1}`` (y así, ``c_1=c_2=1``), pero la acción ``\xi\to M\xi`` es tan costosa de calcular como resolver el problema original!

En general, para construir un buen precondicionador ``M`` para una matriz ``A`` hay que lograr un adecuado equilibrio entre ambas condiciones, y es necesario saber *de dónde proviene* la matriz ``A`` (por ejemplo, de la discretización de una ecuación diferencial) y usar esa información para diseñar ``M``.

El método del gradiente conjugado precondicionado está implementado en MATLAB/OCTAVE a través de la función `pcg`. Es importante tener en cuenta que para MATLAB/OCTAVE, el precondicionador ``M`` es una aproximación de ``A``, y no de ``A^{-1}``, como hemos supuesto. En este caso, las condiciones dadas para un buen precondicionador deben reemplazarse por las siguientes:

- Es "fácil" resolver sistemas con matriz ``M`` (es decir, la acción ``\xi\to M^{-1}\xi`` es "fácil" de calcular).
- ``M\approx A`` en el siguiente sentido: existen constantes ``c_1,c_2`` tales que
```math
c_1 \leq \frac{x^TMx}{x^TAx} \leq c_2, \qquad \forall x\neq0, \qquad\text{con } \frac{c_2}{c_1}\approx1.
```
"""

# ╔═╡ e5f6a7b8-dddd-4a1b-8c2d-000000000002
md"""
**Nota** (de las diapositivas de la Clase 11). Puede ser útil comparar, lado a lado, el algoritmo del gradiente conjugado sin precondicionar y el precondicionado:

| | GC | GCP |
|---|---|---|
| Inicialización | ``r=b-Ax``, ``\rho=r^Tr`` | ``r=b-Ax``, ``z=Mr``, ``\rho=r^Tz`` |
| Dirección | ``d=r+\beta d`` | ``d=z+\beta d`` |
| Actualizar ``r`` | ``r=r-\alpha w`` | ``r=r-\alpha w``, luego ``z=Mr`` |
| Escalar | ``\rho=r^Tr=\|r\|_2^2`` | ``\rho=r^Tz=\|r\|_M^2`` |
| Parada | ``\|r\|_2<\text{tol}\,\|b\|_2`` | ``\|r\|_M<\text{tol}\,\|b\|_M`` |
| Memoria | 4 vectores | 5 vectores |
| Costo por iteración | 1 mat-vec (``Ad``) | 1 mat-vec más 1 acción ``M\xi`` |
| Tasa | ``\left(\dfrac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\right)^k`` | ``\left(\dfrac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\right)^k`` |

El GC es el caso particular ``M=I`` del GCP: ahí ``z=r`` y todo colapsa al algoritmo sin precondicionar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000032
md"""
### 5.6.2 Precondicionadores genéricos

Existen algunos precondicionadores *genéricos*, es decir, que se construyen para cualquier ``A``, pero en estos casos, no deberíamos esperar un cambio *dramático* en la convergencia de un método iterativo aplicado al sistema precondicionado respecto del método aplicado al sistema sin precondicionar.

Consideremos el sistema
```math
Ax=b.
```
Si ``M`` es una matriz invertible, entonces
```math
Ax=b \iff MAx=Mb \iff x=(I-MA)x+Mb.
```
Teniendo en cuenta que la iteración
```math
x_{k+1} = (I-MA)x_k+Mb,
```
converge si y solo si ``\rho(I-MA)<1``, a veces se elige la matriz ``M`` de manera que
```math
I-MA=B,
```
donde ``B`` es la matriz de iteración de un método iterativo estacionario.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000033
md"""
**Precondicionamiento diagonal.**

Para el método de Jacobi, la matriz de iteración viene dada por
```math
B_J = -D^{-1}(L+U),
```
donde ``D``, ``L`` y ``U`` son la parte diagonal, la parte triangular inferior estricta, y la parte triangular superior estricta de ``A``, respectivamente. Buscamos ``M_J`` tal que ``I-M_JA=B_J``, o equivalentemente,
```math
M_JA = I+D^{-1}(L+U) = D^{-1}(\underbrace{D+L+U}_{=A}).
```
Luego,
```math
M_J = D^{-1}.
```
Esta elección para la matriz ``M`` se conoce como *precondicionamiento diagonal o de Jacobi*. Notemos que la acción ``\xi\to M_J\xi`` es fácil de realizar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000034
md"""
**Precondicionador de Gauss-Seidel simétrico.**

Recordemos que la matriz de iteración del método de Gauss-Seidel viene dada por
```math
B_{GS} = -(L+D)^{-1}U.
```
Por otro lado, el *método de Gauss-Seidel hacia atrás* tiene la matriz de iteración
```math
B_{BGS} = -(U+D)^{-1}L.
```
Finalmente, el *método de Gauss-Seidel simétrico* consiste en aplicar una vez Gauss-Seidel (hacia adelante) y una vez Gauss-Seidel hacia atrás. Así, la matriz de iteración en este caso viene dada por
```math
B_{SGS} = B_{BGS}B_{GS} = (U+D)^{-1}L(L+D)^{-1}U.
```
Si ``A`` es simétrica, se tiene que ``U=L^T`` y
```math
B_{SGS} = (L^T+D)^{-1}L(L+D)^{-1}L^T.
```
En este caso, el precondicionador de Gauss-Seidel simétrico está definido por
```math
M_{SGS} = (L^T+D)^{-1}D(L+D)^{-1},
```
ya que
```math
\begin{align*}
I-M_{SGS}A &= I-(L^T+D)^{-1}D(L+D)^{-1}(L+D+L^T)\\
&= (L^T+D)^{-1}(L^T+D-D(I+(L+D)^{-1}L^T))\\
&= (L^T+D)^{-1}(L^T-D(L+D)^{-1}L^T)\\
&= (L^T+D)^{-1}(\underbrace{L^T((L+D)^{-1}L^T)^{-1}-D}_{=L})(L+D)^{-1}L^T\\
&= B_{SGS}.
\end{align*}
```
Notar que ``\xi\to M_{SGS}\xi`` consiste en resolver un sistema triangular inferior, multiplicar una matriz diagonal por un vector y finalmente resolver un sistema triangular superior.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000035
md"""
**Precondicionadores de factorización incompleta.**

Supongamos que ``A=R^TR+E``, donde ``R`` es una matriz triangular superior y ``E`` es "pequeña". Podemos tomar
```math
M = (R^TR)^{-1} = R^{-1}(R^T)^{-1},
```
observando que la acción ``\xi\to M\xi`` consiste en resolver un sistema triangular inferior (con sustitución hacia adelante) seguido de un sistema triangular superior (con sustitución hacia atrás).

La función `cholinc` de MATLAB/OCTAVE calcula factorizaciones incompletas de matrices.
"""

# ╔═╡ e5f6a7b8-dddd-4a1b-8c2d-000000000003
md"""
**Nota** (de las diapositivas de la Clase 11). Un resumen de los precondicionadores genéricos de esta sección:

| Precondicionador | ``M`` | Acción ``\xi\to M\xi`` | Preparación |
|---|---|---|---|
| Ninguno | ``I`` | gratis | ninguna |
| Diagonal (Jacobi) | ``D^{-1}`` | ``n`` divisiones | ninguna |
| Gauss-Seidel simétrico | ``(L^T+D)^{-1}D(L+D)^{-1}`` | 2 sustituciones triangulares y 1 producto diagonal | ninguna |
| Cholesky incompleta | ``(R^TR)^{-1}`` | 2 sustituciones triangulares | factorizar una vez |
| Basado en la estructura | (multigrid, etc.) | ``O(n)`` | depende del problema |

Los precondicionadores genéricos mejoran la constante ``c_2/c_1``, pero rara vez el orden de su dependencia en ``h``: si ``\kappa_2(A)=O(h^{-2})``, en general sigue valiendo ``c_2/c_1=Ch^{-2}`` con ``C<1``, es decir, se gana una constante pero el exponente de ``h`` sigue ahí. Para lograr ``c_2/c_1=O(1)`` independientemente de ``h`` hace falta un precondicionador que use *de dónde viene* la matriz, como el de la sección siguiente.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000036
md"""
## 5.7 Un precondicionador basado en la estructura: multigrid en 1D

Los precondicionadores de la sección anterior son *genéricos*: se construyen a partir de las entradas de ``A``, sin usar *de dónde proviene* la matriz. Como anticipamos, los precondicionadores más efectivos suelen aprovechar el origen del sistema. En esta sección presentamos, sobre un problema modelo unidimensional, las ideas del método *multigrid* (o de *multimalla*), que da lugar tanto a un método iterativo con costo óptimo como a uno de los mejores precondicionadores conocidos para matrices provenientes de la discretización de ecuaciones diferenciales.

**Observación 5.32** (Qué hace falta saber para leer esta sección). El problema modelo se obtiene discretizando una ecuación diferencial, y por eso la primera subsección explica de dónde sale la matriz. Esa explicación es *autocontenida*: solo se usa la definición de derivada, y no se necesita ningún resultado sobre ecuaciones diferenciales. De ahí en adelante, todo lo que hacemos es álgebra lineal con una matriz concreta: el lector que prefiera saltear la deducción puede tomar el sistema (5.16) como dato y seguir leyendo desde la propiedad de suavizado sin más.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000037
md"""
### 5.7.1 El problema modelo: de dónde sale la matriz

Supongamos que queremos determinar una función ``u`` definida en el intervalo ``[0,1]``, de la que sabemos dos cosas: cuánto vale en los extremos,
```math
u(0) = u(1) = 0,
```
y cuánto vale su derivada segunda en cada punto interior,
```math
-u''(x) = f(x), \qquad x\in(0,1),
```
donde ``f`` es una función conocida. Un problema de este tipo (una ecuación que involucra derivadas de la incógnita, más datos en los extremos del intervalo) se llama *problema de valores de borde*. Es el modelo más simple de una familia de problemas que aparecen constantemente en aplicaciones: la temperatura de estado estacionario de una barra cuyos extremos se mantienen a temperatura 0 y que recibe calor a razón ``f``, o la forma que toma una cuerda elástica tensa, sujeta en sus extremos y sometida a una carga ``f``, entre muchos otros.

Para lo que sigue *no necesitamos saber nada sobre este tipo de problemas*: ni si tienen solución, ni cómo son sus soluciones, ni ninguna técnica para resolverlos. Solo vamos a transformarlo en un sistema lineal, y a partir de ahí trabajaremos exclusivamente con ese sistema.

**De la derivada al cociente incremental.** La idea es reemplazar las derivadas por cocientes incrementales, que es de donde provienen. Recordemos que
```math
u'(x) = \lim_{h\to0}\frac{u(x+h)-u(x)}{h},
```
de modo que, para ``h`` pequeño, ``u'(x)\approx\frac{u(x+h)-u(x)}{h}``: la derivada es aproximadamente la pendiente del segmento que une dos puntos cercanos del gráfico de ``u``.

Introducimos entonces una malla uniforme de ``[0,1]`` con paso ``h=\frac{1}{n+1}`` y nodos
```math
x_i = ih, \qquad i=0,1,\dots,n+1,
```
de modo que ``x_0=0`` y ``x_{n+1}=1`` son los extremos del intervalo, y ``x_1,\dots,x_n`` son los nodos interiores. Llamamos ``u_i`` a la aproximación de ``u(x_i)`` que vamos a calcular. Por el dato del borde, ``u_0=u_{n+1}=0``.

Para aproximar ``u''`` aplicamos la idea dos veces. Las pendientes a la derecha y a la izquierda del nodo ``x_i`` son, respectivamente,
```math
\frac{u_{i+1}-u_i}{h} \approx u'\left(x_i+\frac h2\right), \qquad \frac{u_i-u_{i-1}}{h} \approx u'\left(x_i-\frac h2\right),
```
y ``u''`` es la derivada de ``u'``, es decir, la razón de cambio de esas pendientes:
```math
u''(x_i) \approx \frac1h\left(\frac{u_{i+1}-u_i}{h}-\frac{u_i-u_{i-1}}{h}\right) = \frac{u_{i-1}-2u_i+u_{i+1}}{h^2}. \tag{5.15}
```
Esta aproximación se llama *diferencia centrada de segundo orden*; desarrollando ``u(x_i\pm h)`` por Taylor alrededor de ``x_i`` se verifica (ejercicio) que el error cometido es ``O(h^2)``, o sea que la aproximación mejora rápidamente al achicar ``h``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000038
md"""
**El sistema lineal.** Imponiendo ``-u''(x_i)=f(x_i)`` en cada nodo interior y usando (5.15), obtenemos
```math
\frac{-u_{i-1}+2u_i-u_{i+1}}{h^2} = f_i, \qquad i=1,\dots,n, \qquad f_i:=f(x_i),
```
con la convención ``u_0=u_{n+1}=0``. Esto es un sistema de ``n`` ecuaciones lineales con ``n`` incógnitas ``u_1,\dots,u_n``:
```math
A_hu=f, \qquad A_h = \frac{1}{h^2}\begin{pmatrix}2&-1&&\\-1&2&-1&\\&\ddots&\ddots&\ddots\\&&-1&2\end{pmatrix}\in\mathbb{R}^{n\times n}, \tag{5.16}
```
donde ``u=(u_1,\dots,u_n)^T`` y ``f=(f_1,\dots,f_n)^T``. A partir de aquí olvidamos el problema continuo: nuestro objetivo es resolver el sistema (5.16).

**Observación 5.33.** La matriz ``\text{tridiag}(-1,2,-1)`` no es una curiosidad de este problema: salvo constantes, es la misma que aparece al calcular splines cúbicos interpolantes, mencionados al comienzo del Capítulo 2, y reaparece en incontables aplicaciones. Es el ejemplo canónico de matriz *rala*, *sdp* y *mal condicionada*, que es exactamente la combinación que motivó este capítulo.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000039
md"""
**Autovalores y autovectores.** La matriz ``A_h`` es claramente simétrica. Lo notable de este problema modelo es que conocemos *exactamente* sus autovalores y autovectores. Para ``k=1,\dots,n`` definimos ``\varphi^{(k)}\in\mathbb{R}^n`` por
```math
\varphi_i^{(k)} = \text{sen}\left(\frac{ki\pi}{n+1}\right) = \text{sen}(k\pi x_i), \qquad i=1,\dots,n, \tag{5.17}
```
es decir, ``\varphi^{(k)}`` es el vector que se obtiene *muestreando* la función ``\text{sen}(k\pi x)`` en los nodos interiores de la malla. Usando la identidad ``\text{sen}(\theta-\psi)+\text{sen}(\theta+\psi)=2\,\text{sen}\,\theta\cos\psi`` con ``\theta=k\pi x_i`` y ``\psi=k\pi h``, se verifica que
```math
(A_h\varphi^{(k)})_i = \frac{1}{h^2}\left(-\varphi_{i-1}^{(k)}+2\varphi_i^{(k)}-\varphi_{i+1}^{(k)}\right) = \frac{2}{h^2}\left(1-\cos(k\pi h)\right)\varphi_i^{(k)}
```
(la fórmula vale también en ``i=1`` y en ``i=n``, porque ``\varphi_0^{(k)}=\text{sen}(0)=0`` y ``\varphi_{n+1}^{(k)}=\text{sen}(k\pi)=0``). Es decir, ``\varphi^{(k)}`` es autovector de ``A_h`` con autovalor
```math
\lambda_k = \frac{2}{h^2}(1-\cos(k\pi h)) = \frac{4}{h^2}\,\text{sen}^2\left(\frac{k\pi h}{2}\right), \qquad k=1,\dots,n. \tag{5.18}
```
Como todos los ``\lambda_k`` son positivos, ``A_h`` es sdp. Además, los autovalores son distintos dos a dos, de modo que los ``n`` autovectores ``\varphi^{(1)},\dots,\varphi^{(n)}`` son ortogonales y forman una base de ``\mathbb{R}^n``: todo vector de ``\mathbb{R}^n`` se escribe de manera única como combinación lineal de los ``\varphi^{(k)}``. Esta observación es la herramienta central de toda la sección.

**Vectores suaves y vectores oscilatorios.** Conviene *ver* los vectores de ``\mathbb{R}^n``: a un vector ``v=(v_1,\dots,v_n)^T`` le asociamos los puntos ``(x_i,v_i)`` del plano, unidos por segmentos. Con esa imagen, la familia ``\{\varphi^{(k)}\}`` se interpreta de manera muy natural: el índice ``k`` mide cuántas oscilaciones hace el vector a lo largo de la malla, y lo llamamos la *frecuencia* del modo. Los modos con ``k`` chico se llaman *suaves* (sus componentes varían poco de un nodo al siguiente) y los de ``k`` grande, *oscilatorios* (Figura 5.2). Notemos que la frecuencia es relativa a la malla: ``k`` "chico" o "grande" significa chico o grande *comparado con* ``n``.

**El número de condición.** Los autovalores (5.18) crecen con ``k``. Para ``h`` pequeño,
```math
\lambda_1 = \frac{4}{h^2}\text{sen}^2\left(\frac{\pi h}{2}\right) \approx \pi^2, \qquad \lambda_n = \frac{4}{h^2}\text{sen}^2\left(\frac{n\pi h}{2}\right) \approx \frac{4}{h^2},
```
de modo que
```math
\kappa_2(A_h) = \frac{\lambda_n}{\lambda_1} \approx \frac{4}{\pi^2h^2} = O(h^{-2}) = O(n^2). \tag{5.19}
```
El número de condición crece como ``h^{-2}``: al refinar la malla, el sistema se vuelve cada vez peor condicionado. Por (5.19) y la Observación 5.20, el método del gradiente conjugado necesita del orden de ``\sqrt{\kappa_2(A_h)}=O(h^{-1})=O(n)`` iteraciones para reducir el error a una tolerancia fija, y como cada iteración cuesta ``O(n)`` (la matriz es rala), el costo total es ``O(n^2)``. El objetivo de multigrid es resolver (5.16) con costo *óptimo* ``O(n)``.

**Figura 5.2** (descripción, no reproducida gráficamente). Dos autovectores de ``A_h`` para ``n=15``, dibujados como los puntos ``(x_i,\varphi_i^{(k)})`` unidos por segmentos. Arriba, ``\varphi^{(2)}``, un modo suave. Abajo, ``\varphi^{(14)}``, uno oscilatorio, que cambia de signo en casi cada nodo.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000003a
md"""
### 5.7.2 La propiedad de suavizado

La primera pieza de multigrid es un método iterativo estacionario simple, que usaremos *no* para resolver el sistema sino para *suavizar* el error. Consideremos el método de *Jacobi amortiguado* con parámetro ``\omega\in(0,1]``: siendo ``D=\frac{2}{h^2}I`` la diagonal de ``A_h``, la iteración es
```math
u_{m+1} = u_m+\omega D^{-1}(f-A_hu_m) = u_m+\omega\frac{h^2}{2}r_m, \tag{5.20}
```
donde ``r_m=f-A_hu_m`` es el residuo. (Con ``\omega=1`` se recupera el método de Jacobi del Capítulo 4: en efecto, la iteración de Jacobi ``u_{m+1}=D^{-1}[f-(L+U)u_m]`` se reescribe, sumando y restando ``Du_m``, como ``u_{m+1}=u_m+D^{-1}(f-A_hu_m)``. El parámetro ``\omega`` simplemente acorta el paso.) Si ``u_*`` es la solución exacta de (5.16) y ``e_m=u_m-u_*`` es el error, restando ``u_*=u_*+\omega D^{-1}(f-A_hu_*)`` de (5.20) obtenemos
```math
e_{m+1} = R_\omega e_m, \qquad R_\omega := I-\omega D^{-1}A_h = I-\omega\frac{h^2}{2}A_h.
```

La clave está en que ``R_\omega`` tiene los *mismos* autovectores que ``A_h``: como ``A_h\varphi^{(k)}=\lambda_k\varphi^{(k)}``, resulta
```math
R_\omega\varphi^{(k)} = \mu_k\varphi^{(k)}, \qquad \mu_k = 1-\omega\frac{h^2}{2}\lambda_k = 1-2\omega\,\text{sen}^2\left(\frac{k\pi h}{2}\right), \qquad k=1,\dots,n.
```
Como los ``\varphi^{(k)}`` son base de ``\mathbb{R}^n``, podemos escribir el error inicial como ``e_0=\sum_{k=1}^nc_k\varphi^{(k)}``, y entonces
```math
e_m = R_\omega^me_0 = \sum_{k=1}^nc_k\,\mu_k^m\,\varphi^{(k)}. \tag{5.21}
```
O sea: cada modo del error se amortigua, de manera independiente de los demás, por un factor ``\mu_k`` en cada iteración. Analicemos ``\mu_k`` según la frecuencia:

- **Modos suaves** (``k`` chico): ``\mu_k=1-2\omega\,\text{sen}^2(k\pi h/2)\approx1-2\omega\left(\frac{k\pi h}{2}\right)^2\approx1``. Estos modos casi no se reducen.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000003b
md"""
- **Modos oscilatorios** (``k`` grande, digamos ``\frac{n+1}{2}\leq k\leq n``): aquí ``\frac{k\pi h}{2}\in[\frac\pi4,\frac\pi2]``, con lo cual ``\text{sen}^2(k\pi h/2)\in[\frac12,1]`` y los ``\mu_k`` quedan lejos de 1.

Esto explica, de paso, por qué los métodos estacionarios del Capítulo 4 convergen tan lentamente para esta matriz: el radio espectral de ``R_\omega`` es ``\max_k|\mu_k|\approx1``, y está gobernado por los modos suaves.

Definimos el *factor de suavizado* como el peor factor de amortiguamiento sobre los modos oscilatorios,
```math
\mu_{\text{suav}}(\omega) := \max_{\frac{n+1}{2}\leq k\leq n}|\mu_k| = \max_{\frac\pi4\leq\theta\leq\frac\pi2}\left|1-2\omega\,\text{sen}^2\theta\right|.
```
Este máximo es fácil de calcular: llamando ``s=\text{sen}^2\theta``, cuando ``\theta`` recorre ``[\frac\pi4,\frac\pi2]`` el valor de ``s`` recorre ``[\frac12,1]``, y ``1-2\omega s`` es una función lineal decreciente de ``s``; su módulo alcanza el máximo en alguno de los extremos, de donde
```math
\mu_{\text{suav}}(\omega) = \max\left\{|1-\omega|,\ |1-2\omega|\right\}.
```
Para ``\omega\in(0,1]`` el primer término crece y el segundo decrece al aumentar ``\omega``, de modo que el mínimo se alcanza cuando ambos coinciden: ``1-\omega=2\omega-1``, es decir,
```math
\omega = \frac23, \qquad \mu_{\text{suav}}\left(\frac23\right) = \frac13,
```
y esto ocurre *independientemente de* ``h``.

Es decir: con ``\omega=\frac23``, cada barrido de Jacobi amortiguado reduce todos los modos oscilatorios del error al menos en un factor ``\frac13``, sin importar cuán fina sea la malla (Figura 5.3). En cambio, los modos suaves sobreviven casi intactos. Por eso decimos que la relajación *no resuelve* el sistema, pero *suaviza* el error: luego de unos pocos barridos, ``e_m`` es un vector suave, como ilustra la Figura 5.4.

**Figura 5.3** (descripción, no reproducida gráficamente). Factor de amortiguamiento ``\mu_k=1-2\omega\,\text{sen}^2(k\pi h/2)`` del Jacobi amortiguado con ``\omega=\frac23``, en función de la frecuencia ``\theta=\frac{k\pi h}{2}``. En la región de altas frecuencias (``\theta\in[\frac\pi4,\frac\pi2]``) se cumple ``|\mu_k|\leq\frac13``; en cambio, los modos suaves (``\theta`` cerca de 0) tienen ``\mu_k\approx1`` y casi no se reducen.

**Figura 5.4** (descripción, no reproducida gráficamente). Efecto del suavizado (``n=15``, ``\omega=\frac23``). El error inicial mezcla modos suaves y oscilatorios; tres barridos de Jacobi amortiguado eliminan casi por completo las componentes oscilatorias y dejan un error suave. Notar que el error *no* se hizo chico: solo cambió de forma. Los dos gráficos están en la misma escala vertical.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000003c
md"""
### 5.7.3 La idea de la corrección en malla gruesa

¿Cómo atacamos los modos suaves, que la relajación no reduce? La observación central de multigrid es la siguiente: *un modo que es suave en la malla fina se ve más oscilatorio en una malla más gruesa*.

Esto se ve de inmediato en la fórmula (5.17). Supongamos que ``n`` es impar y consideremos la *malla gruesa* formada por los nodos de índice par, ``x_{2j}`` con ``j=1,\dots,n_c``, donde ``n_c=\frac{n-1}2``: es una malla del mismo intervalo, con paso ``2h`` y la mitad de los puntos. Si muestreamos el modo ``\varphi^{(k)}`` de la malla fina en esos nodos, obtenemos el vector de componentes
```math
\varphi_{2j}^{(k)} = \text{sen}\left(\frac{k(2j)\pi}{n+1}\right) = \text{sen}\left(\frac{kj\pi}{n_c+1}\right), \qquad j=1,\dots,n_c,
```
donde usamos que ``n_c+1=\frac{n+1}2``. Pero esa es exactamente la fórmula (5.17) del modo número ``k`` de la malla gruesa. En otras palabras: el modo ``k`` de la malla fina se transforma en el modo ``k`` de la malla gruesa. Como la malla gruesa tiene la mitad de modos (``n_c\approx n/2`` en vez de ``n``), la frecuencia *relativa* del modo se duplicó: un ``k`` que era chico comparado con ``n`` puede ser grande comparado con ``n_c`` (Figura 5.5).

**Figura 5.5** (descripción, no reproducida gráficamente). El mismo modo, visto en dos mallas. Arriba, ``\varphi^{(4)}`` en la malla fina de 15 nodos: es un modo de frecuencia baja respecto de esa malla. Abajo, sus valores en los 7 nodos de la malla gruesa: ahora cambia de signo de un nodo al siguiente, o sea que es un modo *oscilatorio* respecto de la malla gruesa, y por lo tanto la relajación sí lo reduce.

Por lo tanto, si trasladamos el problema del error a una malla gruesa, allí ese error se vuelve accesible a la relajación; y además la malla gruesa tiene la mitad de incógnitas, con lo cual es más barata de tratar.

Falta un ingrediente: la incógnita del problema es ``u``, no el error. Pero el error satisface un sistema con la misma matriz. En efecto, si ``u`` es una aproximación cualquiera y ``e=u_*-u``, entonces
```math
A_he = A_hu_*-A_hu = f-A_hu =: r, \tag{5.22}
```
donde ``r`` es el residuo, que sí podemos calcular. La igualdad (5.22) se llama *ecuación residual*: es un sistema con la misma matriz ``A_h`` y con dato el residuo. Resolverla exactamente equivaldría a resolver el problema original; pero como después de suavizar sabemos que ``e`` es suave, podemos aproximarlo resolviendo una versión de (5.22) en la malla gruesa, que es un problema de la misma forma pero con la mitad de incógnitas. Una vez obtenida esa corrección aproximada, la llevamos de vuelta a la malla fina y corregimos ``u\leftarrow u+e``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000003d
md"""
### 5.7.4 Operadores de transferencia entre mallas

Para hacer precisa la idea anterior necesitamos operadores que trasladen vectores entre la malla fina ``\Omega_h`` (con ``n`` puntos interiores) y la malla gruesa ``\Omega_{2h}`` (con ``n_c=\frac{n-1}2`` puntos interiores y paso ``2h``). Para que el proceso pueda repetirse suponemos que ``n+1=2^L`` para algún ``L\in\mathbb{N}``; en particular ``n`` es impar y ``n_c`` es entero. El punto grueso ``j`` coincide con el punto fino ``2j``, es decir, ``x_j^{2h}=x_{2j}^h`` (Figura 5.6).

**Figura 5.6** (descripción, no reproducida gráficamente). Malla fina y malla gruesa para ``n=15`` (``n_c=7``). Los nodos de la malla gruesa son los de índice par de la fina; los de índice impar existen solo en la malla fina.

**Prolongación (interpolación).** El operador de prolongación ``I_{2h}^h:\mathbb{R}^{n_c}\to\mathbb{R}^n`` lleva un vector de la malla gruesa a la fina por interpolación lineal:
```math
(I_{2h}^hv)_{2j} = v_j, \qquad (I_{2h}^hv)_{2j+1} = \frac12(v_j+v_{j+1}),
```
con la convención ``v_0=v_{n_c+1}=0``. En los nodos que comparten ambas mallas se copia el valor; en los nodos intermedios se promedian los vecinos.

**Figura 5.7** (descripción, no reproducida gráficamente). El operador de prolongación. Arriba, un vector ``v\in\mathbb{R}^{n_c}`` de la malla gruesa. Abajo, el vector ``I_{2h}^hv\in\mathbb{R}^n``: en los nodos que ambas mallas comparten se copia el valor de ``v``, y en los nodos intermedios se pone el promedio de los dos vecinos. Los valores nuevos caen exactamente sobre los segmentos que unen los viejos: por eso decimos que ``I_{2h}^h`` interpola linealmente.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000003e
md"""
**Restricción.** El operador de restricción ``I_h^{2h}:\mathbb{R}^n\to\mathbb{R}^{n_c}`` lleva un vector de la malla fina a la gruesa. Usamos la *ponderación completa* (*full weighting*):
```math
(I_h^{2h}w)_j = \frac14(w_{2j-1}+2w_{2j}+w_{2j+1}), \qquad j=1,\dots,n_c.
```
Es decir, el valor grueso es un promedio pesado del valor fino correspondiente y de sus dos vecinos. (Podría pensarse en simplemente copiar, ``(I_h^{2h}w)_j=w_{2j}``; el promedio funciona mejor porque no ignora la información de los nodos impares, y porque hace valer la relación (5.23) de más abajo.)

El peso ``\frac14(1,2,1)`` se suele representar con el siguiente diagrama, llamado *estarcido* (*stencil*) del operador: en la malla fina, los pesos ``\frac14``, ``\frac12``, ``\frac14`` sobre ``w_{2j-1}``, ``w_{2j}``, ``w_{2j+1}`` respectivamente, confluyen en el nodo grueso ``(I_h^{2h}w)_j``.

Estos dos operadores no son independientes: se verifica directamente que
```math
I_h^{2h} = \frac12(I_{2h}^h)^T. \tag{5.23}
```
Esta relación *variacional* entre restricción y prolongación es la que garantiza que el operador de dos mallas resulte simétrico, algo esencial para poder usar multigrid como precondicionador del gradiente conjugado.

**Figura 5.8** (descripción, no reproducida gráficamente). El operador de restricción. Arriba, un vector ``w`` de la malla fina, formado por una parte suave más una oscilación de la forma ``(-1)^i``. Abajo, su restricción ``I_h^{2h}w``: el promedio pesado ``\frac14(1,2,1)`` cancela exactamente la componente oscilatoria y reproduce la parte suave (aquí, multiplicada por ``\cos^2(\pi h/2)\approx0.96``). Ambos gráficos están en la misma escala vertical. Si en lugar de promediar se copiara el valor, ``(I_h^{2h}w)_j=w_{2j}``, la oscilación *no* se cancelaría y aparecería en la malla gruesa disfrazada de componente suave.

**Operador grueso.** Finalmente, en la malla gruesa usamos la matriz de la misma familia (5.16), pero con ``n_c`` incógnitas y paso ``2h``:
```math
A_{2h} = \frac{1}{(2h)^2}\begin{pmatrix}2&-1&&\\-1&2&-1&\\&\ddots&\ddots&\ddots\\&&-1&2\end{pmatrix}\in\mathbb{R}^{n_c\times n_c}.
```
Puede probarse (ejercicio) que esta elección coincide con el *operador de Galerkin* ``A_{2h}=I_h^{2h}A_hI_{2h}^h``, otra manera estándar de construir el problema grueso.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000003f
md"""
### 5.7.5 El esquema de dos mallas y el ciclo en V

Reuniendo las piezas, una iteración del *esquema de dos mallas* para ``A_hu=f``, partiendo de una aproximación ``u``, consiste en:

1. *Pre-suavizado.* Aplicar ``\nu_1`` barridos de Jacobi amortiguado a ``A_hu=f``.
2. *Residuo.* Calcular ``r=f-A_hu``.
3. *Restricción.* Trasladar el residuo a la malla gruesa: ``r_{2h}=I_h^{2h}r``.
4. *Corrección gruesa.* Resolver ``A_{2h}e_{2h}=r_{2h}``.
5. *Prolongación y corrección.* Interpolar y corregir: ``u\leftarrow u+I_{2h}^he_{2h}``.
6. *Post-suavizado.* Aplicar ``\nu_2`` barridos de Jacobi amortiguado a ``A_hu=f``.

El pre-suavizado elimina los modos oscilatorios del error; la corrección gruesa se ocupa de los modos suaves; el post-suavizado limpia las nuevas componentes oscilatorias que introduce la interpolación.

Queda un punto: el paso (4) requiere resolver un sistema con la matriz ``A_{2h}`` de la malla gruesa. Pero *ese* sistema tiene exactamente la misma estructura que el original, solo que con la mitad de las incógnitas. Entonces, en lugar de resolverlo de forma exacta, podemos aplicarle *recursivamente* el mismo esquema de dos mallas, bajando a una malla aún más gruesa, y así sucesivamente, hasta llegar a una malla tan pequeña (por ejemplo, un solo punto interior) que el sistema se resuelva de manera directa. Esta recursión es el *ciclo en V* (o *V-cycle*), esquematizado en la Figura 5.9.

El costo de un ciclo en V es proporcional al trabajo en la malla fina más el de todas las mallas gruesas. Como cada nivel tiene aproximadamente la mitad de puntos que el anterior, el costo total es
```math
cn\left(1+\frac12+\frac14+\cdots\right) \leq 2cn = O(n),
```
es decir, *lineal* en el número de incógnitas. Además, el factor de reducción del error de un ciclo en V está acotado por una constante ``\rho<1`` *independiente de* ``h``: es el Teorema 5.34, que enunciamos más abajo. Combinando ambos hechos: para reducir el error a una tolerancia fija hace falta un número de ciclos independiente de ``n``, y cada ciclo cuesta ``O(n)``, de modo que multigrid resuelve el problema de Poisson con costo óptimo ``O(n)``, muy por debajo del ``O(n^2)`` del gradiente conjugado sin precondicionar.

**Figura 5.9** (descripción, no reproducida gráficamente). Un ciclo en V con cuatro niveles (``\Omega_h,\Omega_{2h},\Omega_{4h},\Omega_{8h}``). En la rama descendente se suaviza y se restringe el residuo hacia mallas cada vez más gruesas; en el nivel más grueso se resuelve en forma exacta; en la rama ascendente se prolonga la corrección y se vuelve a suavizar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000040
md"""
### 5.7.6 Implementación

A continuación damos una implementación del ciclo en V en MATLAB/OCTAVE. El suavizador es el Jacobi amortiguado (5.20) (matriz-libre, aprovechando la estructura tridiagonal); `restringir` y `prolongar` implementan ``I_h^{2h}`` e ``I_{2h}^h``, y `vciclo` es la función recursiva. Seguimos la convención de nombres en minúscula del resto del texto.

```octave
function u = jacobiamortiguado(u,f,h,nu)
% function u = jacobiamortiguado(u,f,h,nu)
% Aplica nu barridos de Jacobi amortiguado (omega=2/3) al
% sistema de Poisson 1D A_h u = f sobre la malla interior de paso h.
% u, f son vectores columna de largo n (nodos interiores).
n = length(u);
omega = 2/3;
for m = 1:nu
   % residuo r = f - A_h u, con A_h = tridiag(-1,2,-1)/h^2
   Au = (2*u - [u(2:n);0] - [0;u(1:n-1)]) / h^2;
   u = u + omega * (h^2/2) * (f - Au);
end

function rc = restringir(r)
% function rc = restringir(r)
% Ponderacion completa (full weighting) de la malla fina a la gruesa.
% r: vector de largo n (n impar) en la malla fina.
% rc: vector de largo nc=(n-1)/2 en la malla gruesa.
n = length(r);
nc = (n-1)/2;
rc = zeros(nc,1);
for j = 1:nc
   rc(j) = ( r(2*j-1) + 2*r(2*j) + r(2*j+1) ) / 4;
end

function v = prolongar(vc)
% function v = prolongar(vc)
% Interpolacion lineal de la malla gruesa a la fina.
% vc: vector de largo nc en la malla gruesa.
% v: vector de largo n=2*nc+1 en la malla fina.
nc = length(vc);
n = 2*nc+1;
v = zeros(n,1);
vext = [0; vc; 0];              % agrega los ceros del borde
for j = 1:nc
   v(2*j) = vc(j);              % nodos comunes: se copia
end
for j = 0:nc
   v(2*j+1) = ( vext(j+1) + vext(j+2) ) / 2;   % nodos intermedios: promedio
end

function u = vciclo(u,f,h,nu1,nu2)
% function u = vciclo(u,f,h,nu1,nu2)
% Un ciclo en V para el sistema de Poisson 1D A_h u = f.
% nu1, nu2: numero de barridos de pre- y post-suavizado.
n = length(u);
if n == 1
   u = f * h^2 / 2;             % malla mas gruesa: A_h es escalar 2/h^2
   return
end
u  = jacobiamortiguado(u,f,h,nu1);           % (1) pre-suavizado
Au = (2*u - [u(2:n);0] - [0;u(1:n-1)]) / h^2;
r  = f - Au;                                  % (2) residuo
rc = restringir(r);                           % (3) restriccion
ec = zeros(length(rc),1);
ec = vciclo(ec,rc,2*h,nu1,nu2);               % (4) correccion gruesa (recursiva)
u  = u + prolongar(ec);                       % (5) prolongacion y correccion
u  = jacobiamortiguado(u,f,h,nu2);            % (6) post-suavizado
```

Iterando la función `vciclo` hasta que el residuo sea suficientemente pequeño se obtiene un método iterativo (multigrid) que resuelve el sistema con costo óptimo.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000041
md"""
### 5.7.7 Multigrid como precondicionador del gradiente conjugado

Volvamos a nuestro tema. Un ciclo en V define, de manera natural, una acción lineal ``\xi\mapsto M\xi``: dado un vector ``\xi``, aplicamos *un* ciclo en V al sistema ``A_he=\xi`` partiendo de ``e=0``, y el resultado es ``M\xi``. Es decir, ``M`` es el operador "un ciclo en V con dato inicial nulo". Este ``M`` tiene justamente las dos propiedades que pedíamos de un buen precondicionador:

- La acción ``\xi\mapsto M\xi`` es *barata*: cuesta ``O(n)``, por lo visto arriba.
- ``M\approx A_h^{-1}`` en el sentido de la Sección 5.1.6, con ``\frac{c_2}{c_1}`` acotado *independientemente de* ``h``. Esto no es evidente, y lo deducimos a continuación del hecho de que un ciclo en V reduce el error en un factor ``\rho<1`` independiente de ``h``.

Para que ``M`` sea sdp (y así poder usar el gradiente conjugado precondicionado), basta con simetrizar el ciclo: tomar ``\nu_1=\nu_2`` y usar un suavizador simétrico. El Jacobi amortiguado ya es simétrico (pues ``D^{-1}`` lo es), y la relación variacional (5.23) entre restricción y prolongación garantiza que el operador de dos mallas resultante sea simétrico. Supondremos de aquí en más que se hicieron esas elecciones, de modo que ``M`` es sdp.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000042
md"""
**De la convergencia del ciclo en V a la cota del número de condición.** Veamos con cuidado cómo se encadenan los dos hechos, porque el argumento es instructivo y usa todo lo visto en las dos secciones anteriores. Consta de tres pasos.

*Paso 1: el ciclo en V es un método iterativo estacionario con matriz ``M``.* Un ciclo en V aplicado a ``A_hu=f`` partiendo de una aproximación ``u`` produce una nueva aproximación que depende de manera afín de ``u`` y de ``f``; además, es *consistente*, en el sentido de que si ``u`` ya es la solución exacta, el ciclo la deja fija (el residuo es nulo, y entonces todos los suavizados y correcciones no cambian nada). Por lo tanto el ciclo se escribe como
```math
u \longmapsto u+M(f-A_hu) \tag{5.24}
```
para cierta matriz ``M``, y tomando ``u=0`` y ``f=\xi`` vemos que ese ``M`` es exactamente el operador "un ciclo en V con dato inicial nulo" definido arriba. Es decir, multigrid es un método estacionario como los del Capítulo 4, cuya matriz de iteración es
```math
E := I-MA_h,
```
ya que si ``u_*`` es la solución exacta, restando ``u_*=u_*+M(f-A_hu_*)`` de (5.24) se obtiene que el nuevo error es ``E`` por el error viejo.

*Paso 2: la convergencia del ciclo en V acota el espectro de ``MA_h``.* Aquí usamos el resultado central de la teoría de multigrid, cuya demostración excede los alcances de este curso.

**Teorema 5.34** (Convergencia del ciclo en V). *Consideremos el ciclo en V para el sistema (5.16), con Jacobi amortiguado (``\omega=\frac23``) como suavizador, ``\nu_1=\nu_2=\nu\geq1`` barridos de pre- y post-suavizado, y los operadores de transferencia ``I_h^{2h}``, ``I_{2h}^h`` definidos arriba. Sea ``E=I-MA_h`` su matriz de iteración. Entonces existe una constante ``\rho<1``, independiente de ``h`` y de la cantidad de niveles, tal que*
```math
\|Ev\|_{A_h} \leq \rho\|v\|_{A_h}, \qquad \forall v\in\mathbb{R}^n. \tag{5.25}
```
*Además, ``\rho`` se puede hacer tan chico como se quiera aumentando el número ``\nu`` de barridos de suavizado.*

*Referencias.* Una presentación accesible, con el análisis de dos mallas para este mismo problema modelo, puede verse en [BHM]; el tratamiento completo (incluyendo la independencia respecto de la cantidad de niveles, que es la parte difícil, y la extensión a problemas más generales) está en [TOS].

Lo esencial de (5.25) es que ``\rho`` *no depende de* ``h``: por más que refinemos la malla, cada ciclo en V reduce el error en norma ``A_h`` por el mismo factor. Es el mismo tipo de enunciado que vimos para los métodos estacionarios del Capítulo 4, pero allí el factor análogo tendía a 1 al refinar.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000043
md"""
Observemos ahora que ``MA_h`` es autoadjunta respecto del producto escalar ``\langle u,v\rangle_{A_h}:=u^TA_hv``: si ``M`` y ``A_h`` son simétricas,
```math
\langle MA_hu,v\rangle_{A_h} = (MA_hu)^TA_hv = u^TA_hMA_hv = \langle u,MA_hv\rangle_{A_h}.
```
Lo mismo vale entonces para ``E=I-MA_h``. En consecuencia ``E`` es diagonalizable con autovalores reales, y su radio espectral está acotado por su norma: por (5.25),
```math
|\lambda| \leq \rho \qquad \text{para todo } \lambda\in\sigma(E)=\sigma(I-MA_h).
```
Como ``\sigma(I-MA_h)=\{1-\nu : \nu\in\sigma(MA_h)\}``, esto dice precisamente que
```math
\sigma(MA_h) \subset [1-\rho,\ 1+\rho]. \tag{5.26}
```
Es decir: *el sistema precondicionado tiene todos sus autovalores agrupados alrededor de 1*, que es exactamente la situación ideal del Ejemplo 5.17 y la que persigue todo precondicionador.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000044
md"""
*Paso 3: de ``\sigma(MA_h)`` a ``\frac{c_2}{c_1}`` y a ``\kappa_2(SA_hS)``.* Sea ``S=M^{1/2}``. Por un lado, ``SA_hS`` y ``MA_h`` son semejantes, ya que ``SA_hS=S^{-1}(S^2A_h)S=S^{-1}(MA_h)S``, y por lo tanto tienen el mismo espectro. Por otro lado, las mejores constantes ``c_1,c_2`` de la definición de precondicionador son también autovalores de ``MA_h``: escribiendo ``x=A_h^{1/2}z`` resulta
```math
\frac{x^TMx}{x^TA_h^{-1}x} = \frac{z^TA_h^{1/2}MA_h^{1/2}z}{z^Tz},
```
de modo que el cociente recorre exactamente el intervalo ``[\lambda_{\min}(B),\lambda_{\max}(B)]`` con ``B:=A_h^{1/2}MA_h^{1/2}``, y ``B`` es semejante a ``MA_h`` (pues ``B=A_h^{1/2}(MA_h)A_h^{-1/2}``). Juntando todo,
```math
c_1 = \lambda_{\min}(MA_h), \qquad c_2 = \lambda_{\max}(MA_h), \qquad \frac{c_2}{c_1} = \kappa_2(SA_hS),
```
y por (5.26),
```math
\kappa_2(SA_hS) = \frac{c_2}{c_1} \leq \frac{1+\rho}{1-\rho}, \tag{5.27}
```
que *no depende de* ``h``, porque ``\rho`` no depende de ``h``.

**Qué significa en la práctica.** Reemplazando (5.27) en la cota de convergencia del gradiente conjugado precondicionado, el error en norma ``A_h`` se reduce en cada iteración por el factor
```math
\frac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1} \leq \frac{\sqrt{\frac{1+\rho}{1-\rho}}-1}{\sqrt{\frac{1+\rho}{1-\rho}}+1}.
```
Por ejemplo, con un valor típico ``\rho=0.1`` se obtiene ``\frac{c_2}{c_1}\leq1.22`` y un factor de reducción menor que ``0.05`` por iteración: unas 5 iteraciones alcanzan para ganar seis cifras, *cualquiera sea* ``n``. Comparar con el gradiente conjugado sin precondicionar, que para esta matriz necesita ``O(n)`` iteraciones.

A diferencia de los precondicionadores genéricos, cuya mejora sobre el sistema sin precondicionar es modesta, el precondicionador multigrid logra que el número de condición del sistema precondicionado quede *acotado independientemente de* ``h``. En consecuencia, el gradiente conjugado precondicionado con multigrid converge en un número de iteraciones que *no crece* al refinar la malla: es el escenario ideal anticipado al comienzo de la sección, y se logra precisamente por haber usado la información sobre *de dónde proviene* la matriz ``A_h``.

**Observación 5.35.** Aunque presentamos multigrid sobre un problema modelo unidimensional por simplicidad, las mismas ideas (suavizado, corrección en malla gruesa, transferencia entre mallas y recursión) se extienden a dos y tres dimensiones y a problemas más generales, y es en esos casos donde multigrid muestra toda su potencia como resolvedor y como precondicionador de costo óptimo.
"""

# ╔═╡ e5f6a7b8-dddd-4a1b-8c2d-000000000004
md"""
**Nota** (de las diapositivas de la Clase 12, resultados numéricos para el problema modelo (5.16)). La siguiente tabla muestra la cantidad de iteraciones necesarias para alcanzar ``\|r_k\|_2/\|b\|_2<10^{-10}``, según el tamaño ``n`` y el método:

| ``n`` | GC | GCP-Jacobi | Multigrid solo | GCP-MG (``\nu=1``) | GCP-MG (``\nu=2``) |
|---|---|---|---|---|---|
| 31 | 16 | 16 | 15 | 8 | 6 |
| 127 | 64 | 64 | 16 | 8 | 6 |
| 511 | 256 | 256 | 16 | 8 | 6 |
| 2047 | 1024 | 1024 | 16 | 17 | 6 |
| 4095 | 2049 | 2049 | 16 | 8 | 5 |

Algunas observaciones sobre estos números: el GC sin precondicionar duplica sus iteraciones al duplicar ``n``, es decir, crecen como ``O(n)``, tal como predice la Observación 5.20. La columna GCP-Jacobi coincide exactamente con la del GC, porque para esta matriz ``D=\frac2{h^2}I`` es múltiplo de la identidad, y por lo tanto el precondicionador diagonal no hace nada. En cambio, tanto multigrid solo como GCP-MG mantienen una cantidad de iteraciones que no crece con ``n``, de acuerdo con el Teorema 5.34, y GCP-MG necesita aproximadamente la mitad de ciclos que multigrid solo: el gradiente conjugado actúa como acelerador de multigrid. Esta tabla ilustra numéricamente lo pedido en el Problema 5.14; se sugiere reproducirla (o una tabla análoga) antes de mirarla.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000045
md"""
## 5.8 Gradiente conjugado para problemas no simétricos

Supongamos que tenemos una matriz invertible ``A``, no necesariamente simétrica, y consideramos el sistema
```math
Ax=b.
```

**CGNR.**

El sistema ``Ax=b`` es equivalente al sistema normal dado por
```math
A^TAx = A^Tb,
```
donde ahora sí, la matriz del sistema, es decir, ``A^TA``, es sdp. Por lo tanto, podemos aplicar el método del gradiente conjugado a este sistema, observando que no es necesario calcular el producto ``A^TA``, sino que es suficiente proveer la rutina ``\xi\to A^T(A\xi)``. Por otro lado, es importante mencionar que este método puede no resultar conveniente ya que
```math
\kappa_2(A^TA) \approx \kappa_2(A)^2,
```
y así, el número de condición del sistema normal empeora si ``\kappa_2(A)`` es grande.

Puesto que GC aplicado al sistema normal minimiza ``\|\cdot-x_*\|_{A^TA}`` sobre ``x_0+\mathcal{K}_k`` en el nivel ``k`` y puesto que
```math
\|x-x_*\|_{A^TA}^2 = (x-x_*)^TA^TA(x-x_*) = (A(x-x_*))^TA(x-x_*) = \|A(x-x_*)\|_2^2 = \|Ax-b\|_2^2,
```
este método se conoce como CGNR (*Conjugate Gradient on the Normal equations to minimize the Residual*).

**CGNE.**

El sistema ``Ax=b`` es equivalente al siguiente:
```math
\begin{cases}AA^Ty=b\\x=A^Ty.\end{cases}
```
En este caso, podemos aplicar GC para resolver el sistema ``AA^Ty=b``, puesto que ``AA^T`` es sdp. Notemos que en cada paso ``k``, se minimiza ``\|\cdot-y_*\|_{AA^T}`` sobre ``y_0+\mathcal{K}_k``. Tomando ``x=A^Ty``, se tiene que
```math
\|y-y_*\|_{AA^T}^2 = (y-y_*)^TAA^T(y-y_*) = (A^T(y-y_*))^TA^T(y-y_*) = \|A^T(y-y_*)\|_2^2 = \|x-x_*\|_2^2,
```
y por esto, este método se conoce como CGNE (*Conjugate Gradient on the Normal equations to minimize the Error*).
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000046
md"""
## 5.9 GMRes

### 5.9.1 La propiedad de minimización y sus consecuencias

El método GMRes (*Generalized Minimum Residual*) es un método de subespacio de Krylov para matrices que no son sdp. A diferencia de CGNR y CGNE, la iteración de GMRes no requiere el cálculo de la acción de ``A^T`` sobre un vector, lo cual puede ser ventajoso en algunos casos. Sin embargo, requiere el almacenamiento de una base para ``\mathcal{K}_k``, y por lo tanto el almacenamiento necesario es mayor.

La ``k``-ésima iteración de GMRes es la solución al problema de mínimos cuadrados siguiente
```math
x_k = \operatorname*{argmin}_{x\in x_0+\mathcal{K}_k}\|b-Ax\|_2, \tag{5.28}
```
donde ``\mathcal{K}_k:=\text{span}\{r_0,Ar_0,\dots,A^{k-1}r_0\}``

Dejamos para más adelante la implementación efectiva (y eficiente) de cada iteración y nos enfocamos a continuación en las consecuencias que se desprenden de esta definición.

El análisis que sigue es parecido al que se hizo para GC. Observemos que
```math
x\in x_0+\mathcal{K}_k \iff x=x_0+\sum_{j=0}^{k-1}\gamma_jA^jr_0,
```
y luego
```math
r = b-Ax = b-Ax_0-\sum_{j=0}^{k-1}\gamma_jA^{j+1}r_0 = r_0-\sum_{j=1}^{k}\gamma_{j-1}A^jr_0.
```
Es decir, ``x\in x_0+\mathcal{K}_k`` si y solo si ``r=b-Ax=p(A)r_0``, con ``p\in\overline{\mathcal{P}}_k``, un polinomio *residual* (``p(0)=1``) de grado ``\leq k``.

Queda entonces demostrado el siguiente teorema:
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000047
md"""
**Teorema 5.36.** *Sea ``A`` una matriz no singular y sea ``x_k`` la ``k``-ésima iteración de GMRes. Entonces*
```math
\|r_k\|_2 = \min_{p\in\overline{\mathcal{P}}_k}\|p(A)r_0\|_2,
```
*que implica*
```math
\|r_k\|_2 \leq \|p(A)r_0\|_2 \leq \|p(A)\|_2\|r_0\|_2, \tag{5.29}
```
*para todo* ``p\in\overline{\mathcal{P}}_k``.

Inmediatamente tenemos el siguiente corolario.

**Corolario 5.37.** *Sea ``A`` una matriz no singular y sea ``x_k`` la ``k``-ésima iteración de GMRes. Entonces, cualquiera sea ``p\in\overline{\mathcal{P}}_k``*
```math
\frac{\|r_k\|_2}{\|r_0\|_2} \leq \|p(A)\|_2.
```

Podemos ahora aplicar el corolario para demostrar que GMRes termina en un número finito de pasos.

**Teorema 5.38.** *Sea ``A\in\mathbb{R}^{n\times n}`` una matriz no singular, entonces GMRes encontrará la solución, a lo sumo, en ``n`` pasos.*

**Demostración.** El polinomio característico de ``A`` es ``p(z)=\det(A-zI)``, y es de grado ``n``. Como ``A`` es no singular, resulta ``p(0)=\det(A)\neq0`` y entonces consideramos el polinomio residual
```math
\bar p(z) := \frac{p(z)}{p(0)} \in \overline{\mathcal{P}}_n.
```
El Teorema de Cayley-Hamilton dice que el polinomio característico de una matriz ``A`` aplicado a la misma matriz ``A`` da la matriz nula. En este caso, ``p(A)=0``, y por lo tanto ``\bar p(A)=0``, que a su vez, por el Teorema 5.36, implica
```math
r_n = b-Ax_n = 0,
```
y ``x_n`` es la solución.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000048
md"""
En la sección 5.1.2 hemos usado que para ``A`` simétrica hay una base ortonormal de autovectores, es decir ``A=V\Lambda V^T``, con ``V`` ortogonal y ``\Lambda`` diagonal, conteniendo los autovalores de ``A`` en la diagonal. Esto sirvió para construir polinomios que permitieron obtener cotas para el decrecimiento del error. No podemos hacer lo mismo si ``A`` no es simétrica. Sin embargo, si suponemos ``A`` diagonalizable, podemos usar (5.29) e información sobre los autovalores de ``A`` para obtener cotas de decrecimiento similares a las obtenidas para GC. Supongamos que ``A`` es diagonalizable, y entonces existe una matriz no singular (posiblemente compleja) ``V`` tal que
```math
A = V\Lambda V^{-1}.
```
Aquí ``\Lambda`` es una matriz diagonal (posiblemente compleja) con los autovalores de ``A`` sobre la diagonal, y las columnas de ``V`` son los autovectores correspondientes. Además, ``A^k=V\Lambda^kV^{-1}``. Si ``A`` es normal¹ las columnas de ``V`` son ortonormales y ``V^{-1}=\bar V^T=V^H``.

*(¹Se dice que una matriz ``A`` es normal si ``A`` conmuta con ``A^T``, i.e., ``AA^T=A^TA``.)*

Se puede usar la estructura de matrices diagonalizables para probar el siguiente resultado, análogo al Corolario 5.14 para GC.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000049
md"""
**Teorema 5.39.** *Sea ``A=V\Lambda V^{-1}`` una matriz no singular diagonalizable, con ``\Lambda`` diagonal y sea ``x_k`` la ``k``-ésima iteración de GMRes. Entonces, para todo ``p\in\overline{\mathcal{P}}_k``, resulta*
```math
\frac{\|r_k\|_2}{\|r_0\|_2} \leq \kappa_2(V)\max_{z\in\sigma(A)}|p(z)|. \tag{5.30}
```

**Demostración.** Sea ``p\in\overline{\mathcal{P}}_k``. Entonces ``p(A)=Vp(\Lambda)V^{-1}`` y luego
```math
\|p(A)\|_2 \leq \|V\|_2\|V^{-1}\|_2\|p(\Lambda)\|_2 = \kappa_2(V)\max_{z\in\sigma(A)}|p(z)|.
```

Observemos que si ``A`` es normal, entonces ``V`` es ortogonal y ``\|V\|_2=\|V^{-1}\|_2=1``, por lo que ``\kappa_2(V)=1``. En general es difícil calcular ``\kappa_2(V)``, y se espera que sea grande si ``V`` está cerca de ser no diagonalizable.

Enunciamos ahora un par de consecuencias sencillas de (5.30).

**Corolario 5.40.** *Sea ``A`` no singular y diagonalizable. Si ``A`` tiene solo ``k`` autovalores distintos, entonces GMRes terminará en a lo sumo ``k`` iteraciones.*

**Corolario 5.41.** *Sea ``A`` una matriz no singular. Si ``b`` es una combinación lineal de ``k`` autovectores de ``A``, entonces el método GMRes aplicado a la resolución de ``Ax=b``, partiendo de ``x_0=0`` terminará en a lo sumo ``k`` iteraciones.*

**Criterio de parada.**

Como ocurre con GC, el método GMRes se considera usualmente un método iterativo. La tasa de convergencia para matrices diagonalizables involucra ``\kappa_2(V)``, pero igualmente se asemeja a GC. También se han hallado tasas de convergencia para matrices no diagonalizables, pero dado que el conjunto de tales matrices tiene medida cero en el espacio de las matrices complejas de ``n\times n``, no las analizaremos en este curso.

Como ocurre con GC terminaremos las iteraciones de GMRes cuando
```math
\|r_k\|_2 \leq \eta\|b\|_2.
```
Según el Lema 4.10, ``\dfrac{\|e_k\|}{\|x_*\|}\leq\kappa(A)\dfrac{\|r_k\|}{\|b\|}``. Por lo tanto, el criterio de parada garantizará que
```math
\frac{\|e_k\|_2}{\|x_*\|_2} \leq \kappa_2(A)\eta,
```
es decir, el error relativo estará acotado por la tolerancia ``\eta`` multiplicada por el número de condición de ``A``.
"""

# ╔═╡ e5f6a7b8-0000-4a1b-8c2d-000000000001
md"""
**Nota** (de las diapositivas de la Clase 13). Puede ser útil comparar, lado a lado, el gradiente conjugado y GMRes:

| | GC | GMRes |
|---|---|---|
| Hipótesis | ``A`` sdp | ``A`` no singular |
| Espacio | ``x_0+\mathcal{K}_k`` | ``x_0+\mathcal{K}_k`` |
| Minimiza | ``\|x-x_*\|_A`` (error en norma ``A``) | ``\|b-Ax\|_2`` (residuo en norma 2) |
| Polinomios | ``\|x_k-x_*\|_A=\min_p\|p(A)(x_0-x_*)\|_A`` | ``\|r_k\|_2=\min_p\|p(A)r_0\|_2`` |
| Terminación | ``\leq n`` pasos | ``\leq n`` pasos |
| Recurrencia | dos términos | no, hay que guardar toda la base |
| Memoria por iteración | ``O(n)``, fija | ``O(nk)``, crece con ``k`` |

Lo que se pierde sin simetría: para ``A`` sdp, ``A=V\Lambda V^T`` con ``V`` ortogonal, y eso daba ``\|p(A)\|_2=\max_{\lambda\in\sigma(A)}|p(\lambda)|`` (Corolario 5.14 vía la demostración del Corolario 5.14). Sin simetría esa igualdad ya no vale, y por eso aparece el factor ``\kappa_2(V)`` en el Teorema 5.39.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000004a
md"""
### 5.9.2 Un resultado de convergencia para GMRes

El siguiente teorema simple nos dice que si ``A`` *está cerca* de la matriz identidad ``I``, el residuo tiende a cero de manera lineal.

**Teorema 5.42.** *Sea ``A`` una matriz tal que ``\|I-A\|_2\leq\rho<1``. Entonces*
```math
\|r_k\|_2 \leq \rho^k\|r_0\|_2, \qquad k=1,2,\dots
```

**Demostración.** Consideremos ``p_k(z)=(1-z)^k``, entonces ``p_k\in\overline{\mathcal{P}}_k`` y ``p_k(A)=(I-A)^k``, por lo tanto, el Corolario 5.37 implica que
```math
\|r_k\|_2 \leq \|(I-A)^k\|_2\|r_0\|_2 \leq \|I-A\|_2^k\|r_0\|_2 \leq \rho^k\|r_0\|_2.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000004b
md"""
En el Lema 4.10 vimos que
```math
\frac{1}{\kappa(A)}\frac{\|A\tilde x-Ax\|}{\|Ax\|} \leq \frac{\|\tilde x-x\|}{\|x\|} \leq \kappa(A)\frac{\|A\tilde x-Ax\|}{\|Ax\|}.
```
La base para demostrar esto era que, en general, para ``S\in\mathbb{R}^{n\times n}`` invertible,
```math
\|y\| \leq \|S^{-1}Sy\| \leq \|S^{-1}\|\|Sy\|, \qquad\text{y}\qquad \|Sy\|\leq\|S\|\|y\|, \qquad\forall y\in\mathbb{R}^n.
```
En consecuencia,
```math
\frac{\|Sy\|}{\|S\|} \leq \|y\| \leq \|S^{-1}\|\|Sy\|, \qquad\text{y}\qquad \frac{1}{\|S^{-1}\|\|S\tilde y\|} \leq \frac{1}{\|\tilde y\|} \leq \frac{\|S\|}{\|S\tilde y\|}, \qquad\forall y,\tilde y\in\mathbb{R}^n.
```
Finalmente
```math
\frac{1}{\|S\|\|S^{-1}\|}\frac{\|Sy\|}{\|S\tilde y\|} \leq \frac{\|y\|}{\|\tilde y\|} \leq \|S\|\|S^{-1}\|\frac{\|Sy\|}{\|S\tilde y\|}, \qquad\forall y,\tilde y\in\mathbb{R}^n,
```
o lo que es lo mismo
```math
\frac{1}{\kappa(S)}\frac{\|Sy\|}{\|S\tilde y\|} \leq \frac{\|y\|}{\|\tilde y\|} \leq \kappa(S)\frac{\|Sy\|}{\|S\tilde y\|}, \qquad\forall y,\tilde y\in\mathbb{R}^n. \tag{5.31}
```
Aplicando esta última desigualdad para ``S=A``, ``y=e_k:=x_*-x_k``, ``\tilde y=e_0:=x_*-x_0``, obtenemos
```math
\frac{\|e_k\|_2}{\|e_0\|} \leq \kappa_2(A)\frac{\|Ae_k\|_2}{\|Ae_0\|_2} = \kappa_2(A)\frac{\|r_k\|_2}{\|r_0\|_2} \leq \kappa_2(A)\rho^k.
```
En síntesis, si ``\|I-A\|_2=\rho<1``, el método GMRes produce soluciones aproximadas que satisfacen
```math
\|e_k\|_2 \leq \kappa_2(A)\|e_0\|\rho^k, \qquad k=1,2,\dots
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000004c
md"""
### 5.9.3 Precondicionamiento para GMRes

El precondicionamiento para GMRes es diferente del que se hace para GC. No debemos ahora preocuparnos porque el sistema precondicionado sea sdp y todo se reduce a hallar una matriz ``M`` que sea aproximada a ``A^{-1}`` para que se logre ``\|I-MA\|_2=:\rho_I<1`` o ``\|I-AM\|_2=:\rho_D<1``.

**Precondicionador a izquierda.**

Si disponemos de una matriz ``M`` que sea aproximada a ``A^{-1}`` de manera que
```math
\|I-MA\|_2 = \rho_I < 1,
```
podemos utilizar GMRes para resolver el sistema equivalente
```math
MAx = Mb,
```
y las iteraciones de GMRes satisfarán
```math
\frac{\|e_k\|_2}{\|e_0\|_2} \leq \kappa_2(MA)\frac{\|r_k\|_2}{\|r_0\|_2} \leq \kappa_2(MA)\rho_I^k,
```
aquí ``r_k=Mb-MAx_k``. Si además ``\kappa_2(MA)<<\kappa_2(A)``, el residuo relativo es más aproximado al error relativo y el criterio de parada será más preciso.

**Precondicionador a derecha.**

Si disponemos de una matriz ``M`` que sea aproximada a ``A^{-1}`` de manera que
```math
\|I-AM\|_2 = \rho_D < 1,
```
entonces consideramos el siguiente problema equivalente a ``Ax=b``:
```math
AMy = b, \qquad x = My.
```
La idea es entonces resolver ``(AM)y=b`` utilizando GMRes y luego obtener ``x`` multiplicando ``M`` por ``y``. En este caso, el residuo satisface
```math
\frac{\|r_k\|_2}{\|r_0\|_2} \leq \rho_D^k,
```
donde ``r_k=b-AMy_k``, que coincide con ``b-Ax_k`` y luego el error relativo cumple
```math
\frac{\|e_k\|_2}{\|e_0\|_2} \leq \kappa_2(A)\frac{\|r_k\|_2}{\|r_0\|_2} \leq \kappa_2(A)\rho_D^k.
```
Si bien se gana en velocidad de convergencia, puede no ganarse nada en cómo se relaciona el residuo con el error.

Un aspecto importante de GMRes precondicionado es que no hace falta modificar el algoritmo, basta con aplicar el algoritmo al sistema ``MAx=Mb`` o a ``AMy=b`` y luego calcular ``x=My``.
"""

# ╔═╡ e5f6a7b8-0000-4a1b-8c2d-000000000002
md"""
**Nota** (de las diapositivas de la Clase 13). Comparación lado a lado de las dos formas de precondicionar GMRes:

| | A izquierda: ``MAx=Mb`` | A derecha: ``AMy=b``, ``x=My`` |
|---|---|---|
| Condición | ``\|I-MA\|_2=\rho_I<1`` | ``\|I-AM\|_2=\rho_D<1`` |
| Qué minimiza GMRes | ``\|M(b-Ax)\|_2`` (residuo precondicionado) | ``\|b-Ax\|_2`` (residuo verdadero) |
| Cota del residuo | ``\leq\rho_I^k`` | ``\leq\rho_D^k`` |
| Cota del error relativo | ``\leq\kappa_2(MA)\rho_I^k`` | ``\leq\kappa_2(A)\rho_D^k`` |
| Paso final | ninguno | ``x=My`` |

La ventaja de precondicionar a izquierda es que, si además ``\kappa_2(MA)\ll\kappa_2(A)``, el residuo que ve el criterio de parada se parece más al error relativo verdadero. La ventaja de precondicionar a derecha es que el residuo que GMRes minimiza es el residuo verdadero ``b-Ax``, aunque la cota de error siga teniendo la constante ``\kappa_2(A)`` original. En ambos casos no hace falta modificar el algoritmo: basta con darle la rutina ``\xi\mapsto M(A\xi)`` o ``\xi\mapsto A(M\xi)``, según corresponda (a diferencia del GCP, donde hubo que reescribir el algoritmo para que ``S=M^{1/2}`` desapareciera).
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000004d
md"""
### 5.9.4 Implementación de GMRes

Comenzamos recordando que ``x_k`` resuelve un problema de mínimos cuadrados:
```math
x_k = \operatorname*{argmin}_{x\in x_0+\mathcal{K}_k}\|b-Ax\|_2.
```
Supongamos que ``\{v_1,v_2,\dots,v_k\}`` es una base ortonormal de ``\mathcal{K}_k``, por lo que
```math
z\in\mathcal{K}_k \iff z=\sum_{\ell=1}^ky_\ell v_\ell = \underbrace{[v_1\ v_2\ \dots\ v_k]}_{V_k}\begin{bmatrix}y_1\\y_2\\\vdots\\y_k\end{bmatrix},
```
con ``V_k\in\mathbb{R}^{n\times k}``, ``y\in\mathbb{R}^k``. El problema resulta entonces equivalente a hallar ``y\in\mathbb{R}^k`` que minimice
```math
\|b-A(x_0+V_ky)\|_2 = \|r_0-AV_ky\|_2.
```
Es decir que queremos
```math
\text{minimizar } \|r_0-AV_ky\|_2, \qquad\text{para } y\in\mathbb{R}^k. \tag{5.32}
```
Este es un problema estándar de mínimos cuadrados que podría resolverse, por ejemplo, con una factorización QR de ``AV_k``. El problema es que sería necesario calcular el producto de ``A`` por ``V_k`` en cada iteración del algoritmo.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000004e
md"""
Si utilizamos ortogonalización de Gram-Schmidt, podemos representar la solución de (5.32) eficientemente. La ortogonalización de Gram-Schmidt aplicada al espacio de Krylov ``\mathcal{K}_k`` se llama *proceso de Arnoldi*. Dados ``x_k``, ``b``, un mapeo ``\xi\to A\xi``, y una dimensión ``k``, el proceso de Arnoldi resulta:

**Algoritmo de Arnoldi.**

```octave
Dada A en R^{nxn}, x0, b en R^n y k en N.
r0 = b - A*x0;
v1 = r0 / ||r0||_2;
Para i = 1, 2, ..., k-1, hacer
    w{i+1} = A * vi;
    z{i+1} = w{i+1} - sum_{j=1}^{i} (w{i+1}' * vj) * vj;
    v{i+1} = z{i+1} / ||z{i+1}||_2;
Fin Para
```
"""

# ╔═╡ e73a4fc2-f9fd-4fbd-ab6c-4b98464e1b6f
"""
	proceso_de_arnoldi()

Ingenuo.
"""
function proceso_de_arnoldi(A, b, x0, atol=1e-6)
	nrow, ncol = size(A)
	
	r = Vector{Float64}(undef, nrow)
	w = similar(r)
	v = Vector{Vector{Float64}}(undef, nrow)
	
	r .= b - A * x0
	v[1] = normalize(r)

	@inbounds @simd for in in 1:nrow-1
		w .= A * v[i]
		r .= w[i] - sum(w' * v[i]) * v[i]
		v[i+1] = normalize(r)
	end

	return v
end

# ╔═╡ 522b9b6b-6401-4b75-844b-9c8db64a0ffd
md"""
El vector ``w_{i+1}=Av_i`` pertenece a ``\mathcal{K}_{i+1}``. Si en algún paso ``w_{i+1}\in\mathcal{K}_i``, la dimensión del espacio no aumenta, resulta ``z_{i+1}=0`` y el algoritmo falla (ocurre un *breakdown*). El siguiente lema dice que en este caso, resulta que la solución del problema pertenece a ``x_0+\mathcal{K}_i``, y por lo tanto ya no hace falta buscar más vectores ``v_i``.

**Lema 5.43.** *Sea ``A`` no singular y sean ``v_j``, ``j=1,2,\dots``, los vectores generados por el algoritmo de Arnoldi. Sea ``i`` el menor entero para el cual*
```math
z_{i+1} = Av_i-\sum_{j=1}^i(Av_i^Tv_j)v_j = 0.
```
*Entonces* ``x=A^{-1}b\in x_0+\mathcal{K}_i``.

**Demostración.** Por la construcción del algoritmo, ``Av_j\in\mathcal{K}_{j+1}`` para ``j=1,2,\dots,i``, y por hipótesis ``\{v_1,v_2,\dots,v_i\}`` es una base ortonormal de ``\mathcal{K}_i``. Como ``z_{i+1}=0`` resulta que ``Av_i\in\mathcal{K}_i``, y luego ``A\mathcal{K}_i\subset\mathcal{K}_i``. Si armamos la matriz ``V_i=[v_1\ v_2\ \dots\ v_i]``, resulta que
```math
AV_i = V_iH,
```
con ``H`` una matriz de ``i\times i``, que es no singular pues ``A`` es invertible.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000004f
md"""
Tomando ``\beta=\|r_0\|_2`` y ``e_1=(1,0,\dots,0)^T\in\mathbb{R}^i`` tenemos
```math
\|r_i\|_2 = \|b-Ax_i\|_2 = \|r_0-A(x_i-x_0)\|_2.
```
Como ``x_i-x_0\in\mathcal{K}_i`` existe ``y\in\mathbb{R}^i`` tal que ``x_i-x_0=V_iy``, y como ``r_0=\beta V_ie_1``,
```math
\|r_i\|_2 = \|\beta V_ie_1-AV_iy\|_2 = \|V_i\beta e_1-V_iHy\|_2 = \|V_i(\beta e_1-Hy)\|_2 = \|\beta e_1-Hy\|_2,
```
la última igualdad es consecuencia de que ``V_i`` es ortogonal. La elección ``y=\beta H^{-1}e_1`` es la que da ``\|r_i\|_2=0`` que es el mínimo residuo sobre ``x_0+\mathcal{K}_i``.

Si el proceso de Arnoldi no falla, podemos usarlo para implementar GMRes de manera eficiente. Llamemos
```math
h_{ji} := (Av_i)^Tv_j, \qquad 1\leq i\leq k,\ 1\leq j\leq k+1.
```
Por la construcción de Arnoldi, la matriz ``H\in\mathbb{R}^{(k+1)\times k}`` cuyas entradas son las ``h_{ji}`` recién definidas resulta una matriz Hessenberg superior². En efecto, por la construcción, ``v_{i+1}\perp v_j`` para ``j=1,2,\dots,i``, y ``Av_j\in\text{span}\{v_1,\dots,v_{j+1}\}``, luego ``v_{i+1}\perp Av_j`` para ``j=1,2,\dots,i-1``.

*(²Una matriz ``H`` se dice Hessenberg superior si ``h_{ij}=0``, para todo ``i>j+1``.)*
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000050
md"""
Además, a menos que el proceso termine prematuramente, construye matrices ``V_k`` con columnas ortonormales, tal que
```math
\underbrace{A}_{n\times n}\underbrace{V_k}_{n\times k} = \underbrace{V_{k+1}}_{n\times(k+1)}\underbrace{H_k}_{(k+1)\times k}.
```
Luego, para algún ``y_k\in\mathbb{R}^k``,
```math
r_k = b-Ax_k = r_0-A(x_k-x_0) = V_{k+1}\underbrace{(\beta e_1-H_ky_k)}_{\in\mathbb{R}^{k+1}}.
```
Es decir, ``x_k=x_0+V_ky_k``, con ``y_k=\operatorname*{argmin}_{y\in\mathbb{R}^k}\|\beta e_1-H_ky\|_2``. Más aún, ``r_k=b-Ax_k`` puede calcularse sin necesidad de conocer ``x_k``. En efecto, por la ortogonalidad de ``V_{k+1}`` resulta
```math
\|r_k\|_2 = \|V_{k+1}(\beta e_1-H_ky_k)\|_2 = \|\beta e_1-H_ky_k\|_2.
```
Recordemos que el objetivo es hallar, dado ``\varepsilon>0``, un vector ``x\in\mathbb{R}^n`` tal que
```math
\|b-Ax\|_2 \leq \varepsilon\|b\|_2.
```
La entrada del algoritmo es una iteración inicial ``x_0\in\mathbb{R}^n``, el vector del lado derecho ``b\in\mathbb{R}^n``, un mapeo ``\xi\to A\xi`` que representamos con ``A`` y la tolerancia ``\varepsilon``. Limitamos además el máximo de iteraciones por ``k_{\max}``. La salida del algoritmo será el vector ``x`` que satisfaga ``\|b-Ax\|_2\leq\varepsilon\|b\|_2`` o bien que se obtenga después de ``k_{\max}`` iteraciones.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000051
md"""
**Método GMRes.** (Versión preliminar 1)

```octave
Dada A en R^{nxn}, x0, b en R^n, epsilon > 0 y kmax en N.
r0 = b - A*x0;
beta = rho = ||r0||_2;
v1 = r0 / rho;
k = 0;
Mientras rho > epsilon*||b||_2 y k < kmax, hacer
    k = k + 1;
    w = A * vk;
    Para j = 1, ..., k, hacer
        h_{j,k} = w' * vj;      % = (A vk)' vj
    Fin Para
    v{k+1} = w - sum_{j=1}^{k} h_{j,k} * vj;      % (*!)
    h_{k+1,k} = ||v{k+1}||_2;
    v{k+1} = v{k+1} / h_{k+1,k};
    yk = argmin_{y en R^k} ||beta*e1 - Hk*y||_2;      % (*)
    rho = ||beta*e1 - Hk*yk||_2;      % (**)
Fin Mientras
xk = x0 + Vk*yk;
```
"""

# ╔═╡ bcbbe045-49c7-437e-a612-5438943056fc
md"""
**Observación 5.44.** Vale la pena hacer algunas observaciones:

- El vector ``x_k`` se calcula solo al final, y no se necesita durante las iteraciones.
- Este método requiere que se almacene la base ortonormal del espacio de Krylov a medida que se ejecuta la iteración. Esto quiere decir que para realizar ``k`` iteraciones de GMRes se deben almacenar ``k`` vectores de longitud ``n``. Para problemas con ``n`` grande, esto puede ser un problema. Una manera de implementar esto es restringir ``k_{\max}`` a lo que permita la memoria disponible, ejecutar GMRes y chequear el residuo ``b-Ax_k`` al finalizar, si ``k=k_{\max}`` (si ``k<k_{\max}`` entonces ya se logró el objetivo). Si la norma del residuo es mayor que ``\varepsilon\|b\|_2``, entonces se vuelve a llamar a GMRes con ``x_0=x_{k_{\max}}``. No hay una teoría de convergencia para el algoritmo con *restart*, y en efecto, puede disminuir la velocidad de convergencia. Sin embargo, cuando funciona (muy a menudo), reduce significativamente el costo de almacenamiento de la iteración.

- Este algoritmo puede implementarse fácilmente en MATLAB/OCTAVE y el paso (*) puede resolverse directamente con la barra invertida "``\backslash``" en ``O(k^3)`` operaciones. Esto está bien cuando ``k<<n``, que es lo más usual. Otra opción, más recomendable es usar `qrhouse` y `qrsolve` que permiten resolver los pasos (*) y (**) simultáneamente.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000052
md"""
- Un problema más grave con este algoritmo es que los vectores ``v_j`` pueden no ser ortogonales por errores de cancelación. Si esto ocurre, el residuo y la solución aproximada pueden ser poco precisos. Un remedio parcial a este problema consiste en reemplazar la ortogonalización de Gram-Schmidt clásica de (*!) por la siguiente modificación

```octave
Para j = 1, ..., k
    w = w - (w' * vj) * vj;
Fin Para
v{k+1} = w;
```

Los dos algoritmos coinciden cuando la precisión de los cálculos es infinita, pero la forma modificada mantiene mejor la ortogonalidad en la práctica.

Para ver la influencia de esta mejora, basta tomar ``\delta=10^{-7}``, considerar la matriz ``A=\begin{bmatrix}1&1&1\\\delta&\delta&0\\\delta&0&\delta\end{bmatrix}`` y ortogonalizar usando MATLAB/OCTAVE utilizando Gram-Schmidt y su modificación. Hacer y ver cuánto se diferencia ``V^TV`` de la matriz identidad en cada caso.

Con la modificación propuesta a Gram-Schmidt, el algoritmo resulta:

**Método GMRes.** (Versión preliminar 2)

```octave
Dada A en R^{nxn}, x0, b en R^n, epsilon > 0 y kmax en N.
r0 = b - A*x0;
beta = rho = ||r0||_2;
v1 = r0 / rho;
k = 0;
Mientras rho > epsilon*||b||_2 y k < kmax, hacer
    k = k + 1;
    v{k+1} = A * vk;
    Para j = 1, ..., k, hacer
        h_{j,k} = v{k+1}' * vj;
        v{k+1} = v{k+1} - h_{j,k} * vj;
    Fin Para
    h_{k+1,k} = ||v{k+1}||_2;
    v{k+1} = v{k+1} / h_{k+1,k};
    yk = argmin_{y en R^k} ||beta*e1 - Hk*y||_2;      % (*)
    rho = ||beta*e1 - Hk*yk||_2;      % (**)
Fin Mientras
xk = x0 + Vk*yk;
```

Aún con la modificación a Gram-Schmidt, se puede perder la ortogonalidad en las columnas de ``V_k``. Se puede chequear la pérdida de ortogonalidad y reortogonalizar si es necesario, o usar métodos más estables para crear la matriz ``V``. Estas implementaciones más complejas son necesarias cuando ``A`` es mal condicionada o si ``k_{\max}`` es grande. No entraremos en más detalles sobre esto en este curso.

**Observación 5.45.** Una última observación es que la matriz ``H_{k+1}`` *contiene* a la ``H_k`` en el sentido de que se agrega una fila en cada iteración. Una manera de hacer más eficientes los pasos (*)-(**) se obtiene manteniendo la factorización QR de ``H_k`` y calculando solo lo necesario para obtener la factorización QR de ``H_{k+1}``. Con una implementación eficiente que tenga en cuenta esto, se logra que cada iteración de GMRes necesite ``O(nk)`` operaciones. Veremos esto en la siguiente sección.
"""

# ╔═╡ e5f6a7b8-0000-4a1b-8c2d-000000000003
md"""
**Nota** (de las diapositivas de la Clase 14). Vale la pena comparar explícitamente el costo por iteración de GC y de GMRes, para ver de dónde sale la diferencia de comportamiento asintótico:

| | GC (``A`` sdp) | GMRes (iteración ``k``) |
|---|---|---|
| Productos mat-vec | 1 | 1 |
| Productos escalares | 2 | ``k+1`` (ortogonalización más norma) |
| Operaciones axpy | 3 | ``k`` |
| Rotaciones de Givens | - | ``k``, de ``2\times2``: ``O(k)`` |
| Memoria | 4 vectores | ``k+1`` vectores, más ``O(k^2)`` para ``H`` |
| Costo total hasta la iteración ``k`` | ``O(nk)`` | ``O(nk^2)`` |

La diferencia la da la recurrencia de dos términos, que GC obtiene gratis por ser ``A`` sdp: sin ella, cada iteración de GMRes es más cara que la anterior (hay que ortogonalizar contra todos los vectores previos). Por eso, en GMRes es todavía más importante precondicionar bien: con pocas iteraciones es la única manera de mantener acotados la memoria y el costo.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000053
md"""
### 5.9.5 Implementación con rotaciones de Givens

En los pasos (*)-(**) de la versión preliminar 2 de GMRes tenemos que calcular
```math
y_k = \operatorname*{argmin}_{y\in\mathbb{R}^k}\|\beta e_1-H_ky\|_2, \qquad\text{y}\qquad \rho = \|\beta e_1-H_ky_k\|_2.
```
Por un lado, como ``H_k`` es upper Hessenberg de ``(k+1)\times k``, podemos hallar su factorización QR usando las rotaciones de Givens, como se explica en la Sección 3.6, ``Q`` será una composición o multiplicación de rotaciones de Givens de ``(k+1)\times(k+1)`` y ``R`` una matriz triangular superior de ``(k+1)\times k``.

Suponiendo que ``R_k=Q_k^TH_k``,
```math
\|\beta e_1-H_ky_k\|_2 = \|Q_k^T(\beta e_1-H_ky_k)\|_2 = \|\underbrace{Q_k^T\beta e_1}_{g}-\underbrace{Q_k^TH_k}_{R_k}y_k\|_2,
```
tenemos que:

- el mínimo se encuentra en ``y_k=(R_k)_{1:k,1:k}\backslash g_{1:k}``, resolviendo un sistema triangular;
- y el mínimo valor es ``\rho=|g_{k+1}|``.

Por último, es importante notar que no hace falta realizar completamente dichas factorizaciones en cada paso, por el siguiente motivo. Pensemos inductivamente:

``H_{k+1}`` se obtiene a partir de ``H_k`` agregándole una columna y completando con ceros:
```math
H_{k+1} = \left[\begin{array}{c|c}&h_{1,k+1}\\H_k&\vdots\\&h_{k+1,k+1}\\\hline0\ 0\ \cdots\ 0&h_{k+2,k+1}\end{array}\right] = \begin{bmatrix}H_k&h_{1:k+1,k+1}\\0^T&h_{k+2,k+1}\end{bmatrix}.
```
Luego, si ``H_k=Q_kR_k``,
```math
H_{k+1} = \begin{bmatrix}Q_kR_k&h_{1:k+1,k+1}\\0^T&h_{k+2,k+1}\end{bmatrix} = \begin{bmatrix}Q_k&0\\0^T&1\end{bmatrix}\begin{bmatrix}R_k&Q_k^Th_{1:k+1,k+1}\\0^T&h_{k+2,k+1}\end{bmatrix}.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000054
md"""
Esto es casi una factorización QR de ``H_{k+1}``, sólo nos falta hacer un cero donde está ``h_{k+2,k+1}``. Para ello, insertamos una rotación de Givens ``G_{k+1}:=G_{k+1}(c_{k+1},s_{k+1})`` y obtenemos
```math
\begin{align*}
H_{k+1} &= \begin{bmatrix}Q_k&0\\0^T&1\end{bmatrix}G_{k+1}^TG_{k+1}\begin{bmatrix}R_k&Q_k^Th_{1:k+1,k+1}\\0^T&h_{k+2,k+1}\end{bmatrix}\\
&= \begin{bmatrix}Q_k&0\\0^T&1\end{bmatrix}G_{k+1}^TR_{k+1}.
\end{align*}
```
Para ir terminando, observamos que ``G_{k+1}`` sólo afecta las filas ``k+1`` y ``k+2`` de la matriz a la que multiplica a su derecha, y además la última fila de ``R_k`` consta de todos ceros. Por lo tanto
```math
R_{k+1} = \left[\begin{array}{c|c}&\\R_k&G_{k+1}Q_k^Th_{1:k+2,k+1}\\&\\\hline0^T&\end{array}\right],
```
es decir, ``R_{k+1}`` se obtiene agregando una columna a ``R_k`` con ``G_{k+1}Q_k^Th_{k+1}``.

Por último, observamos que
```math
Q_{k+1} = G_1^TG_2^T\dots G_{k+1}^T = (G_{k+1}G_k\dots G_2G_1)^T. \tag{5.33}
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000055
md"""
Así llegamos a la siguiente versión eficiente del algoritmo de GMRes, donde en lugar de almacenar ``H_k`` sólo almacenamos ``R_k`` (almacenamos ``R_k`` en el lugar que ocuparía ``H_k``) y para acceder a ``Q_k`` vamos guardando los coeficientes ``(c_1,s_1),\dots,(c_k,s_k)`` de las matrices de Givens que dan lugar a ``Q_k`` según la fórmula (5.33)

**Método GMRes.** (Versión final)

```octave
Dada A en R^{nxn}, x0, b en R^n, epsilon > 0 y kmax en N.
r0 = b - A*x0;
beta = rho = ||r0||_2;
v1 = r0 / rho;
k = 0;
Mientras rho > epsilon*||b||_2 y k < kmax, hacer
    k = k + 1;
    v{k+1} = A * vk
    Para j = 1, ..., k, hacer
        h_{j,k} = v{k+1}' * vj;
        v{k+1} = v{k+1} - h_{j,k} * vj;
    Fin Para
    h_{k+1,k} = ||v{k+1}||_2;
    v{k+1} = v{k+1} / h_{k+1,k};
    Si k > 1, entonces
        h_{1:k+1,k} = Qk_1' * h_{1:k+1,k};      % reemplazar col_k(H) por Q_{k-1}' col_k(H)   (dagger)
    nu = sqrt(h_{k,k}^2 + h_{k+1,k}^2);
    c_k = h_{k,k} / nu;   s_k = -h_{k+1,k} / nu;      % k-esima rotacion de Givens Gk
    h_{k,k} = nu;   h_{k+1,k} = 0;      % aplicar Gk a col_k(H)
    g = Gk(c_k,s_k) * g;
    yk = trisup(H_{1:k,1:k}, g_{1:k});
    rho = |g_{k+1}|;      % (**)
Fin Mientras
xk = x0 + Vk*yk.
```

Por último observamos que en el paso (dagger) tenemos que aplicar ``Q_{k-1}^T`` al vector ``h_k``. Esto se implementa a través del sucesivo cálculo de rotaciones de Givens de ``2\times2`` a subvectores de ``h_k``, de la siguiente manera:
```math
h\leftarrow Q_{k-1}^Th = G_{k-1}G_{k-2}\dots G_2G_1h \iff \begin{cases}\text{Paso (1)} & h_{1:2}=G_1h_{1:2}\\\text{Paso (2)} & h_{2:3}=G_2h_{2:3}\\&\vdots\\\text{Paso }(k-1) & h_{k-1:k}=G_{k-1}h_{k-1:k}\end{cases}
```
Esto se puede programar de la siguiente forma
```octave
for j=1:k-1
    G = [c(j) -s(j); s(j) c(j)];
    H(j:j+1,k) = G*H(j:j+1,k);
end
```

Queda como ejercicio la implementación de un algoritmo que incorpore el *restart*.

Existen otros métodos iterativos basados en espacios de Krylov, pero escapan a este curso, y los vistos alcanzan para tener una idea sobre las bases de estos métodos.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000056
md"""
## 5.10 Problemas
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000057
md"""
### 5.1)

Sea ``A\in\mathbb{R}^{n\times n}`` sdp. Demostrar que existe una única matriz ``B\in\mathbb{R}^{n\times n}`` sdp tal que ``B^2=A``. En este caso, la matriz ``B`` se denota por ``A^{1/2}`` y se llama la *raíz cuadrada* de ``A``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000058
md"""
### 5.2)

Sea ``A`` una matriz sdp y ``b\in\mathbb{R}^n``, y sea ``\phi:\mathbb{R}^n\to\mathbb{R}`` definida por
```math
\phi(x) = \frac12x^TAx-b^Tx.
```
Demostrar que si ``x^*`` es la solución exacta de ``Ax=b``, entonces
```math
\phi(y) = \phi(x^*)+\frac12\|y-x^*\|_A^2.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000059
md"""
### 5.3)

Para resolver el sistema lineal ``Ax=b`` con el método del descenso más pronunciado y con el método del gradiente conjugado en MATLAB/OCTAVE, programar las funciones `[x, iter, X] = descmaspron(A,b,x0,tol,max_iter)` y `[x, iter, X] = gradconj(A,b,x0,tol,max_iter)`, respectivamente. Aquí, `x0` es el iterando inicial, `tol` es una tolerancia para la norma 2 del residuo, `max_iter` es la cantidad máxima de iteraciones permitidas, `x` es la aproximación obtenida con las `iter` iteraciones realizadas, y las columnas de la matriz `X` son las sucesivas aproximaciones generadas por el método. Realizar la programación de manera que el primer argumento `A` pueda ser tanto la matriz del sistema en forma explícita como una rutina que devuelve el producto matriz-vector para cada vector dado.

*Importante:* Puede calcular solo un producto matriz-vector por iteración.
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000001
md"""
**Resolución.**

Usando el algoritmo final del descenso más pronunciado (Observación 5.5) y el algoritmo final del gradiente conjugado (Sección 5.1.5), con `A` pudiendo ser tanto una matriz como una función que calcula `A*v`, y guardando las sucesivas aproximaciones en `X`:
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000002
function descmaspron(A, b, x0; tol=1e-8, max_iter=1000)
	Av(v) = A isa Function ? A(v) : A * v
	x = copy(x0)
	r = b - Av(x)
	X = reshape(copy(x), :, 1)
	iter = 0
	while norm(r) > tol && iter < max_iter
		Ar = Av(r)              # unico producto matriz-vector de la iteracion
		α = (r' * r) / (r' * Ar)
		x = x + α * r
		r = r - α * Ar          # = b - A*x, sin volver a multiplicar por A
		X = hcat(X, x)
		iter += 1
	end
	return x, iter, X
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000003
function gradconj(A, b, x0; tol=1e-8, max_iter=1000)
	Av(v) = A isa Function ? A(v) : A * v
	x = copy(x0)
	r = b - Av(x)
	ρ = r' * r
	X = reshape(copy(x), :, 1)
	d = copy(r)
	ρold = 0.0
	iter = 0
	while sqrt(ρ) > tol && iter < max_iter
		if iter > 0
			β = ρ / ρold
			d = r + β * d
		end
		w = Av(d)                # unico producto matriz-vector de la iteracion
		α = ρ / (d' * w)
		x = x + α * d
		r = r - α * w
		ρold = ρ
		ρ = r' * r
		X = hcat(X, x)
		iter += 1
	end
	return x, iter, X
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000004
md"""
En ambos casos, `x` es la última aproximación, `iter` la cantidad de iteraciones realizadas, y las columnas de `X` (con `x0` como primera columna) son las sucesivas aproximaciones. En cada iteración se hace un único producto matriz-vector: en `descmaspron` se calcula `Ar` una vez y se reutiliza tanto para `α` como para actualizar `r`; en `gradconj`, análogamente con `Ad` (llamado `w`).
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000005a
md"""
### 5.4)

Dada la matriz ``A=\begin{pmatrix}2&1\\1&3\end{pmatrix}`` y el vector ``b=\begin{pmatrix}1\\-2\end{pmatrix}``, dibujar las curvas de nivel de la energía del sistema ``Ax=b`` y sus gradientes. Observar las direcciones de los mismos. Intentar determinar visualmente la posición del mínimo de dicha función.

#### 5.4.a)
Determinar la solución del sistema mediante el método del descenso más pronunciado, tratando de visualizar las aproximaciones y las direcciones que se toman en cada paso. Considerar una tolerancia de ``10^{-3}``.

#### 5.4.b)
Determinar la solución del sistema mediante el método del gradiente conjugado, visualizando las aproximaciones y las direcciones que se toman en cada paso hasta llegar a la solución.
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000005
md"""
**Resolución.**

La solución exacta de ``Ax=b`` es ``x_*=(1,-1)`` (de ``2x+y=1``, ``x+3y=-2``). Los autovalores de ``A`` son
```math
\lambda_{1,2} = \frac{5\mp\sqrt5}2, \qquad \lambda_1\approx1.382, \quad \lambda_2\approx3.618, \qquad \kappa_2(A)=\frac{\lambda_2}{\lambda_1}=\frac{3+\sqrt5}2\approx2.618.
```

Como ``A`` no es diagonal, las curvas de nivel de ``\phi`` son elipses centradas en ``x_*``, pero *rotadas* respecto de los ejes coordenados: el semieje mayor apunta en la dirección del autovector asociado a ``\lambda_1`` (el autovalor más chico), aproximadamente ``(1,-0.618)``, y el semieje menor en la dirección del autovector asociado a ``\lambda_2``, aproximadamente ``(1,1.618)`` (perpendicular a la anterior). El gradiente ``\nabla\phi(x)=Ax-b=-r`` en cada punto es perpendicular a la curva de nivel que pasa por ese punto; el descenso más pronunciado avanza en la dirección ``-\nabla\phi=r``.
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000006
begin
	A_54 = [2.0 1.0; 1.0 3.0]
	b_54 = [1.0, -2.0]
	xstar_54 = A_54 \ b_54
	ev_54 = eigen(A_54)
	φ_54(x, y) = 0.5 * (2x^2 + 2x * y + 3y^2) - (x - 2y)
	xs_54 = range(-1, 3, length=150)
	ys_54 = range(-3, 1, length=150)
	plot_54_base = Plots.contour(xs_54, ys_54, φ_54, levels=25, aspect_ratio=:equal,
		xlabel="x", ylabel="y", legend=:topleft, title="Curvas de nivel de φ")
	Plots.scatter!(plot_54_base, [xstar_54[1]], [xstar_54[2]], color=:red, markersize=6, label="x*")
	plot_54_base
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000007
md"""
**5.4.a)** Descenso más pronunciado, con ``x_0=(0,0)`` y tolerancia ``10^{-3}`` sobre ``\|r_k\|_2``:
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000008
begin
	xsd_54, itsd_54, Xsd_54 = descmaspron(A_54, b_54, [0.0, 0.0]; tol=1e-3, max_iter=1000)
	plot_54_sd = deepcopy(plot_54_base)
	Plots.plot!(plot_54_sd, Xsd_54[1, :], Xsd_54[2, :], marker=:circle, color=:blue, label="SD")
	plot_54_sd
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000009
md"""
El descenso más pronunciado necesita **``$(itsd_54)`` iteraciones**. Se ve el zigzagueo característico del método: cada residuo es ortogonal al anterior (``r_k\perp r_{k-1}``, Observación 5.6), y por la simetría particular de este ejemplo, las iteraciones impares caen exactamente sobre la recta ``y=-1`` mientras que las pares se acercan a ``x_*`` desde la diagonal ``y=-x`` (columnas de `Xsd_54` arriba).

**5.4.b)** Gradiente conjugado, con el mismo ``x_0`` y tolerancia:
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-00000000000a
begin
	xcg_54, itcg_54, Xcg_54 = gradconj(A_54, b_54, [0.0, 0.0]; tol=1e-3, max_iter=1000)
	plot_54_cg = deepcopy(plot_54_base)
	Plots.plot!(plot_54_cg, Xcg_54[1, :], Xcg_54[2, :], marker=:square, color=:green, label="GC")
	plot_54_cg
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-00000000000b
md"""
Como ``A`` tiene ``k=2`` autovalores distintos, el Teorema 5.21 garantiza que el gradiente conjugado converge en a lo sumo 2 pasos, cualquiera sea ``x_0``; acá converge en **``$(itcg_54)`` iteraciones** (columnas de `Xcg_54`). Notar que el primer paso coincide con el de descenso más pronunciado (siempre ocurre, porque ``d_1=r_0`` en ambos métodos), pero el gradiente conjugado llega a la solución exacta en el segundo paso, en vez de necesitar las ``$(itsd_54)`` iteraciones de (a).
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000005b
md"""
### 5.5)

Consideramos el problema de hallar la función ``u`` que satisface
```math
\begin{cases}-u''(x)=f(x), & x\in(0,1)\\u(0)=u(1)=0,\end{cases}
```
donde ``f(x)=\pi^2\,\text{sen}(\pi x)``. Se puede verificar fácilmente que la solución de este problema es ``u(x)=\text{sen}(\pi x)``. Dada una partición del intervalo ``[0,1]`` en ``N+1`` intervalos de igual longitud ``[x_0,x_1],[x_1,x_2],[x_2,x_3],\dots,[x_N,x_{N+1}]`` donde ``x_0=0`` y ``x_{N+1}=1``, la discretización de este problema por diferencias finitas se reduce a resolver el sistema ``Ax=b`` donde ``A`` es la matriz tridiagonal cuadrada de tamaño ``N`` que tiene 2 sobre la diagonal principal y ``-1`` en las otras posiciones no-nulas, y el vector ``b=h^2(f(x_1),f(x_2),\dots,f(x_N))^T``, donde ``h=\frac{1}{N+1}`` es la longitud de cada subintervalo.

El vector solución ``\mathbf{x}`` aproxima los valores de la solución ``u`` en los puntos ``x_1,x_2,\dots,x_N``, es decir, ``\mathbf{x}\approx(u(x_1),u(x_2),\dots,u(x_N))^T``.

Considerando distintos valores de ``N`` (``N=10,50,100,\dots``):

#### 5.5.a)
Cargar la matriz ``A`` en formato sparse, utilizando la función `spdiags` de MATLAB/OCTAVE. Cargar también el vector de términos independientes ``b``.

#### 5.5.b)
Utilizar los métodos de Jacobi, Gauss-Seidel y SOR para resolver el sistema lineal. Realizar un gráfico del error (en escala logarítmica) de cada método en función de la cantidad de iteraciones. Obtener conclusiones.

*Sugerencia:* Para calcular el error, utilizar como solución exacta la solución del sistema dada por el comando ``\backslash`` de MATLAB/OCTAVE.

#### 5.5.c)
Analizar para diversos valores de ``\alpha`` el comportamiento del método de Richardson. Comparar los resultados obtenidos en el item anterior.

#### 5.5.d)
Investigar el comportamiento de los métodos del descenso más pronunciado y del gradiente conjugado para resolver el sistema lineal.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000005c
md"""
### 5.6)

Crear una función `matvec` que para cada vector ``x`` devuelva el resultado de ``Ax``, sin realizar el producto de manera explícita, donde
```math
A = \begin{bmatrix}2&-1&&&\\-1&2&-1&&0\\&-1&2&-1&\\&&\ddots&\ddots&\ddots\\&0&&-1&2&-1\\&&&&-1&2\end{bmatrix}.
```

*Ayuda:* Tratar de interpretar la acción de ``A`` sobre un vector dado ``x``, y "transformarlo" según esta acción sin utilizar la matriz ``A``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000005d
md"""
### 5.7)

Consideramos el problema de hallar la función ``u`` que satisface el siguiente problema de valores de borde
```math
\begin{cases}-\Delta u(x,y)=f(x,y), & (x,y)\in(0,1)\times(0,1)\\u(x,0)=u(x,1)=0, & x\in[0,1]\\u(0,y)=u(1,y)=0, & y\in[0,1]\end{cases}
```
donde ``f(x,y)=-2[x(x-1)+y(y-1)]``. Si consideramos una malla uniforme de tamaño ``h=\frac{1}{N+1}``, sobre la que los nodos interiores se enumeran de izquierda a derecha y de abajo hacia arriba, la discretización de este problema por diferencias finitas se reduce a resolver el sistema ``Ax=b`` donde ``A`` es la matriz cuadrada de tamaño ``N^2`` expresada en bloques como
```math
A = \begin{pmatrix}B&-I_N&&&\\-I_N&B&-I_N&&\\&-I_N&B&-I_N&\\&&\ddots&\ddots&\ddots\\&&&-I_N&B&-I_N\\&&&&-I_N&B\end{pmatrix},
```
donde la matriz ``B`` de tamaño ``N`` está dada por ``B=\begin{pmatrix}4&-1&&&\\-1&4&-1&&\\&-1&4&\ddots&\\&&\ddots&\ddots&-1\\&&&-1&4\end{pmatrix}``, e ``I_N`` es la matriz identidad de orden ``N``. Las componentes del vector de términos independientes ``b`` son ``h^2`` por el valor de ``f`` en los correspondientes puntos de la grilla (enumerados de izquierda a derecha y de abajo hacia arriba).

El vector solución ``\mathbf{x}`` aproxima los valores de la solución ``u`` en los respectivos puntos de la grilla.

Considerando distintos valores de ``N`` (``N=10,50,100,\dots``), repetir los items (a), (b), (c) y (d) del Ejercicio 5.5.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000005e
md"""
### 5.8)

Consideremos la siguiente cota (5.6) del Teorema 5.19 para las iteraciones del método de gradiente conjugado:
```math
\|x_k-x_*\|_A \leq 2\|x_0-x_*\|_A\left(\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\right)^k.
```
Este ejercicio consiste en ver que esta cota puede ser demasiado pesimista en algunos casos. Veamos el siguiente ejemplo: Supongamos que los autovalores de una matriz simétrica ``A`` se encuentran contenidos en ``[1,1.5]\cup[399,400]=:I``.

#### 5.8.a)
Dar una estimación de ``\kappa_2(A)`` y calcular ``\dfrac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}``.

#### 5.8.b)
Tomar el polinomio ``p_{3k}(z)=\dfrac{(1.25-z)^k(400-z)^{2k}}{(1.25)^k(400)^{2k}}`` y considerar la estimación
```math
\|x_{3k}-x\| \leq \|x_0-x\|_A\max_{z\in I}|p_{3k}(z)|.
```
¿Cuál de las dos cotas es mejor? Graficar el comportamiento de las cotas como función de ``k`` y sacar conclusiones.
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-00000000000c
md"""
**Resolución.**

**5.8.a)** Como ``\sigma(A)\subset[1,1.5]\cup[399,400]``, la única información disponible da ``\lambda_{\min}\geq1`` y ``\lambda_{\max}\leq400``, así que la estimación (en el peor caso posible) es
```math
\kappa_2(A) \leq \frac{400}{1} = 400, \qquad \frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1} \leq \frac{\sqrt{400}-1}{\sqrt{400}+1} = \frac{19}{21} \approx 0.9048.
```

**5.8.b)** La cota general del Teorema 5.19, evaluada en ``3k`` iteraciones, da ``\text{cota A}(k) := 2\left(\frac{19}{21}\right)^{3k}``, mientras que la cota del enunciado (que es exactamente el Corolario 5.14 aplicado a este ``p_{3k}`` particular) da ``\text{cota B}(k) := \max_{z\in[1,1.5]\cup[399,400]}|p_{3k}(z)|``:
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-00000000000d
begin
	κ_58 = 400.0
	r_58 = (sqrt(κ_58) - 1) / (sqrt(κ_58) + 1)
	cotaA_58(k) = 2 * r_58^(3k)
	function cotaB_58(k)
		p(z) = ((1.25 - z)^k * (400.0 - z)^(2k)) / (1.25^k * 400.0^(2k))
		z1 = range(1, 1.5, length=2000)
		z2 = range(399, 400, length=2000)
		max(maximum(abs.(p.(z1))), maximum(abs.(p.(z2))))
	end
	ks_58 = 1:30
	plot_58 = Plots.plot(ks_58, cotaA_58.(ks_58), yscale=:log10, marker=:circle,
		label="cota general (Teo. 5.19)", xlabel="k", ylabel="cota para ||x_{3k}-x*||_A",
		legend=:bottomleft)
	Plots.plot!(plot_58, ks_58, cotaB_58.(ks_58), marker=:square, label="cota especifica (p_{3k})")
	plot_58
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-00000000000e
md"""
La cota B (específica) es muchísimo mejor, y la diferencia crece sin límite con ``k``: por ejemplo, `cotaA_58(10)` da ≈ `round(cotaA_58(10),sigdigits=4)` mientras que `cotaB_58(10)` da ``\approx$(round(cotaB_58(10),sigdigits=4))``, varios órdenes de magnitud de diferencia. Esto es exactamente el fenómeno de la Observación 5.18 y el Ejemplo 5.17 (autovalor aislado): el polinomio ``p_{3k}`` se construye a medida, con un factor que se anula en el centro del bloque aislado ``[399,400]`` combinado con un factor tipo Tchebychev sobre el intervalo agrupado ``[1,1.5]``, y así explota la estructura del espectro. La cota general del Teorema 5.19, en cambio, solo usa ``\kappa_2(A)`` y por lo tanto es ciega a que el espectro está agrupado salvo por un único autovalor aislado. La conclusión práctica, como ya señalaba la Observación 5.18, es que si se conoce algo más sobre ``\sigma(A)`` que su cociente ``\kappa_2(A)``, conviene explotarlo: la cota basada solo en ``\kappa_2(A)`` puede ser enormemente pesimista.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000005f
md"""
### 5.9)

Sea ``A\in\mathbb{R}^{n\times n}``. Demostrar que si ``S\in\mathbb{R}^{n\times n}`` es invertible, entonces ``\sigma(S^2A)=\sigma(SAS)``.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000060
md"""
### 5.10)

Sean ``A`` y ``S`` sdp, y ``M=S^2``. Demostrar que si
```math
c_1 \leq \frac{x^TMx}{x^TA^{-1}x} \leq c_2, \qquad \forall x\in\mathbb{R}^n\setminus\{0\},
```
entonces
```math
\kappa_2(SAS) \leq \frac{c_2}{c_1}.
```

*Observación:* Cuando se buscan precondicionadores ``M`` que sean "cercanos" a ``A^{-1}``, se entiende que sean cercanos en este sentido:
```math
c_1 \leq \frac{x^TMx}{x^TA^{-1}x} \leq c_2, \qquad \forall x\in\mathbb{R}^n\setminus\{0\}, \qquad\text{con } \frac{c_2}{c_1} \text{ pequeño}.
```
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000061
md"""
### 5.11)

#### 5.11.a)
Modificar la función `gradconj` del Ejercicio 5.3 para implementar el método del gradiente conjugado precondicionado a través de la función `[x, iter, X] = gradconjp(A,b,precon,x0,tol,max_iter)`. Aquí, `precon` es el precondicionador a usar, dado como matriz en forma explícita o como una rutina que devuelve el producto matriz-vector para cada vector dado.

#### 5.11.b)
Comparar el comportamiento del método del gradiente conjugado usual y del método precondicionado con el precondicionador de Cholesky incompleto para la resolución del sistema del Ejercicio 5.5 y/o 5.7.
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-00000000000f
md"""
**Resolución.**

**5.11.a)** Basta agregar el precondicionador a `gradconj` del Ejercicio 5.3, siguiendo el algoritmo del gradiente conjugado precondicionado (Sección 5.1.6):
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000010
function gradconjp(A, b, precon, x0; tol=1e-8, max_iter=1000)
	Av(v) = A isa Function ? A(v) : A * v # por si es rutina
	Mv(v) = precon isa Function ? precon(v) : precon * v #idem
	
	x = copy(x0)
	r = b - Av(x)
	z = Mv(r)
	ρ = r' * z
	X = reshape(copy(x), :, 1)
	d = copy(z)
	ρold = 0.0
	R = [norm(r)]
	
	iter = 0
	while norm(r) > tol && iter < max_iter
		if iter > 0
			β = ρ / ρold
			d = z + β * d
		end
		w = Av(d)
		α = ρ / (d' * w)
		x = x + α * d
		r = r - α * w
		push!(R, norm(r))
		z = Mv(r)
		ρold = ρ
		ρ = r' * z
		X = hcat(X, x)
		iter += 1
	end
	return x, iter, X, R
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000011
md"""
Con `precon = v -> v` (es decir, ``M=I``) se recupera exactamente `gradconj` del Ejercicio 5.3.
"""

# ╔═╡ c36e1dac-583a-4ace-b6d6-854d7e60f034
md"""
**5.11.b)** Para la matriz tridiagonal del Ejercicio 5.5 (el laplaciano discreto 1D), la factorización de Cholesky *no tiene fill-in*: al ser tridiagonal, el factor exacto ``R`` ya es bidiagonal, así que la factorización incompleta coincide con la completa, y el precondicionador resultante es esencialmente ``M=A^{-1}``: la comparación queda degenerada (GCP converge casi de inmediato). El caso interesante es el del Ejercicio 5.7 (laplaciano 2D de 5 puntos), donde la factorización de Cholesky sí tiene *fill-in* dentro de la banda, y la versión incompleta (que descarta ese fill-in conservando el patrón de dispersión de ``A``) produce un precondicionador genuinamente no trivial:
"""

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000012
"""
	icholesky(A)
"""
function icholesky(A)
	n = size(A, 1)
	L = zeros(n, n)
	for k in 1:n
		L[k, k] = sqrt(A[k, k] - L[k, 1:k-1]' * L[k, 1:k-1])
		for i in k+1:n
			if A[i, k] != 0
				L[i, k] = (A[i, k] - L[i, 1:k-1]' * L[k, 1:k-1]) / L[k, k]
			end
		end
	end
	return LowerTriangular(L)
end

# ╔═╡ 3ab565b6-8975-49e3-b6fe-de57ebc518d0
function laplaciano2d_511(n)
	dn(x; v=-1.0) = x => fill(v, n^2-abs(x))
	A = diagm(dn(0, v=4.0), dn.([1, -1, n, -n])...)

	for i in 1:(n-1)
		A[n*i+1, n] = 0
		A[n, n*i+1] = 0
	end
	
	return A
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000014
begin
	Ns_511 = [8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52]
	filas_511 = map(Ns_511) do N
		A_511 = laplaciano2d_511(N)
		n_511 = N * N
		b_511 = ones(n_511)
		L_511 = icholesky(A_511)
		precon_511 = v -> L_511' \ (L_511 \ v)
		precon_jacobi = diagm(0 => diag(A_511))
		
		_, it_gc_511, _ = gradconj(A_511, b_511, zeros(n_511); tol=1e-8, max_iter=5000)
		
		_, it_gcp_511, _, res_ichol = gradconjp(A_511, b_511, precon_511, zeros(n_511); tol=1e-8, max_iter=5000)

		_, it_jacobi_511, _, res_jaco = gradconjp(A_511, b_511, precon_jacobi, zeros(n_511); tol=1e-8, max_iter=5000)
		
		
		(N = N, 
		 	n = n_511, 
		 	GC = it_gc_511, 
		 	GCP_icholesky = it_gcp_511,
		 	res_ichol = res_ichol,
		 	res_jaco = res_jaco,
		 	GCP_jacobi = it_jacobi_511)
	end
	filas_511
end

# ╔═╡ d0243450-c789-4523-b239-eab2dff062ad
let
	Ns_511 = [8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52]
	filas_511 = map(Ns_511) do N
		A_511 = diagm(
			-1 => fill(-1, N-1), 
			0 => fill(2, N), 
			1 => fill(-1, N-1))
		
		n_511 = N
		b_511 = ones(n_511)
		L_511 = icholesky(A_511)
		precon_511 = v -> L_511' \ (L_511 \ v)
		precon_jacobi = diagm(0 => diag(A_511))
		
		_, it_gc_511, _ = gradconj(A_511, b_511, zeros(n_511); tol=1e-8, max_iter=5000)
		
		_, it_gcp_511, _ = gradconjp(A_511, b_511, precon_511, zeros(n_511); tol=1e-8, max_iter=5000)

		_, it_jacobi_511, _ = gradconjp(A_511, b_511, precon_jacobi, zeros(n_511); tol=1e-8, max_iter=5000)
		
		
		(N = N, 
		 	n = n_511, 
		 	GC = it_gc_511, 
		 	GCP_icholesky = it_gcp_511,
		 	GCP_jacobi = it_jacobi_511)
	end
	filas_511
end

# ╔═╡ c67999bc-c433-4430-87a5-06bd59202136
let
	fig = Figure(size=(700,600))
	ax = Axis(fig[1,1])

	inc = [itm.n for itm in filas_511]
	GC = [itm.GC for itm in filas_511]
	GCP_chol = [itm.GCP_icholesky for itm in filas_511]
	GCP_jaco = [itm.GCP_jacobi for itm in filas_511]

	lines!(ax, inc, GC, label = "Grad. conjugado")
	scatter!(ax, inc, GCP_jaco, label = "Grad. conj. `jaco`", color=:red)
	lines!(ax, inc, GCP_chol, label = "Grad. conj. `ichol`")
	
	axislegend()
	fig
end

# ╔═╡ 99db9aa9-58c2-478e-b55c-b2c2c48e4176
let
	fig = Figure(size=(700,500))
	ax = Axis(fig[1,1])

	inc = [itm.n for itm in filas_511]
	res_ichol = [itm.res_ichol for itm in filas_511]
	res_jaco = [itm.res_jaco for itm in filas_511]

	for res in res_ichol
		lines!(ax, eachindex(res), res)
	end


	for res in res_jaco
		lines!(ax, eachindex(res), -res)
	end
	
	
	
	fig
end

# ╔═╡ e5f6a7b8-ffff-4a1b-8c2d-000000000015
md"""
Comparando la cantidad de iteraciones de `gradconj` contra `gradconjp` con este precondicionador (calculado una sola vez, antes de iterar), para alcanzar ``\|r_k\|_2<10^{-8}`` (ver la tabla `filas_511` arriba, con ``N`` puntos por lado y ``n=N^2`` incógnitas): para ``N`` chico, el precondicionador no ayuda demasiado (incluso necesita alguna iteración más, porque el sistema ya es chico y bien condicionado), pero a medida que ``N`` crece, las iteraciones de GC crecen aproximadamente en proporción a ``N`` (consistente con ``\kappa_2(A)=O(N^2)`` y la Observación 5.20), mientras que las de GCP con Cholesky incompleta crecen mucho más despacio. Este es el comportamiento típico de un precondicionador genérico (Sección 5.1.6): mejora la constante y aplana la curva, pero no logra el ``O(1)`` independiente de ``N`` que sí logra multigrid (Sección 5.1.7 y Problema 5.14), precisamente porque la factorización incompleta no usa *de dónde viene* la matriz, solo su patrón de ceros.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000062
md"""
### 5.12) (Multigrid como método iterativo.)

Implementar en MATLAB/OCTAVE las funciones `jacobiamortiguado`, `restringir`, `prolongar` y `vciclo` de la Sección 5.1.7.

Para los experimentos, tomar ``n=2^L-1`` y ``h=\frac{1}{n+1}``, elegir un vector ``u_*\in\mathbb{R}^n`` cualquiera (por ejemplo, ``(u_*)_i=x_i(1-x_i)\,\text{sen}(\pi x_i)``, o un vector aleatorio) y definir ``f=A_hu_*``. De este modo se conoce *exactamente* la solución del sistema discreto ``A_hu=f``, y se puede medir el error verdadero y no solo el residuo.
"""

# ╔═╡ ca9731d4-fb74-4e3c-a9a8-e5cab1bdf254
"""
	function u = jacobiamortiguado(u, f, h, ν)

Aplica ν barridos de Jacobi amortiguado (ω = 2/3) al
sistema de Poisson 1D `A_h u = f` sobre la malla interior de paso `h`.

 - `u`, `f` son vectores columna de largo `n` (nodos interiores).
"""
function jacobiamortiguado(u, f, h, nu)
    n = length(u); ω = 2/3
    for _ in 1:nu
        Au = (2u .- [u[2:n];0] .- [0;u[1:n-1]]) / h^2
        u .+= ω * (h^2/2) .* (f .- Au)
    end
    return u
end

# ╔═╡ 20786201-ff40-462a-87d5-421b75787a0e
function restringir(r)
    n = length(r); nc = (n-1)÷2
    return [(r[2j-1] + 2r[2j] + r[2j+1])/4 for j in 1:nc]
end

# ╔═╡ f606ece1-f7a4-4fd7-a4e9-be28cf2c3be3
function prolongar(vc)
    nc = length(vc); n = 2nc+1
    vext = [0; vc; 0]
    v = zeros(n)
    v[2:2:n-1] = vc
    v[1:2:n] = [(vext[j+1]+vext[j+2])/2 for j in 0:nc]
    return v
end

# ╔═╡ 88d558cf-2e3e-444e-8bf2-0901bf0a77f0
function vciclo(u, f, h, nu1, nu2)
    n = length(u)
    n == 1 && return f * h^2 / 2        # malla más gruesa: A_h es escalar 2/h²
    
    u = jacobiamortiguado(u, f, h, nu1)         # (1) pre-suavizado
    Au = (2u .- [u[2:n];0] .- [0;u[1:n-1]]) / h^2
    r = f - Au                                   # (2) residuo
    rc = restringir(r)                           # (3) restricción
    ec = vciclo(zeros(length(rc)), rc, 2h, nu1, nu2)   # (4) corrección (recursiva)
    u += prolongar(ec)                           # (5) prolongación y corrección
    u = jacobiamortiguado(u, f, h, nu2)          # (6) post-suavizado
    return u
end

# ╔═╡ fc934079-d098-4bab-b648-e3312c0902f8
begin 
	function Ah(N)
		    h = 1/(N+1)
		    SymTridiagonal(fill(2/h^2, N), fill(-1/h^2, N-1))
		end
	
	function uexacta(N)
		h = 1/(N+1)
		[(i*h)*(1-i*h)*sin(π*i*h) for i in 1:N]
	end
end

# ╔═╡ 19adc32b-6e30-4349-a635-77757a8efd70
md"""
#### 5.12.a) *Test de consistencia.*
Verificar que si se parte de ``u=u_*``, un ciclo en V devuelve ``u_*`` (a menos de error de redondeo). Explicar por qué esto tiene que ser así, y por qué es una buena manera de detectar errores de programación.
"""

# ╔═╡ cbb452d2-fe7a-40b6-bc31-eb36ab32b121
let 	
	L = 10
	N = 2^L-1; h = 1/(N+1)
	A = Ah(N)
	u_exacto = uexacta(N)
	f = A * u_exacto
	u = vciclo(copy(u_exacto), f, h, 1, 1)
	@info norm(u - u_exacto)
end

# ╔═╡ e1de6880-500f-4386-a599-5731b3cfadd1
md"""
#### 5.12.b) *Test de la relación variacional.*
Construir las matrices de ``I_{2h}^h`` e ``I_h^{2h}`` aplicando las funciones a los vectores de la base canónica, y verificar numéricamente la identidad (5.23).
"""

# ╔═╡ 85a8dbf3-cdc7-43e8-8b5b-aae5b2d0d080
let
	N = 15
	nc = (N-1) ÷ 2
	canonico(n, j) = (e = zeros(n); e[j] = 1.0; e)
	
	I2h_h = hcat([prolongar(canonico(nc, j)) for j in 1:nc]...)
	Ih_2h = hcat([restringir(canonico(N, j)) for j in 1:N]...)
	@info norm(Ih_2h - 0.5*I2h_h')
end

# ╔═╡ 6dbfe708-ef86-4295-b6f7-1e99c8edbca9
md"""
#### 5.12.c)
Partiendo de ``u=0``, iterar ``u\leftarrow\text{vciclo}(u,f,h,\nu,\nu)`` hasta que ``\frac{\|f-A_hu\|_2}{\|f\|_2}<10^{-10}``, con ``\nu=1``. Registrar la cantidad de ciclos necesarios para ``L=5,6,\dots,11``. ¿Crece con ``n``?
"""

# ╔═╡ d2859751-1875-4e25-b383-1386195bd61c
for L in 5:11
    N = 2^L - 1; h = 1/(N+1)
    A = Ah(N)
    f = A * uexacta(N)
    u = zeros(N)
    k = 0
    while norm(f - A*u)/norm(f) > 1e-10
        u = vciclo(u, f, h, 1, 1)
        k += 1
    end
    println("L=$(rpad(L, 2))  n=$(rpad(N, 4))  ciclos=$k")
end

# ╔═╡ d38355e2-5150-40a3-99ed-0988ed3f1013
md"""
#### 5.12.d)
Calcular, en cada ciclo, el cociente ``\dfrac{\|r_{k+1}\|_2}{\|r_k\|_2}``. Verificar que se estabiliza en un valor ``\rho`` que casi no depende de ``h``, y comparar los valores de ``\rho`` obtenidos con ``\nu=1,2,3``. Contrastar con el Teorema 5.34.
"""

# ╔═╡ 38ac5ce1-82f5-4e7d-903b-20e9e01950bc
let
    function rho_medido(N, nu; maxciclos=50, tol=1e-8)
        h = 1/(N+1)
        A = Ah(N)
        f = A * uexacta(N)
        nf = norm(f)
        u = zeros(N)
        prev = norm(f - A*u)
        ratio = NaN
        for _ in 1:maxciclos
            u = vciclo(u, f, h, nu, nu)
            cur = norm(f - A*u)
            cur/nf < tol && break    # parar ANTES de que el residuo caiga al piso de
                                       # precisión de punto flotante (~1e-11): más allá de
                                       # eso, el cociente deja de medir ρ y empieza a medir
                                       # ruido de redondeo (probarlo: con maxciclos grande
                                       # y sin este corte, salen cocientes >1, inestables)
            ratio = cur/prev
            prev = cur
        end
        return ratio
    end

    nus = 1:1:50
    img = []
    
    for nu in nus, L in [8]
        N = 2^L - 1
        ρm = rho_medido(N, nu)
        push!(img, ρm)
        println("ν=$nu  L=$L  ρ≈", round(rho_medido(N, nu), digits=4))
    end

    fig = Figure(size=(700,400))
    ax = Axis(fig[1,1])

    lines!(ax, nus, img)

    fig
end

# ╔═╡ c5d20dfe-eb1f-41ce-8cc3-dd877c32aa32
md"""
#### 5.12.e) *Comparación con el suavizador solo.*
Repetir el item (c) usando únicamente `jacobiamortiguado` (es decir, sin corrección en malla gruesa), con una tolerancia más modesta como ``10^{-6}`` y valores chicos de ``L``. Verificar que la cantidad de barridos crece como ``O(n^2)``, y explicar la diferencia a la luz de la propiedad de suavizado.
"""

# ╔═╡ e18b7fd5-1937-4df9-bff7-672690541421
let
    function barridos_necesarios(N, nu; maxbarridos=200_000)
        h = 1/(N+1)
        A = Ah(N)
        f = A * uexacta(N)
        u = zeros(N)
        k = 0
        while norm(f - A*u)/norm(f) > 1e-10 && k < maxbarridos
            u = jacobiamortiguado(u, f, h, nu)
            k += 1
        end
        return k
    end
    
    for L in 4:7
        N = 2^L - 1
        println("L=$L  n=$N  barridos de Jacobi solo=", barridos_necesarios(N, 1))
    end
end

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000063
md"""
### 5.13) (Multigrid como precondicionador del gradiente conjugado.)

Usando las funciones del Problema 5.12 y la implementación del gradiente conjugado precondicionado de la Sección 5.1.6, tomar como precondicionador la acción
```math
\xi \longmapsto M\xi = \texttt{vciclo}(0,\xi,h,\nu,\nu),
```
es decir, un ciclo en V con dato inicial nulo.

#### 5.13.a)
Construir la matriz ``M`` aplicando esa acción a los vectores de la base canónica (para ``n`` chico, digamos ``n=31``). Verificar numéricamente que ``M`` es simétrica cuando ``\nu_1=\nu_2``, y que *deja de serlo* si se toma ``\nu_1\neq\nu_2``. ¿Por qué importa?

#### 5.13.b)
Calcular ``\sigma(MA_h)`` y verificar que está contenido en un intervalo alrededor de 1. Comparar ``\kappa_2(MA_h)=\frac{c_2}{c_1}`` con la cota ``\frac{1+\rho}{1-\rho}`` de (5.27), usando el ``\rho`` medido en el Problema 5.12(d).

#### 5.13.c)
Resolver el sistema con gradiente conjugado precondicionado, con la misma tolerancia que antes, para ``L=5,6,\dots,12``. Verificar que la cantidad de iteraciones *no crece* con ``n``. Comparar ``\nu=1`` con ``\nu=2``.

#### 5.13.d)
Verificar que la cota de error del gradiente conjugado precondicionado, con el ``\frac{c_2}{c_1}`` medido en (b), es consistente con la cantidad de iteraciones observada.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000064
md"""
### 5.14) (Comparación de los cuatro métodos.)

Para el mismo sistema de los dos problemas anteriores y para ``L=5,6,\dots,12``, armar una tabla con la cantidad de iteraciones que necesita cada uno de los siguientes métodos para alcanzar una tolerancia fija:

(i) gradiente conjugado sin precondicionar;

(ii) gradiente conjugado con precondicionamiento diagonal (Jacobi);

(iii) multigrid solo, como método iterativo (Problema 5.12);

(iv) gradiente conjugado precondicionado con multigrid (Problema 5.13).

#### 5.14.a)
Graficar, en escala logarítmica en ambos ejes, la cantidad de iteraciones en función de ``n``. Identificar cuáles curvas tienen pendiente ``\approx1`` y cuáles son horizontales, y relacionarlo con las cotas teóricas ``O(\sqrt{\kappa_2(A_h)})=O(n)`` y ``O(1)``.

#### 5.14.b)
Observar qué ocurre con el precondicionamiento diagonal. *Sugerencia:* mirar la diagonal de ``A_h``. ¿Qué dice esto sobre cuándo conviene usar un precondicionador genérico?

#### 5.14.c)
Comparar (iii) con (iv): ¿cuántos ciclos en V usa cada uno? Dado que el costo de una iteración de (iv) es esencialmente el de un ciclo en V más el de un producto mat-vec, ¿cuál de los dos resulta más barato? Interpretar el resultado: el gradiente conjugado actúa como un *acelerador* de multigrid.

#### 5.14.d)
Estimar el costo total de cada método en operaciones (recordando que un ciclo en V cuesta ``O(n)`` y que una iteración de gradiente conjugado también) y verificar que solo (iii) y (iv) alcanzan el costo óptimo ``O(n)``.

#### 5.14.e) (Opcional.)
Repetir las mediciones con la matriz del Problema 5.18, que incluye un término de convección y *no* es simétrica. ¿Qué falla? ¿Qué habría que cambiar?
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000065
md"""
### 5.15)

Programar el algoritmo de GMRes con las rotaciones de Givens para que sea lo más eficiente posible.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000066
md"""
### 5.16)

Programar el algoritmo de GMRes con *Restart*.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000067
md"""
### 5.17)

Comparar el comportamiento del método del gradiente conjugado y de GMRes para resolver el sistema lineal dado en el ejercicio 5.5.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000068
md"""
### 5.18)

Consideramos el problema de hallar la función ``u`` que satisface
```math
\begin{cases}-u''(x)+u'=f(x), & x\in(0,1)\\u(0)=u(1)=0,\end{cases}
```
donde ``f(x)=\pi^2\,\text{sen}(\pi x)+\pi\cos(\pi x)``. Se puede verificar fácilmente que la solución de este problema es ``u(x)=\text{sen}(\pi x)``. Dada una partición del intervalo ``[0,1]`` en ``N+1`` intervalos de igual longitud ``[x_0,x_1],[x_1,x_2],[x_2,x_3],\dots,[x_N,x_{N+1}]`` donde ``x_0=0`` y ``x_{N+1}=1``, la discretización de este problema por diferencias finitas se reduce a resolver el sistema ``Ax=b`` donde
```math
A = \frac{1}{h^2}\begin{bmatrix}2&-1&&&\\-1&2&-1&&0\\&-1&2&-1&\\&&\ddots&\ddots&\ddots\\&0&&-1&2&-1\\&&&&-1&2\end{bmatrix} + \frac{1}{2h}\begin{bmatrix}0&1&&&\\-1&0&1&&0\\&-1&0&1&\\&&\ddots&\ddots&\ddots\\&0&&-1&0&1\\&&&&-1&0\end{bmatrix},
```
y el vector ``b=(f(x_1),f(x_2),\dots,f(x_N))^T``, donde ``h=\frac{1}{N+1}`` es la longitud de cada subintervalo.

El vector solución ``\mathbf{x}`` aproxima los valores de la solución ``u`` en los puntos ``x_1,x_2,\dots,x_N``, es decir, ``\mathbf{x}\approx(u(x_1),u(x_2),\dots,u(x_N))^T``.

Considerando distintos valores de ``N`` (``N=10,50,100,\dots``), comparar el comportamiento de GMRes, CGNE, CGNR para resolver el sistema lineal.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-000000000069
md"""
### 5.19)

Demostrar el Corolario 5.40.
"""

# ╔═╡ e5f6a7b8-cccc-4a1b-8c2d-00000000006a
md"""
### 5.20)

Demostrar el Corolario 5.41.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CairoMakie = "~0.15.15"
Plots = "~1.41.7"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.0"
manifest_format = "2.1"
project_hash = "e7a597ae0b2549f26610f06c0e10a3232ae9bd71"

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
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

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
git-tree-sha1 = "7c2c19b5a26e601634bf718490b89d59685f122e"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.1"
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
git-tree-sha1 = "1fa950ebc3e37eccd51c6a8fe1f92f7d86263522"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.7+0"

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
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
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

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
registries = "General"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "4ac548adcad90c1d5d677af13568a748af4c952b"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.7"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
registries = "General"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

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

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
registries = "General"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

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

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
registries = "General"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "7a58e45171b63ed4782f2d36fdee8713a469e6e0"
registries = "General"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.2+0"

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
git-tree-sha1 = "5bad39456d9f0166184fce2248783dd9862645c1"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.0"
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

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "64bbbb7d1499297751b536dd39c58b20750ab1db"
registries = "General"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.5.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "4d777f73c46b46b8b5276206059cf8a195499314"
registries = "General"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.27"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f8eb8f7ba13ea75083531647fc8faeda8d541f07"
registries = "General"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.27+0"

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

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

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
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

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

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
registries = "General"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "633b5a34494e711f694ccbc88a6e00102f10238c"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.9.0"

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

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "df7566479bd64f20bd16b09960145e70160ffb3b"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.12"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

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
version = "1.11.103+0"

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
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.1"

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

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
registries = "General"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

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
git-tree-sha1 = "1bcebd887dd33f1108210b3954049b1bb8af0e7a"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.15+0"

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

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
registries = "General"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "f20e945b895d2009c6c28d8bbf40a5cd846f7c2f"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "83bd514e8ff16b5858ac54c53fa0bcf6002a3b00"
registries = "General"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.7"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "144895f6166994730ee7ff8113b981fc360638f1"
registries = "General"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "159d253ab126d5b29230cf53521899bea4ef4648"
registries = "General"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
registries = "General"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
registries = "General"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
registries = "General"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
registries = "General"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

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
git-tree-sha1 = "4db094d5e079abbda658acfe1c4d098430417717"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.8"

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

[[deps.Showoff]]
deps = ["Dates"]
git-tree-sha1 = "8238217340ad0aaabe11afe39c1098b5bc9f4c8e"
registries = "General"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.1.1"

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
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

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
git-tree-sha1 = "39e70e0ab5d7f89833a62ab7c79df15d4fc417c1"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.22"
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
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
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

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
registries = "General"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
registries = "General"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
registries = "General"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

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

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
registries = "General"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
registries = "General"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

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

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
registries = "General"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

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

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "dcb316b3ce0941f195537dda56bea4517fcd3ff5"
registries = "General"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.4+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
registries = "General"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
registries = "General"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

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

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
registries = "General"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
registries = "General"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
registries = "General"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
registries = "General"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
registries = "General"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
registries = "General"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
registries = "General"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
registries = "General"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "2e59214e017a55cb87474a00fa76035c82ac0e17"
registries = "General"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.47.0+2"

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

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
registries = "General"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
registries = "General"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

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

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
registries = "General"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "28e57478e8a160d346a19c28b3fffb9273bcc9c2"
registries = "General"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.134+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
registries = "General"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
registries = "General"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
registries = "General"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

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
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
registries = "General"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

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

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
registries = "General"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═e5f6a7b8-cccc-4a1b-8c2d-000000000001
# ╠═e5f6a7b8-cccc-4a1b-8c2d-000000000002
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000003
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000004
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000005
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000006
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000007
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000008
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000009
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000000a
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000000b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000000c
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000000d
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000000e
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000000f
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000010
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000011
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000012
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000013
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000014
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000015
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000016
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000017
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000018
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000019
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000001a
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000001b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000001c
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000001d
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000001e
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000001f
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000020
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000021
# ╟─e5f6a7b8-dddd-4a1b-8c2d-000000000001
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000022
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000023
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000024
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000025
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000026
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000027
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000028
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000029
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000002a
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000002b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000002c
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000002d
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000002e
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000002f
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000030
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000031
# ╟─e5f6a7b8-dddd-4a1b-8c2d-000000000002
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000032
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000033
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000034
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000035
# ╟─e5f6a7b8-dddd-4a1b-8c2d-000000000003
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000036
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000037
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000038
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000039
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000003a
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000003b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000003c
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000003d
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000003e
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000003f
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000040
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000041
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000042
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000043
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000044
# ╟─e5f6a7b8-dddd-4a1b-8c2d-000000000004
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000045
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000046
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000047
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000048
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000049
# ╟─e5f6a7b8-0000-4a1b-8c2d-000000000001
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000004a
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000004b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000004c
# ╟─e5f6a7b8-0000-4a1b-8c2d-000000000002
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000004d
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000004e
# ╟─e73a4fc2-f9fd-4fbd-ab6c-4b98464e1b6f
# ╟─522b9b6b-6401-4b75-844b-9c8db64a0ffd
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000004f
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000050
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000051
# ╟─bcbbe045-49c7-437e-a612-5438943056fc
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000052
# ╟─e5f6a7b8-0000-4a1b-8c2d-000000000003
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000053
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000054
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000055
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000056
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000057
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000058
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000059
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000001
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000002
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000003
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000004
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000005a
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000005
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000006
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000007
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000008
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000009
# ╠═e5f6a7b8-ffff-4a1b-8c2d-00000000000a
# ╟─e5f6a7b8-ffff-4a1b-8c2d-00000000000b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000005b
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000005c
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000005d
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000005e
# ╟─e5f6a7b8-ffff-4a1b-8c2d-00000000000c
# ╠═e5f6a7b8-ffff-4a1b-8c2d-00000000000d
# ╟─e5f6a7b8-ffff-4a1b-8c2d-00000000000e
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000005f
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000060
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000061
# ╟─e5f6a7b8-ffff-4a1b-8c2d-00000000000f
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000010
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000011
# ╟─c36e1dac-583a-4ace-b6d6-854d7e60f034
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000012
# ╠═3ab565b6-8975-49e3-b6fe-de57ebc518d0
# ╠═e5f6a7b8-ffff-4a1b-8c2d-000000000014
# ╠═d0243450-c789-4523-b239-eab2dff062ad
# ╠═c67999bc-c433-4430-87a5-06bd59202136
# ╠═99db9aa9-58c2-478e-b55c-b2c2c48e4176
# ╟─e5f6a7b8-ffff-4a1b-8c2d-000000000015
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000062
# ╠═ca9731d4-fb74-4e3c-a9a8-e5cab1bdf254
# ╠═20786201-ff40-462a-87d5-421b75787a0e
# ╠═f606ece1-f7a4-4fd7-a4e9-be28cf2c3be3
# ╠═88d558cf-2e3e-444e-8bf2-0901bf0a77f0
# ╠═fc934079-d098-4bab-b648-e3312c0902f8
# ╟─19adc32b-6e30-4349-a635-77757a8efd70
# ╠═cbb452d2-fe7a-40b6-bc31-eb36ab32b121
# ╟─e1de6880-500f-4386-a599-5731b3cfadd1
# ╠═85a8dbf3-cdc7-43e8-8b5b-aae5b2d0d080
# ╟─6dbfe708-ef86-4295-b6f7-1e99c8edbca9
# ╠═d2859751-1875-4e25-b383-1386195bd61c
# ╟─d38355e2-5150-40a3-99ed-0988ed3f1013
# ╠═38ac5ce1-82f5-4e7d-903b-20e9e01950bc
# ╟─c5d20dfe-eb1f-41ce-8cc3-dd877c32aa32
# ╠═e18b7fd5-1937-4df9-bff7-672690541421
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000063
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000064
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000065
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000066
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000067
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000068
# ╟─e5f6a7b8-cccc-4a1b-8c2d-000000000069
# ╟─e5f6a7b8-cccc-4a1b-8c2d-00000000006a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
