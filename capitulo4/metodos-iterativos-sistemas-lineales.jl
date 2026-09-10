### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000001
using LinearAlgebra, PlutoUI, Luxor

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000002
TableOfContents()

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000003
md"""
# Capítulo 4: Métodos iterativos estacionarios para sistemas lineales
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000004
md"""
En este capítulo describiremos algunos métodos iterativos para sistemas lineales. Los métodos *directos* como el método de eliminación de Gauss y sus variantes, como la descomposición LU, el método de Cholesky y adaptaciones a matrices de tipo banda, etc. son los métodos de elección para muchos problemas. Por favor no dejes que nada de lo dicho en este capítulo te haga olvidar esto.

Hay situaciones que requieren un tratamiento diferente al de los métodos directos, y aquí algunos inconvenientes de los métodos directos:

- El método de eliminación de Gauss o descomposición LU, en general, produce un *rellenado*: las matrices ``L`` y ``U`` pueden tener elementos no nulos en ubicaciones donde la matriz original ``A`` tiene ceros. Si la cantidad de relleno es significativa, aplicar un método directo puede resultar muy costoso. Esto ocurre a menudo en casos en que la matriz es de banda y rala (*sparse*) dentro de la banda.

- A veces no necesitamos resolver el sistema exactamente, sino que es suficiente una solución aproximada. Por ejemplo en el capítulo que sigue, veremos métodos para sistemas no lineales de ecuaciones en los cuales cada iteración involucra resolver un sistema lineal. Frecuentemente es suficiente resolver el sistema lineal dentro de cada iteración sólo aproximadamente. Se busca en estos casos un balance entre el costo de hallar la solución y el error de aproximación obtenido. Los métodos directos no pueden producir este balance, dado que llegan a la solución recién al completar todas las operaciones que los definen; los cálculos intermedios no dan lugar a soluciones aproximadas.

- A veces tenemos una muy buena aproximación de nuestra solución y queremos encontrar una solución más aproximada a bajo costo. Esto ocurre por ejemplo en ecuaciones diferenciales dependientes del tiempo. A menudo la solución del instante de tiempo anterior es cercana a la solución del próximo instante de tiempo, y es definitivamente provechoso, utilizarla como una aproximación inicial y tratar de mejorarla.

- A veces solamente tenemos una función o rutina que nos da el resultado de hacer el producto matriz vector y no disponemos de la matriz completa, o la tendríamos con un procedimiento muy costoso. Por ejemplo en aplicaciones del procesamiento digital de señales, a menudo el caso es que sólo tenemos las entradas y las salidas, o sea cuál es la salida, dada una entrada determinada, sin conocer la transformación de manera explícita. Muchos de los métodos iterativos que veremos, solamente necesitan acceso a una función o rutina ``x\to Ax``, que permita dado ``x`` conocer el resultado ``Ax``, sin necesidad de conocer la matriz completa.

En este capítulo veremos dos tipos de métodos iterativos para la resolución de sistemas lineales. Los primeros se denominan *estacionarios* y consisten en una iteración de la forma ``x_{k+1}=Mx_k+c``, con ``M`` y ``c`` fijos. Los últimos son no estacionarios y están basados en *métodos de descenso*, que sirven para matrices simétricas y definidas positivas (sdp).

En este capítulo usaremos la notación y los resultados básicos de Álgebra Lineal repasados en el Capítulo 2: el producto escalar y el producto tensorial, la escritura de los productos matriciales por filas y por columnas, y la caracterización de las matrices no singulares.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000005
md"""
## Normas de vectores y matrices. Número de condición
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000006
md"""
### Normas de vectores
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000007
md"""
Antes de estudiar los métodos iterativos necesitamos definir *convergencia* en ``\mathbb{R}^n``. Lo hacemos definiendo *normas* de vectores y matrices. Recordemos que si ``X`` es un espacio vectorial, una norma es una función ``\|\cdot\|:X\to\mathbb{R}`` que satisface

(i) ``\|x\|\geq0``, ``\forall x\in X``, y ``\|x\|=0`` sii ``x=0``.

(ii) ``\|\alpha x\|=|\alpha|\|x\|``, ``\forall\alpha\in\mathbb{R}``, ``\forall x\in X``. (Homogeneidad)

(iii) ``\|x+y\|\leq\|x\|+\|y\|``, ``\forall x,y\in X``. (Desigualdad triangular)

También recordamos que dos normas ``\|\cdot\|`` y ``\|\!|\cdot|\!\|`` en un espacio vectorial ``X`` se dicen equivalentes si existe ``c>0`` tal que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000008
md"""
```math
\frac{1}{c}\|x\| \leq \|\!|x|\!\| \leq c\|x\|, \qquad \forall x\in X.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000009
md"""
Si el espacio es de dimensión finita, entonces todas las normas son equivalentes. Esto puede demostrarse viendo que cualquier norma es equivalente a la *norma uno*, que, dada una base ``\{v_1,v_2,\dots,v_n\}`` de ``X`` se define como
"""

# ╔═╡ 2f9d20ff-9e4b-439a-9c00-fc8565d4bf61
@warn "Preguntar si no es necesaria la transitividad"

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000000a
md"""
```math
\|v\|_1 = \sum_{i=1}^n|x_i|, \qquad\text{si } v=\sum_{i=1}^nx_iv_i.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000000b
md"""
Las normas más utilizadas en ``\mathbb{R}^n`` son las *normas p*, ``1\leq p\leq\infty``, que se definen a continuación para ``x=(x_1,x_2,\dots,x_n)^T``:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000000c
md"""
```math
\|x\|_p = \left(\sum_{i=1}^n|x_i|^p\right)^{1/p}, \qquad\text{si } 1\leq p<\infty,
```
```math
\|x\|_\infty = \max_{1\leq i\leq n}|x_i|, \qquad\text{si } p=\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000000d
md"""
Las más usuales entre estas son:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000000e
md"""
```math
\|x\|_1 = \sum_{i=1}^n|x_i|, \qquad \|x\|_2 = \left(\sum_{i=1}^n|x_i|^2\right)^{1/2}, \qquad \|x\|_\infty = \max_{1\leq i\leq n}|x_i|.
```
"""

# ╔═╡ f0437696-4357-42f8-89ea-6e3a254a0c2a
md"""
Se deja como ejercicio verificar que las normas ``\|\cdot\|_1`` y ``\|\cdot\|_\infty`` satisfacen los axiomas de norma. Para la norma ``\|\cdot\|_2``, la verificación de la desigualdad triangular no es inmediata; para ello es conveniente introducir el siguiente concepto.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000000f
md"""
**Productos escalares.**

Un *producto escalar* (o *producto interno*) en un espacio vectorial real ``X`` es una función ``(\cdot,\cdot):X\times X\to\mathbb{R}`` que satisface:

(i) ``(x,y)=(y,x)``, ``\forall x,y\in X``. (Simetría)

(ii) ``(\alpha x+\beta y,z)=\alpha(x,z)+\beta(y,z)``, ``\forall\alpha,\beta\in\mathbb{R}``, ``\forall x,y,z\in X``. (Bilinealidad)

(iii) ``(x,x)\geq0``, ``\forall x\in X``, y ``(x,x)=0`` sii ``x=0``. (Definida positiva)

Definiendo ``\|x\|=(x,x)^{1/2}`` se obtiene una norma gracias a la siguiente desigualdad fundamental.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000010
md"""
**Teorema 4.1** (Desigualdad de Cauchy-Schwarz). *Sea ``(\cdot,\cdot)`` un producto escalar en ``X`` y sea ``\|\cdot\|`` la norma que genera. Entonces*

```math
|(x,y)| \leq \|x\|\|y\|, \qquad \forall x,y\in X.
```

(o, ``(x, y)^2 \leq (x, y)(y, y)``)

*Demostración.* Si ``y=0`` la desigualdad es trivial. Si ``y\neq0``, para todo ``t\in\mathbb{R}`` resulta
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000011
md"""
```math
0 \leq \|x-ty\|^2 = (x-ty,x-ty) = \|x\|^2-2t(x,y)+t^2\|y\|^2.
```
"""

# ╔═╡ cf19f3da-e87a-4624-a896-4b6276409318
md"""
O, de otra manera, si tomamos ``z=x-ty``, entonces ``(z, z)\geq 0`` y, además

```math
\begin{align*}
	0 &\leq (z, z) = (x - ty, x - ty) &\text{(Ax. III)}\\
		&= (x, x - ty) -t (y, x - ty) &\text{(Ax. II)}\\
		&= (x, x) -t(x, y) -t\big[(y, x) -t(y, y)\big] &\text{(Ax. II)}\\
		&= (x, x) -t(x, y) - t(y, x) + t^2 (y, y)\\
		&= (x, x) -t(x, y) - t(x, y) + t^2 (y, y) &\text{(Ax. I)}\\
		&= (x, x) -2t(x, y) + t^2 (y, y)\\
%%%		&= \big[(x, x) -t(y, y)\big]^2
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000012
md"""
Eligiendo ``t=(x,y)/\|y\|^2`` o ``t = (x, y)/(y,y)`` obtenemos
"""

# ╔═╡ 75fb8468-25e0-4506-bf65-1c2d93d58f84
md"""
```math
\begin{align*}
	0 &\leq (x, x) -2t(x, y) + t^2 (y, y)\\
		&= (x, x) - 2\big[(x, y) / (y, y)\big] (x,y) + \big[(x, y) / (y, y)\big]^2(y, y)\\
		&= (x,x) - 2(x,y)^2/(y,y) + (x,y)^2(y,y)/(y,y)^2\\
		&= (x,x) - 2(x,y)^2/(y,y) + (x,y)^2/(y,y)\\
		&= (x,x) - \frac{(x, y)^2}{(y,y)}
\end{align*}
```
"""

# ╔═╡ 28d73e21-98c3-40c1-a2c2-e4f78a605f5a
md"""
De donde 

```math
\begin{align*}
	0 &\leq (x, y) - \frac{(x, y)^2}{(y, y)}\\
	\frac{(x, y)^2}{(y, y)} &\leq (x, y)\\
	(x, y)^2 &\leq (x, x) (y, y)
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000013
md"""
o, como aparece en el libro,
```math
0 \leq \|x\|^2 - \frac{(x,y)^2}{\|y\|^2},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000014
md"""
de donde se sigue la desigualdad. 

**Ejemplo 4.2.** En ``\mathbb{R}^n``, el *producto escalar estándar* es ``(x,y)=x^Ty=\sum_{i=1}^nx_iy_i``, que genera la norma ``\|x\|_2``. La desigualdad de Cauchy-Schwarz toma aquí la forma
"""

# ╔═╡ dfcb451e-7d09-4a8a-8f60-a44569572824
@warn "Aca va con x o cdot?"

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000015
md"""
```math
|x^Ty| \leq \|x\|_2\|y\|_2.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000016
md"""
Más en general, dada ``A\in\mathbb{R}^{n\times n}`` simétrica y definida positiva (sdp), la función
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000017
md"""
```math
(x,y)_A = x^TAy
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000018
md"""
es un producto escalar en ``\mathbb{R}^n`` (EJERCICIO), que genera la *norma-A*
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000019
md"""
```math
\|x\|_A = (x,x)_A^{1/2} = (x^TAx)^{1/2}.
```
"""

# ╔═╡ 5466d779-1be8-40ac-9126-60fc3161784d
md"""
**Demostración**

1. ``\|x\|_A \geq 0``, pues al ser ``A`` sdp, se tiene que ``(x^TAx) \geq 0``, y por lo tanto, ``(x^TAx)^{1/2} \geq 0``.


2. ``\|\alpha x\|_A = \big((\alpha x^T)A(\alpha x)\big)^{1/2}=(\alpha^2x^TAx)^{1/2}=|\alpha|(x^TAx)^{1/2}``.


3.
```math
\begin{align*}
\|x+y\|_A &= \big((x+y)^{T}A(x+y)\big)^{1/2}&\text{(Por definición)}\\
	&= \big((x^T+y^T)A(x+y)\big)^{1/2}\\
	&= \big((x^TA + y^TA)(x+y)\big)^{1/2}\\
	&= \big(x^TAx + x^TAy + y^TAx + y^TAy\big)^{1/2}\\
	&= \big(x^TAx + 2x^TAy + y^TAy\big)^{1/2} & \text{(Por simetría de $A$)}\\
	&= \big(\big[(x^TAx)^{1/2} + (y^TAy)^{1/2}\big]^2\big)^{1/2}\\
	&= (x^TAx)^{1/2} + (y^TAy)^{1/2}\\
	&= \|x\|_A + \|y\|_A &\text{(Por definición)}
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000001a
md"""
En particular, ``\|x\|_I=\|x\|_2``.

**Definición 4.3** (Convergencia en ``\mathbb{R}^n``). Dada una sucesión ``\{x_k\}_{k\in\mathbb{N}}\subset\mathbb{R}^n``. Se dice que ``x_k`` tiende a ``x\in\mathbb{R}^n`` si ``\|x_k-x\|\to0`` cuando ``k\to\infty``, para alguna norma ``\|\cdot\|`` de ``\mathbb{R}^n``.

**Observación 4.4.** Puede demostrarse que todas las normas en ``\mathbb{R}^n`` son equivalentes (hecho que escapa al interés de este curso). Como consecuencia, si una sucesión converge en alguna norma, entonces converge en todas las normas. En particular ``\|x-x_k\|_\infty\to0`` y luego cada componente ``x_k^i\to x^i``, para ``i=1,2,\dots,n``, cuando ``k\to\infty``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000001b
md"""
### Normas de matrices
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000001c
md"""
En esta sección hablaremos sobre las normas en el espacio de las matrices. Prestaremos especial atención a las llamadas *normas inducidas*, que se definen como sigue.

**Definición 4.5.** Dada una norma ``\|\cdot\|`` del espacio de vectores ``\mathbb{R}^n``, se define la norma inducida por esta como
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000001d
md"""
```math
\|A\| = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\|Ax\|}{\|x\|}, \qquad \forall A\in\mathbb{R}^{n\times n}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000001e
md"""
La norma inducida se denota con el mismo símbolo de la norma que la induce, así, por ejemplo
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000001f
md"""
```math
\|A\|_1 = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\|Ax\|_1}{\|x\|_1} = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\sum_{i=1}^n|(Ax)_i|}{\sum_{i=1}^n|x_i|},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000001f1
md"""
```math
\|A\|_2 = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\|Ax\|_2}{\|x\|_2} = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\big(\sum_{i=1}^n|(Ax)_i|^2\big)^{1/2}}{\big(\sum_{i=1}^n|x_i|^2\big)^{1/2}},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000001f2
md"""
```math
\|A\|_\infty = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\|Ax\|_\infty}{\|x\|_\infty} = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\max_{1\leq i\leq n}|(Ax)_i|}{\max_{1\leq i\leq n}|x_i|}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000020
md"""
**En MATLAB/OCTAVE.** La función `norm` calcula diferentes normas de vectores y matrices. Leer su documentación.

**Observación 4.6.** Cualquiera sea la norma ``\|\cdot\|`` de ``\mathbb{R}^n``, la norma inducida de la matriz identidad es 1. En efecto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000021
md"""
```math
\|I\| = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\|Ix\|}{\|x\|} = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}} \frac{\|x\|}{\|x\|} = 1.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000022
md"""
**Observación 4.7.** Utilizando la *homogeneidad* de la norma puede demostrarse que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000023
md"""
```math
\|A\| = \max_{\substack{x\in\mathbb{R}^n\\0<\|x\|\leq1}} \frac{\|Ax\|}{\|x\|} = \max_{\substack{x\in\mathbb{R}^n\\0<\|x\|\leq1}} \|Ax\| = \max_{\substack{x\in\mathbb{R}^n\\\|x\|=1}} \|Ax\|,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000024
md"""
para toda ``A\in\mathbb{R}^{n\times n}`` y toda norma ``\|\cdot\|`` de ``\mathbb{R}^n``. EJERCICIO.

Un par de propiedades importantes de las normas inducidas se enuncian en el siguiente teorema.

**Teorema 4.8** (Consistencia). *Si ``\|\cdot\|`` es la norma de ``\mathbb{R}^{n\times n}`` inducida por la norma ``\|\cdot\|`` de ``\mathbb{R}^n``, entonces*

*(i) ``\|Ax\|\leq\|A\|\|x\|``, ``\forall A\in\mathbb{R}^{n\times n}``, ``\forall x\in\mathbb{R}^n``;*

*(ii) ``\|AB\|\leq\|A\|\|B\|``, ``\forall A,B\in\mathbb{R}^{n\times n}``.*

*Demostración.* EJERCICIO. 

No es difícil verificar (EJERCICIO) que, para ``A\in\mathbb{R}^{n\times n}``, resulta
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000025
md"""
```math
\|A\|_1 = \max_{1\leq j\leq n}\sum_{i=1}^n|a_{ij}| = \text{el máximo de las sumas por columnas},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000251
md"""
```math
\|A\|_\infty = \max_{1\leq i\leq n}\sum_{j=1}^n|a_{ij}| = \text{el máximo de las sumas por filas},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000252
md"""
```math
\|A\|_2 = \max\{\sqrt\lambda : \lambda \text{ autovalor de } A^TA\} = \text{máximo } \textit{valor singular} \text{ de } A.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000026
md"""
En el caso particular de ``A`` simétrica, se tiene que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000027
md"""
```math
\|A\|_2 = \max\{|\lambda| : \lambda\text{ autovalor de } A\}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000028
md"""
Esta última afirmación es consecuencia de que si ``A`` es simétrica, entonces ``A=U\Lambda\bar U^T``, con ``\Lambda`` una matriz diagonal que tiene los autovalores de ``A`` en la diagonal, y ``U`` es *unitaria*, es decir, ``\bar U^TU=I``. Aquí y en adelante ``\bar U`` denota la matriz conjugada de ``U``, es decir ``(\bar U)_{ij}=\bar u_{ij}``, donde ``\bar u_{ij}`` es el complejo conjugado de ``u_{ij}=(U)_{ij}``.

Otra norma usual de ``\mathbb{R}^{n\times n}``, que no es inducida, es la *norma de Frobenius*, definida por
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000029
md"""
```math
\|A\|_F = \left(\sum_{i,j=1}^na_{ij}^2\right)^{1/2} = \sqrt{\text{tr}(A^TA)}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000002a
md"""
A pesar de no ser inducida, esta norma también cumple que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000002b
md"""
```math
\|AB\|_F \leq \|A\|_F\|B\|_F.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000002c
md"""
Además, si ``\{\sigma_i\}_{i=1}^r`` son los valores singulares de ``A``, ``\|A\|_F=\big(\sum_{i=1}^r\sigma_i^2\big)^{1/2}``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000002d
md"""
### Número de condición y error relativo
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000002e
md"""
Dada ``A\in\mathbb{R}^{n\times n}`` y ``b\in\mathbb{R}^n``, consideremos el problema de hallar ``x\in\mathbb{R}^n`` tal que ``Ax=b``. Supongamos que ``\tilde x\in\mathbb{R}^n`` es una solución aproximada del sistema. Considerando una norma dada ``\|\cdot\|`` de ``\mathbb{R}^n``, definimos:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000002f
md"""
```math
\text{error} = x-\tilde x,
```
```math
\text{error relativo} = \frac{\|x-\tilde x\|}{\|x\|},
```
```math
\text{residuo} = A\tilde x-b = \tilde b-b, \quad\text{si } \tilde b=A\tilde x,
```
```math
\text{residuo relativo} = \frac{\|\tilde b-b\|}{\|b\|}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000030
md"""
Queremos hallar una estimación para el error relativo ``\frac{\|x-\tilde x\|}{\|x\|}`` en términos del residuo relativo ``\frac{\|\tilde b-b\|}{\|b\|}``. Observemos lo siguiente
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000031
md"""
```math
\|\tilde x-x\| = \|A^{-1}\tilde b-A^{-1}b\| \leq \|A^{-1}\|\|\tilde b-b\|,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000032
md"""
Además,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000033
md"""
```math
\|Ax\|=\|b\| \implies \|b\|\leq\|A\|\|x\|.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000034
md"""
Luego,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000035
md"""
```math
\frac{\|\tilde x-x\|}{\|x\|} \leq \|A\|\|A^{-1}\|\frac{\|\tilde b-b\|}{\|b\|}. \qquad(4.1)
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000036
md"""
Análogamente
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000037
md"""
```math
\|\tilde b-b\| = \|A\tilde x-Ax\| \leq \|A\|\|\tilde x-x\|,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000038
md"""
y
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000039
md"""
```math
\|x\| = \|A^{-1}b\| \leq \|A^{-1}\|\|b\|.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000003a
md"""
Por lo tanto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000003b
md"""
```math
\frac{1}{\|A\|\|A^{-1}\|}\frac{\|\tilde b-b\|}{\|b\|} \leq \frac{\|\tilde x-x\|}{\|x\|}. \qquad(4.2)
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000003c
md"""
De (4.1)–(4.2) obtenemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000003d
md"""
```math
\frac{1}{\|A\|\|A^{-1}\|}\frac{\|\tilde b-b\|}{\|b\|} \leq \frac{\|\tilde x-x\|}{\|x\|} \leq \|A\|\|A^{-1}\|\frac{\|\tilde b-b\|}{\|b\|}. \qquad(4.3)
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000003e
md"""
**Definición 4.9.** Dada una norma inducida ``\|\cdot\|`` en el espacio de las matrices de ``n\times n``, se define el número de condición para esa norma como
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000003f
md"""
```math
\kappa(A) = \|A\|\|A^{-1}\|, \qquad A\in\mathbb{R}^{n\times n}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000040
md"""
Resumimos lo demostrado hasta ahora en un Lema:

**Lema 4.10.** *Si ``x\in\mathbb{R}^n\setminus\{0\}`` es la solución de ``Ax=b``, entonces, cualquiera sea ``\tilde x\in\mathbb{R}^n``, se cumple que*

```math
\frac{1}{\kappa(A)}\frac{\|\tilde b-b\|}{\|b\|} \leq \frac{\|\tilde x-x\|}{\|x\|} \leq \kappa(A)\frac{\|\tilde b-b\|}{\|b\|},
```

*para ``\tilde b=A\tilde x``.*

**En MATLAB/OCTAVE.** La función `cond` calcula el número de condición para diferentes normas y matrices. Leer su documentación, y la de la función `rcond`.

**Observación 4.11.** El número de condición siempre es mayor o igual a uno:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000041
md"""
```math
\kappa(A) = \|A\|\|A^{-1}\| \geq \|AA^{-1}\| = \|I\| = 1.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000042
md"""
Cuando ``\kappa(A)`` es pequeño, o cercano a uno, ``\kappa(A)\cong1`` resulta
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000043
md"""
```math
\frac{\|\tilde x-x\|}{\|x\|} \cong \frac{\|\tilde b-b\|}{\|b\|}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000044
md"""
Es decir, pequeñas perturbaciones en el lado derecho de la ecuación, producen pequeños cambios en las soluciones. En este caso decimos que el problema ``Ax=b`` está *bien condicionado*.

**Observación 4.12.** Dada la solución ``x`` de ``Ax=b``, si ``\tilde x`` es una aproximación de ``x``, entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000045
md"""
```math
\underbrace{\frac{\|\tilde x-x\|}{\|x\|}}_{\text{desconocido}} \leq \kappa(A)\underbrace{\frac{\|A\tilde x-b\|}{\|b\|}}_{\text{calculable}}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000046
md"""
El lado derecho de la desigualdad anterior se utiliza en la implementación de criterios de parada en algunos métodos iterativos.

**Observación 4.13.** Observemos que multiplicar una matriz por una constante arbitraria no nula, no afecta el número de condición:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000047
md"""
```math
\kappa(\alpha A) = \|\alpha A\|\|(\alpha A)^{-1}\| = |\alpha|\|A\|\frac{1}{|\alpha|}\|A^{-1}\| = \kappa(A).
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000048
md"""
Existe una estrecha relación entre el número de condición y los autovalores cuando se considera la norma ``\|\cdot\|_2`` para matrices simétricas, la misma está dada en el siguiente lema.

**Lema 4.14.** *Si ``A\in\mathbb{R}^{n\times n}`` es simétrica, entonces ``\kappa_2(A)=\frac{\lambda_{\max}}{\lambda_{\min}}``, donde ``\lambda_{\max}`` y ``\lambda_{\min}`` denotan el autovalor de módulo máximo y mínimo de ``A``, respectivamente.*

*Demostración.* Sea ``A\in\mathbb{R}^{n\times n}`` una matriz simétrica, entonces ``\|A\|_2=\rho(A)=\max\{|\lambda| : \lambda\in\sigma(A)\}=:\lambda_{\max}``. Como ``A`` es simétrica, ``A^{-1}`` resulta simétrica y entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000049
md"""
```math
\|A^{-1}\|_2 = \rho(A^{-1}) = \max\{|\lambda| : \lambda\in\sigma(A^{-1})\}
```
```math
= \max\{|\lambda|^{-1} : \lambda^{-1}\in\sigma(A)\} = \max\{|\lambda|^{-1} : \lambda\in\sigma(A)\}
```
```math
= \frac{1}{\min\{|\lambda| : \lambda\in\sigma(A)\}} =: \frac{1}{\lambda_{\min}}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000004a
md"""
Por lo tanto ``\kappa_2(A)=\|A\|_2\|A^{-1}\|_2=\lambda_{\max}\frac{1}{\lambda_{\min}}``. 
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000004b
md"""
### Número de condición para matrices rectangulares
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000004c
md"""
La definición de número de condición dada anteriormente se restringe a matrices cuadradas e invertibles. Existe una extensión natural al caso rectangular, motivada por el estudio de *sistemas sobredeterminados*. Dado ``A\in\mathbb{R}^{n\times m}`` con ``n\geq m`` y ``b\in\mathbb{R}^n``, consideramos el *problema de cuadrados mínimos*
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000004d
md"""
```math
\min_{x\in\mathbb{R}^m} \|Ax-b\|_2.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000004e
md"""
Cuando ``A`` tiene rango completo por columnas, este problema tiene solución única, caracterizada por las ecuaciones normales ``A^TAx=A^Tb``.

**Pseudoinversa y valores singulares.**

La *descomposición en valores singulares* (SVD) de ``A\in\mathbb{R}^{n\times m}`` con ``n\geq m`` es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000004f
md"""
```math
A = U\Sigma V^T,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000050
md"""
con ``U\in\mathbb{R}^{n\times n}`` y ``V\in\mathbb{R}^{m\times m}`` ortogonales, y
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000051
md"""
```math
\Sigma = \begin{bmatrix}\text{diag}(\sigma_1,\sigma_2,\dots,\sigma_r,0,\dots,0)\\0\end{bmatrix} \in \mathbb{R}^{n\times m},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000052
md"""
donde ``r=\text{rango}(A)`` y ``\sigma_1\geq\sigma_2\geq\dots\geq\sigma_r>0`` son los *valores singulares* de ``A``, definidos como las raíces cuadradas positivas de los autovalores de ``A^TA``. Notemos que ``\|A\|_2=\sigma_1``.

Si ``A`` tiene rango completo por columnas (``r=m``), su *pseudoinversa de Moore-Penrose* es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000053
md"""
```math
A^+ = (A^TA)^{-1}A^T \in \mathbb{R}^{m\times n},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000054
md"""
la cual satisface ``A^+A=I_m`` y puede escribirse en términos de la SVD como ``A^+=V\Sigma^+U^T``, donde ``\Sigma^+=\begin{bmatrix}\text{diag}(\sigma_1^{-1},\dots,\sigma_m^{-1}) & 0\end{bmatrix}``. En particular, ``\|A^+\|_2=1/\sigma_m``.

**Definición 4.15.** Sea ``A\in\mathbb{R}^{n\times m}`` con ``n\geq m`` y rango completo por columnas. El *número de condición* de ``A`` es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000055
md"""
```math
\kappa_2(A) = \|A\|_2\|A^+\|_2 = \frac{\sigma_1}{\sigma_m},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000056
md"""
donde ``\sigma_1\geq\sigma_2\geq\dots\geq\sigma_m>0`` son los valores singulares de ``A``.

**Observación 4.16.** Esta definición coincide con la anterior cuando ``A`` es cuadrada e invertible: en ese caso ``A^+=A^{-1}``, y los valores singulares de ``A`` son las raíces cuadradas de los autovalores de ``A^TA``.

La siguiente relación, simple pero de enorme importancia práctica, vincula el número de condición de ``A`` con el de la matriz de coeficientes de las ecuaciones normales.

**Lema 4.17.** *Sea ``A\in\mathbb{R}^{n\times m}`` con ``n\geq m`` y rango completo por columnas. Entonces*

```math
\kappa_2(A^TA) = \kappa_2(A)^2.
```

*Demostración.* Sea ``A=U\Sigma V^T`` la SVD de ``A``. Entonces ``A^TA=V(\Sigma^T\Sigma)V^T``, con ``\Sigma^T\Sigma=\text{diag}(\sigma_1^2,\dots,\sigma_m^2)``. Como ``A^TA`` es simétrica definida positiva, sus autovalores coinciden con sus valores singulares, y por el Lema anterior
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000057
md"""
```math
\kappa_2(A^TA) = \left(\frac{\sigma_1}{\sigma_m}\right)^2 = \kappa_2(A)^2. 
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000058
md"""
**Observación 4.18.** Este resultado tiene una consecuencia práctica importante. Si ``\kappa_2(A)`` es grande, entonces ``\kappa_2(A^TA)=\kappa_2(A)^2`` es enormemente mayor. Resolver el problema de cuadrados mínimos mediante las ecuaciones normales ``(A^TA)x=A^Tb`` obliga a trabajar con una matriz cuyo número de condición es el cuadrado del de ``A``, lo que puede provocar serios problemas de estabilidad numérica. Por esta razón, en la práctica se prefiere utilizar la factorización QR de ``A`` (véase el Capítulo 3), que tiene mejor comportamiento numérico.

**Observación 4.19.** Análogamente al Lema 4.10, el número de condición controla la sensibilidad de la solución de cuadrados mínimos a perturbaciones en ``b``. Si ``x^*=A^+b`` y ``\tilde x=A^+(b+\delta b)``, entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000059
md"""
```math
\|\tilde x-x^*\|_2 \leq \|A^+\|_2\|\delta b\|_2,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000005a
md"""
y en el caso de un sistema compatible (``Ax^*=b``):
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000005b
md"""
```math
\frac{\|\tilde x-x^*\|_2}{\|x^*\|_2} \leq \kappa_2(A)\frac{\|\delta b\|_2}{\|b\|_2}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000005c
md"""
**En MATLAB/OCTAVE.** La función `cond(A)` calcula el número de condición también para matrices rectangulares (mediante la SVD). La función `svd(A)` devuelve la descomposición ``A=U\Sigma V^T``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000005d
md"""
### Espectro y radio espectral
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000005e
md"""
Dada una matriz ``A\in\mathbb{R}^{n\times n}`` se define su *espectro* ``\sigma(A)`` como el conjunto de autovalores de ``A``, así
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000005f
md"""
```math
\sigma(A) = \{\lambda\in\mathbb{C} : \lambda \text{ es autovalor de } A\}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000060
md"""
El *radio espectral* ``\rho(A)`` de ``A`` es el máximo módulo de sus autovalores, más precisamente
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000061
md"""
```math
\rho(A) = \max\{|\lambda| : \lambda\in\sigma(A)\}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000062
md"""
Por lo dicho anteriormente,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000063
md"""
```math
\text{si } A \text{ es simétrica} \implies \rho(A)=\|A\|_2.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000064
md"""
El radio espectral juega un papel fundamental en el estudio de métodos iterativos para sistemas lineales, principalmente debido al siguiente teorema:

**Teorema 4.20** (Caracterización del radio espectral). *Sea ``A\in\mathbb{R}^{n\times n}``, entonces*

```math
\rho(A) = \inf\{\|A\| : \|\cdot\| \text{ es una norma inducida}\}.
```

Para demostrar este teorema utilizaremos el siguiente teorema, cuya demostración postergamos hasta el final de la sección.

**Teorema 4.21** (Descomposición de Schur). *Sea ``A\in\mathbb{C}^{n\times n}``. Sean ``\lambda_1,\lambda_2,\dots\lambda_n`` los autovalores de ``A`` contados con su multiplicidad algebraica, es decir que son todos los ceros de ``p(\lambda):=\det(\lambda I-A)``. Entonces existe una matriz unitaria ``U`` (``\bar U^TU=I``) tal que ``UA\bar U^T=T``, donde ``T`` es triangular superior y ``\text{diag}(T)=(\lambda_1,\lambda_2,\dots,\lambda_n)``.*

**Observación 4.22.** Si los autovalores de ``A`` son reales, entonces existe una matriz real ``U`` ortogonal (``U^TU=I``) tal que ``UAU^T=T``, con ``T`` triangular superior (real) que tiene los autovalores de ``A`` en la diagonal. Además, vale la pena enfatizar que esta descomposición es cierta aunque ``A`` no sea diagonalizable.

A continuación utilizaremos el Teorema 4.21 de descomposición de Schur para demostrar el Teorema 4.20 de caracterización del radio espectral.

*Demostración del Teorema 4.20.* Es fácil ver que ``\rho(A)\leq\|A\|`` para toda norma inducida ``\|\cdot\|`` y cualquier matriz cuadrada ``A`` (EJERCICIO). Lo más difícil es ver que dada ``A\in\mathbb{R}^{n\times n}``, para todo ``\varepsilon>0`` existe una norma inducida ``\|\cdot\|`` tal que ``\|A\|<\rho(A)+\varepsilon``.

Sean entonces ``A\in\mathbb{R}^{n\times n}`` y ``\varepsilon>0`` dados. Por el Teorema 4.21 de descomposición de Schur existe ``U`` unitaria tal que ``UA\bar U^T=T`` con ``T`` triangular superior tal que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000065
md"""
```math
\text{diag}(T) = (\lambda_1,\lambda_2,\dots,\lambda_n), \qquad\text{y}\qquad \lambda_i, \quad i=1,2,\dots,n \quad\text{son los autovalores de } A.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000066
md"""
Entonces ``T`` es de la forma
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000067
md"""
```math
T = \begin{bmatrix}\lambda_1&t_{12}&t_{13}&\dots&t_{1n}\\0&\lambda_2&t_{23}&\dots&t_{2n}\\\vdots&\ddots&\ddots&\ddots&\vdots\\0&\dots&0&\lambda_{n-1}&t_{(n-1)n}\\0&\dots&0&0&\lambda_n\end{bmatrix}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000068
md"""
Para ``z\in\mathbb{R}``, ``z>0``, definimos
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000069
md"""
```math
D_z = \text{diag}(z,z^2,\dots,z^n) = \begin{bmatrix}z&0&0&\dots&0\\0&z^2&0&\dots&0\\\vdots&\ddots&\ddots&\ddots&\vdots\\0&0&\dots&z^{n-1}&0\\0&0&\dots&0&z^n\end{bmatrix}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000006a
md"""
Así,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000006b
md"""
```math
D_zUA\bar U^T = D_zT = \begin{bmatrix}z\lambda_1&zt_{12}&zt_{13}&\dots&zt_{1n}\\0&z^2\lambda_2&z^2t_{23}&\dots&z^2t_{2n}\\\vdots&\ddots&\ddots&\ddots&\vdots\\0&\dots&0&z^{n-1}\lambda_{n-1}&z^{n-1}t_{(n-1)n}\\0&\dots&0&0&z^n\lambda_n\end{bmatrix},
```
"""

# ╔═╡ f6ca739b-2bd9-46f1-88c0-4f277da6f81b
@warn "Multiplicar al reves. TDz"

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000006c
md"""
y
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000006d
md"""
```math
D_zUA\bar U^TD_z^{-1} = D_zTD_z^{-1} = \begin{bmatrix}\lambda_1&z^{-1}t_{12}&z^{-2}t_{13}&\dots&z^{-(n-1)}t_{1n}\\0&\lambda_2&z^{-1}t_{23}&\dots&z^{-(n-2)}t_{2n}\\\vdots&\ddots&\ddots&\ddots&\vdots\\0&\dots&0&\lambda_{n-1}&z^{-1}t_{(n-1)n}\\0&\dots&0&0&\lambda_n\end{bmatrix}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000006e
md"""
Elegimos ahora ``z`` suficientemente grande de manera que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000006f
md"""
```math
\max_{1\leq j\leq n}\sum_{i=1}^{j-1}z^{i-j}|t_{ij}| < \varepsilon,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000070
md"""
y entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000071
md"""
```math
\|D_zUA\bar U^TD_z^{-1}\|_1 = \max_{1\leq j\leq n}\left(|\lambda_j|+\sum_{i=1}^{j-1}z^{i-j}|t_{ij}|\right) < \rho(A)+\varepsilon.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000072
md"""
Faltaría ver que ``\|D_zUA\bar U^TD_z^{-1}\|_1`` es una norma inducida de ``A``. Definimos entonces la norma ``\|x\|=\|D_zUx\|_1`` para ``x\in\mathbb{R}^n`` (que resulta una norma, verificar), y observemos que
"""

# ╔═╡ 7e76323b-8d03-4454-ada1-7e175d930a43
@warn "La norma vector es D-1zU"

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000073
md"""
```math
\|A\| = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}\frac{\|Ax\|}{\|x\|} = \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}\frac{\|D_zUAx\|_1}{\underbrace{\|D_zUx\|_1}_{y}} = \max_{\substack{y\in\mathbb{R}^n\\y\neq0}}\frac{\|D_zUA(D_zU)^{-1}y\|_1}{\|y\|_1} = \|D_zUA\bar U^TD_z^{-1}\|_1. 
```
"""

# ╔═╡ afa824e3-52c4-41b0-a94a-53a793ccaeaf
@warn "Ese paso intermedio es por despejar y"

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000074
md"""
La importancia del *radio espectral* queda plasmada en el siguiente teorema.

**Teorema 4.23.** *Sea ``A\in\mathbb{R}^{n\times n}``, entonces las siguientes afirmaciones son equivalentes:*

*(1) ``\rho(A)<1``.*

*(2) ``\|A\|<1`` para alguna norma inducida.*

*(3) ``\lim_{k\to\infty}\|A^k\|=0``, para toda norma inducida ``\|\cdot\|``.*

*(4) ``\lim_{k\to\infty}A^kx=0``, para todo vector ``x\in\mathbb{R}^n``.*

*Demostración.* (1) ``\Rightarrow`` (2) es inmediato por el Teorema 4.20.

(2) ``\Rightarrow`` (3). Supongamos que ``\|A\|_*<1`` para alguna norma inducida ``\|\cdot\|_*``. Sea ``\|\cdot\|`` una norma inducida cualquiera. Entonces, como todas las normas son equivalentes, por ser el espacio de las matrices de dimensión finita, existe ``c>0`` tal que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000075
md"""
```math
\|A^k\| \leq c\|A^k\|_*, \qquad \forall k\in\mathbb{N}_0.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000076
md"""
Por el Teorema 4.8 de consistencia, resulta que ``\|A^k\|_*\leq\|A\|_*^k`` y luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000077
md"""
```math
\|A^k\| \leq c\|A^k\|_* \leq c\|A\|_*^k \to 0, \quad\text{cuando } k\to\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000078
md"""
(3) ``\Rightarrow`` (4). Sea ``x\in\mathbb{R}^n`` y sea ``\|\cdot\|`` una norma en ``\mathbb{R}^n`` y también su norma inducida en ``\mathbb{R}^{n\times n}``. Entonces, por el Teorema 4.8 de consistencia
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000079
md"""
```math
\|A^kx-0\| = \|A^kx\| \leq \|A^k\|\|x\| \to 0, \quad\text{cuando } k\to\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000007a
md"""
(4) ``\Rightarrow`` (1). Sea ``\lambda\in\sigma(A)`` y sea ``x`` un autovector asociado a ``\lambda``, luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000007b
md"""
```math
Ax=\lambda x, \qquad\text{y}\qquad A^kx=\lambda^kx.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000007c
md"""
Así,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000007d
md"""
```math
\|A^kx\| = |\lambda|^k\|x\|
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000007e
md"""
para cualquier norma ``\|\cdot\|`` de ``\mathbb{R}^n``. Por (4), ``\|A^kx\|\to0`` y por ende ``|\lambda|^k\|x\|\to0``. Como ``x`` es un autovector, es no nulo, y por lo tanto ``|\lambda|^k\to0`` cuando ``k\to\infty``. Por lo tanto ``|\lambda|<1``. Como ``\lambda\in\sigma(A)`` era arbitrario, resulta ``\rho(A)=\max_{\lambda\in\sigma(A)}|\lambda|<1``. 

**Teorema 4.24.** *Sea ``A\in\mathbb{R}^{n\times n}``. Entonces, cualquiera sea la norma inducida ``\|\cdot\|``, resulta*

```math
\rho(A) = \lim_{k\to\infty}\|A^k\|^{1/k}.
```

*Demostración.* Veamos primero que ``\rho(A)\leq\|A^k\|^{1/k}``. Sea ``\lambda\in\sigma(A)`` y sea ``x`` un autovector asociado a ``\lambda``. Entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000007f
md"""
```math
Ax=\lambda x \implies A^kx=\lambda^kx \implies \|A^kx\|=|\lambda|^k\|x\|
```
```math
\implies |\lambda|^k = \frac{\|A^kx\|}{\|x\|} \leq \|A^k\| \implies |\lambda|\leq\|A^k\|^{1/k} \implies \rho(A)\leq\|A^k\|^{1/k}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000080
md"""
Veamos ahora que el límite de ``\|A^k\|^{1/k}`` es ``\rho(A)``. Sea ``\varepsilon>0`` y definimos
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000081
md"""
```math
\tilde A = \frac{1}{\rho(A)+\varepsilon}A \implies \rho(\tilde A) = \frac{1}{\rho(A)+\varepsilon}\rho(A) = \frac{\rho(A)}{\rho(A)+\varepsilon} < 1.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000082
md"""
Por el Teorema anterior, ``\lim_{k\to\infty}\|\tilde A^k\|=0`` y por lo tanto, existe ``k_0\in\mathbb{N}`` tal que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000083
md"""
```math
\|\tilde A^k\| = \frac{\|A^k\|}{(\rho(A)+\varepsilon)^k} < 1, \qquad \forall k\geq k_0.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000084
md"""
Es decir que dado ``\varepsilon>0``, existe ``k_0\in\mathbb{N}`` tal que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000085
md"""
```math
\rho(A) \leq \|A^k\|^{1/k} \leq \rho(A)+\varepsilon, \qquad \forall k\geq k_0.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000086
md"""
O sea, ``\lim_{k\to\infty}\|A^k\|^{1/k}=\rho(A)``. 

**Demostración del Teorema 4.21 de descomposición de Schur.**

*Demostración del Teorema 4.21.* Haremos la demostración por inducción sobre el orden ``n`` de la matriz. Si ``n=1``, entonces ``A=[a_{11}]``, ``\lambda_1=a_{11}`` y el teorema es válido con ``T=[a_{11}]`` y ``U=[1]``.

Supongamos que la afirmación del teorema es válida para ``n-1`` con ``n>1`` y demostrémosla para ``n``.

Sean ``\lambda_1,\lambda_2,\dots,\lambda_n`` los autovalores de ``A`` y sea ``x`` un autovector correspondiente a ``\lambda_1`` con ``\|x\|_2^2=\bar x^Tx=1`` y ``x_1\geq0``; notemos que si ``\lambda_1`` es complejo, entonces algunas componentes de ``x`` pueden ser complejas, en ese caso, si ``x_1`` no fuera un número real no-negativo, dividimos todas las componentes de ``x`` por ``x_1`` y luego lo normalizamos. Utilizaremos la transformación de Householder
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000087
md"""
```math
H = I-u\bar u^T, \qquad\text{con } u=\frac{1}{\sqrt{1+x_1}}(x+e_1),
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000088
md"""
que, como vimos antes, satisface
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000089
md"""
```math
\bar H^T=H=H^{-1}, \text{ (H es unitaria)} \qquad\text{y}\qquad Hx=-e_1, \qquad He_1=-x.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000008a
md"""
Entonces ``HA\bar H^T=HAH`` y
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000008b
md"""
```math
\text{col}_1(HA\bar H^T) = HAHe_1 = -HAx = -\lambda_1Hx = \lambda_1e_1.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000008c
md"""
Es decir,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000008d
md"""
```math
HA\bar H^T = \begin{bmatrix}\lambda_1&v^T\\0&\tilde A\end{bmatrix}, \quad\text{con } v\in\mathbb{R}^{n-1} \text{ y } \tilde A\in\mathbb{R}^{(n-1)\times(n-1)}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000008e
md"""
Además,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000008f
md"""
```math
\lambda\in\sigma(A) \iff \lambda\in\sigma(HA\bar H^T) \iff \lambda=\lambda_1 \text{ ó } \lambda\in\sigma(\tilde A),
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000090
md"""
pues ``\det(\lambda I-HA\bar H^T)=(\lambda-\lambda_1)\det(\lambda I-\tilde A)``. Por lo tanto, ``\sigma(\tilde A)=\{\lambda_2,\dots,\lambda_n\}``.

Por la hipótesis inductiva, existe ``\tilde U`` unitaria (de orden ``n-1``) tal que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000091
md"""
```math
\tilde U\tilde A\bar{\tilde U}^T = \tilde T, \quad\text{con } \tilde T \text{ triangular superior y } \text{diag}(\tilde T)=(\lambda_2,\dots,\lambda_n).
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000092
md"""
Es decir, ``\bar{\tilde U}^T\tilde T\tilde U=\tilde A``, y luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000093
md"""
```math
HA\bar H^T = \begin{bmatrix}\lambda_1&v^T\\0&\bar{\tilde U}^T\tilde T\tilde U\end{bmatrix} = \begin{bmatrix}1&0^T\\0&\bar{\tilde U}^T\end{bmatrix}\begin{bmatrix}\lambda_1&v^T\bar{\tilde U}^T\\0&\tilde T\end{bmatrix}\begin{bmatrix}1&0^T\\0&\tilde U\end{bmatrix},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000094
md"""
y finalmente
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000095
md"""
```math
\underbrace{\begin{bmatrix}1&0^T\\0&\tilde U\end{bmatrix}}_{U} HA\bar H^T \underbrace{\begin{bmatrix}1&0^T\\0&\bar{\tilde U}^T\end{bmatrix}}_{\bar U^T} = \underbrace{\begin{bmatrix}\lambda_1&v^T\bar{\tilde U}^T\\0&\tilde T\end{bmatrix}}_{T}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000096
md"""
## Métodos iterativos estacionarios
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000097
md"""
Queremos diseñar métodos iterativos para resolver el problema
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000098
md"""
```math
Ax=b,
```
"""

# ╔═╡ e9ed9b75-1f05-4c90-bb5c-94c5e9325820
md"""
con ``A\in\mathbb{R}^{n\times n}`` no singular y ``b\in\mathbb{R}^n`` dados.
"""

# ╔═╡ 18d8ba73-e4aa-4a16-a83e-60c4243466a2
md"""
```math
\begin{align*}
	x_{k+1} &= Mx_{k} + c\\
	x &= Mx + c\\
	x - x_{k+1} &= M(x - x_k)\\
	\|x - x_{k+1}\| &= \|M(x - x_k)\|\\
	\|x - x_{k+1}\| &\leq \|M\|\|x - x_k\|\\
	\|e_k\| &\leq \|M\|^k \|e_0\|
\end{align*}
```

Si existe una norma tal que ``\|M\| < 1 `` entonces converge.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000099
md"""
**Método de Richardson.**

Una primera idea surge de observar que, si ``\alpha\neq0``,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000009a
md"""
```math
Ax=b \iff x=x-\alpha(Ax-b) \iff x=(I-\alpha A)x+\alpha b,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000009b
md"""
para algún número ``\alpha\in\mathbb{R}``. A partir de esta igualdad podemos definir la iteración
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000009c
md"""
```math
x_{k+1} = \underbrace{(I-\alpha A)}_{\text{matriz de iteración}}x_k+\alpha b, \qquad k=0,1,2,\dots,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000009d
md"""
con ``x_0\in\mathbb{R}^n`` dado. Esta iteración se conoce como *método de Richardson*.

**Método de Jacobi.**

Otra opción surge de considerar que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000009e
md"""
```math
Ax=b \iff \begin{aligned}a_{11}x_1 &= b_1-(a_{12}x_2+\dots+a_{1n}x_n)\\a_{22}x_2 &= b_2-(a_{21}x_1)-(a_{23}x_3+\dots+a_{2n}x_n)\\a_{33}x_3 &= b_3-(a_{31}x_1+a_{32}x_2)-(a_{34}x_4+\dots+a_{3n}x_n)\\&\ \ \vdots\\a_{nn}x_n &= b_n-(a_{n1}x_1+a_{n2}x_2+\dots+a_{n(n-1)}x_{n-1}).\end{aligned}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000009f
md"""
Esto sugiere, por ejemplo, el siguiente método iterativo. Dado ``x_0``, repetir,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a0
md"""
```math
\begin{aligned}
x_{k+1}^1 &= \big[b_1-(a_{12}x_k^2+\dots+a_{1n}x_k^n)\big]/a_{11}\\
x_{k+1}^2 &= \big[b_2-(a_{21}x_k^1)-(a_{23}x_k^3+\dots+a_{2n}x_k^n)\big]/a_{22}\\
x_{k+1}^3 &= \big[b_3-(a_{31}x_k^1+a_{32}x_k^2)-(a_{34}x_k^4+\dots+a_{3n}x_k^n)\big]/a_{33}, \qquad k=0,1,2,\dots\\
&\ \ \vdots\\
x_{k+1}^n &= \big[b_n-(a_{n1}x_k^1+a_{n2}x_k^2+\dots+a_{n(n-1)}x_k^{n-1})\big]/a_{nn}
\end{aligned}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a1
md"""
Este método iterativo se conoce como *método de Jacobi*.

Si escribimos la matriz ``A=L+D+U`` con ``L`` triangular inferior con ceros en la diagonal, ``U`` triangular superior con ceros en la diagonal y ``D`` diagonal, entonces la iteración del método de Jacobi puede escribirse en forma matricial como
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a2
md"""
```math
x_{k+1} = D^{-1}[b-(L+U)x_k] = \underbrace{-D^{-1}(L+U)}_{\text{matriz de iteración}}x_k+D^{-1}b,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a3
md"""
notando que por ser ``D`` diagonal, ``D^{-1}`` no es costosa de calcular.

**Método de Gauss-Seidel.**

Si miramos el método de Jacobi, y pensamos en una implementación secuencial donde se calcula primero ``x_{k+1}^1``, luego ``x_{k+1}^2``, y así sucesivamente, uno podría pensar que si sabe que el método converge, entonces ``x_{k+1}^1`` es *mejor* que ``x_k^1``, y entonces proponer que para calcular ``x_{k+1}^2`` se utilice ``x_{k+1}^1`` en su lugar. Si la misma idea se aplica para el cálculo de ``x_{k+1}^3``, y los que siguen, resulta el método de Gauss-Seidel:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a4
md"""
```math
\begin{aligned}
x_{k+1}^1 &= \big[b_1-(a_{12}x_k^2+\dots+a_{1n}x_k^n)\big]/a_{11}\\
x_{k+1}^2 &= \big[b_2-(a_{21}x_{k+1}^1)-(a_{23}x_k^3+\dots+a_{2n}x_k^n)\big]/a_{22}\\
x_{k+1}^3 &= \big[b_3-(a_{31}x_{k+1}^1+a_{32}x_{k+1}^2)-(a_{34}x_k^4+\dots+a_{3n}x_k^n)\big]/a_{33}\\
&\ \ \vdots\\
x_{k+1}^n &= \big[b_n-(a_{n1}x_{k+1}^1+a_{n2}x_{k+1}^2+\dots+a_{n(n-1)}x_{k+1}^{n-1})\big]/a_{nn}
\end{aligned}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a5
md"""
Si descomponemos ``A=L+D+U`` como antes, esta iteración puede verse en forma matricial como
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a6
md"""
```math
x_{k+1} = \underbrace{-(L+D)^{-1}U}_{\text{matriz de iteración}}x_k+(L+D)^{-1}b,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a7
md"""
**Método de Sobrerrelajación Sucesiva (SOR).**

Dado un número ``\omega\in\mathbb{R}\setminus\{0\}``, consideramos la siguiente *partición* de ``A``:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a8
md"""
```math
A = \Big(\frac{1}{\omega}D+L\Big)+\Big(\frac{\omega-1}{\omega}D+U\Big),
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000a9
md"""
con ``A=L+D+U`` como antes. Esto nos conduce al método de *SOR* que consiste en
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000aa
md"""
```math
x_{k+1} = \Big(\frac{1}{\omega}D+L\Big)^{-1}\Big[b-\Big(\frac{\omega-1}{\omega}D+U\Big)x_k\Big]
```
```math
= -\underbrace{\Big(\frac{1}{\omega}D+L\Big)^{-1}\Big(\frac{\omega-1}{\omega}D+U\Big)}_{\text{matriz de iteración}}x_k + \Big(\frac{1}{\omega}D+L\Big)^{-1}b.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ab
md"""
Observemos que el caso ``\omega=1`` coincide con el método de Gauss-Seidel.

Los métodos descritos son todos métodos de la forma
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ac
md"""
```math
x_{k+1} = Mx_k+c,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ad
md"""
con ``M`` y ``c`` dados. En el siguiente teorema se establece una condición necesaria y suficiente sobre ``M`` para que tal iteración converja a partir de cualquier ``x_0``.

**Teorema 4.25.** *Sea ``M\in\mathbb{R}^{n\times n}``. Entonces, dado ``c\in\mathbb{R}^n``, la iteración*

```math
x_{k+1} = Mx_k+c, \qquad k=0,1,\dots,
```

*converge para cualquier ``x_0`` inicial si y sólo si ``\rho(M)<1``. Además, si ``\rho(M)<1``, la iteración converge a la única solución de ``x=Mx+c``.*

*Demostración.* (``\Leftarrow``): Si ``\rho(M)<1`` entonces existe una norma inducida ``\|\cdot\|`` tal que ``\|M\|<1``. Veamos que ``I-M`` es invertible para concluir que existe única solución de ``x=Mx+c`` y luego demostrar que ``x_k`` converge a ``x``.

Observemos que como ``\|M\|<1`` entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ae
md"""
```math
\sum_{k=0}^N\|M^k\| \leq \sum_{k=0}^N\|M\|^k \leq \sum_{k=0}^\infty\|M\|^k = \frac{1}{1-\|M\|} < \infty,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000af
md"""
y entonces ``\sum_{k=0}^\infty M^k`` converge absolutamente. Definimos ``N=\sum_{k=0}^\infty M^k`` y vemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b0
md"""
```math
N(I-M) = N-NM = \sum_{k=0}^\infty M^k-\sum_{k=0}^\infty M^kM = \sum_{k=0}^\infty M^k-\sum_{k=1}^\infty M^k = I.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b1
md"""
Por lo tanto ``N=(I-M)^{-1}`` y ``I-M`` es invertible. Luego el sistema
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b2
md"""
```math
x=Mx+c \iff (I-M)x=c
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b3
md"""
tiene solución única que llamamos ``x``. Veamos que cualquiera sea ``x_0\in\mathbb{R}^n``, si definimos
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b4
md"""
```math
x_{k+1}=Mx_k+c, \qquad k=0,1,2,\dots,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b5
md"""
resulta que ``x_k\to x``. Restando la ecuación para ``x`` y la definición de ``x_k`` obtenemos la ecuación del error
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b6
md"""
```math
(x_{k+1}-x) = M(x_k-x) = M^2(x_{k-1}-x) = M^{k+1}(x_0-x).
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b7
md"""
Por lo tanto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b8
md"""
```math
\|x_{k+1}-x\| \leq \|M^{k+1}\|\|x_0-x\| \leq \|M\|^{k+1}\|x_0-x\| \longrightarrow 0 \quad\text{cuando } k\to\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000b9
md"""
(``\Rightarrow``): Si el método ``x_{k+1}=Mx_k+c`` converge a partir de cualquier ``x_0``, debe converger siempre a ``x`` solución de ``x=Mx+c``. Si converge para cualquier ``c\in\mathbb{R}^n`` dado, entonces el sistema ``(I-M)x=c`` tiene solución para todo ``c`` por lo que, fijado el ``c``, la solución es única. Resumiendo, dado ``c\in\mathbb{R}^n``, la iteración ``x_{k+1}=Mx_k+c`` converge, a partir de cualquier ``x_0``, al mismo ``x``, solución única de ``x=Mx+c``. Luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ba
md"""
```math
x_k-x = M^k(x_0-x) \longrightarrow 0, \quad\text{cuando } k\to\infty,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000bb
md"""
cualquiera sea ``x_0``. Esto implica que ``M^ky\to0`` para todo ``y\in\mathbb{R}^n`` y por el Teorema 4.23 resulta ``\rho(M)<1``. 

**Observación 4.26.** Notemos que por el Teorema 4.23, la condición ``\rho(M)<1`` es equivalente a que exista una norma inducida ``\|\cdot\|`` tal que ``\|M\|<1``.

**Observación 4.27.** La norma ``\|M\|<1`` da también una noción de la *reducción del error* en cada iteración, pues
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000bc
md"""
```math
\|x-x_{k+1}\| \leq \|M\|\|x-x_k\|, \qquad k=0,1,2,\dots.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000bd
md"""
Así, la convergencia en la norma ``\|\cdot\|`` será más rápida cuanto más pequeño sea ``\|M\|``.

Para demostrar, bajo ciertas condiciones, que los métodos de Jacobi, Gauss-Seidel, SOR y Richardson convergen, veremos que la matriz de iteración correspondiente tiene alguna norma inducida menor a uno.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000be
md"""
### Convergencia del método de Jacobi
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000bf
md"""
La convergencia del método de Jacobi está garantizada para matrices *estrictamente diagonalmente dominantes*.

**Definición 4.28.** Una matriz ``A\in\mathbb{R}^{n\times n}`` se dice *estrictamente diagonalmente dominante* (edd) si
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c0
md"""
```math
|a_{ii}| > \sum_{j\neq i}|a_{ij}|, \qquad\text{para todo } i=1,2,\dots,n.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c1
md"""
**Teorema 4.29.** *Si ``A`` es una matriz estrictamente diagonalmente dominante (edd), entonces ``A`` es invertible, y dado ``b\in\mathbb{R}^n`` el método de Jacobi*

```math
x_{k+1} = -D^{-1}(L+U)x_k+\underbrace{D^{-1}b}_{c}, \qquad k=0,1,\dots,
```

*converge a la única solución de ``Ax=b`` a partir de cualquier iteración inicial ``x_0``.*

*Además,*

```math
\|x-x_{k+1}\|_\infty \leq \mu_J\|x-x_k\|_\infty, \qquad k=0,1,2,\dots,
```

*con ``\mu_J=\max_{1\leq i\leq n}\frac{\sum_{j\neq i}|a_{ij}|}{|a_{ii}|}``, que satisface ``\mu_J<1`` por ser ``A`` edd.*

*Demostración.* La demostración consiste en verificar que la norma inducida por la norma infinito de la matriz de iteración del método de Jacobi es menor que uno, usando que ``A`` es edd.

Recordemos que la iteración de Jacobi es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c2
md"""
```math
x_{k+1} = D^{-1}[b-(L+U)x_k] = \underbrace{-D^{-1}(L+U)}_{M_J}x_k+\underbrace{D^{-1}b}_{c},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c3
md"""
y observemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c4
md"""
```math
M_J = -D^{-1}(L+U) = -\begin{bmatrix}0&\frac{a_{12}}{a_{11}}&\frac{a_{13}}{a_{11}}&\dots&\frac{a_{1n}}{a_{11}}\\\frac{a_{21}}{a_{22}}&0&\frac{a_{23}}{a_{22}}&\dots&\frac{a_{2n}}{a_{22}}\\\vdots&\ddots&\ddots&\ddots&\vdots\\\vdots&&\ddots&\ddots&\\\frac{a_{n1}}{a_{nn}}&\dots&\dots&\frac{a_{n(n-1)}}{a_{nn}}&0\end{bmatrix},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c5
md"""
o más precisamente
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c6
md"""
```math
(M_J)_{ii}=0, \qquad (M_J)_{ij}=-\frac{a_{ij}}{a_{ii}}, \quad\text{si } i\neq j.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c7
md"""
La norma infinito de esta matriz es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c8
md"""
```math
\|M_J\|_\infty = \max_{1\leq i\leq n}\sum_{j=1}^n|(M_J)_{ij}| = \max_{1\leq i\leq n}\sum_{j\neq i}\frac{|a_{ij}|}{|a_{ii}|} = \max_{1\leq i\leq n}\frac{\sum_{j\neq i}|a_{ij}|}{|a_{ii}|} < 1,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000c9
md"""
por ser ``A`` edd. Esto implica que el problema ``x=M_Jx+c`` tiene solución única para cada ``c\in\mathbb{R}^n``, y como es equivalente a ``Ax=b``, se tiene que ``A`` es no singular. 

Sin usar la maquinaria abstracta de normas de matrices también se puede demostrar la convergencia de este método, de la siguiente manera más elemental.

*Demostración alternativa del Teorema 4.29.* Sea ``A\in\mathbb{R}^{n\times n}`` edd, sea ``b\in\mathbb{R}^n`` dado y sea ``x_0\in\mathbb{R}^n`` cualquier aproximación inicial. Recordemos que por definición, la iteración de Jacobi es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ca
md"""
```math
x_{k+1}^i = \frac{b_i}{a_{ii}}-\sum_{j\neq i}\frac{a_{ij}}{a_{ii}}x_k^j, \qquad i=1,2,\dots,n,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000cb
md"""
y claramente la solución de ``Ax=b`` satisface
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000cc
md"""
```math
x_i = \frac{b_i}{a_{ii}}-\sum_{j\neq i}\frac{a_{ij}}{a_{ii}}x_j, \qquad i=1,2,\dots,n.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000cd
md"""
Luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ce
md"""
```math
\underbrace{x_{k+1}^i-x_i}_{e_{k+1}^i} = -\sum_{j\neq i}\frac{a_{ij}}{a_{ii}}\underbrace{(x_k^j-x_j)}_{e_k^j}, \qquad i=1,2,\dots,n,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000cf
md"""
y por lo tanto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d0
md"""
```math
|e_{k+1}^i| \leq \sum_{j\neq i}\frac{|a_{ij}|}{|a_{ii}|}|e_k^j| \leq \sum_{j\neq i}\frac{|a_{ij}|}{|a_{ii}|}\|e_k\|_\infty \leq \left(\max_{1\leq i\leq n}\sum_{j\neq i}\frac{|a_{ij}|}{|a_{ii}|}\right)\|e_k\|_\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d1
md"""
Si definimos ``\mu_J=\max_{1\leq i\leq n}\frac{\sum_{j\neq i}|a_{ij}|}{|a_{ii}|}`` resulta
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d2
md"""
```math
\|e_{k+1}\|_\infty \leq \mu_J\|e_k\|_\infty,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d3
md"""
y ``\mu_J<1`` por ser ``A`` edd. 
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d4
md"""
### Convergencia del método de Gauss-Seidel
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d5
md"""
Para el método de Gauss-Seidel tenemos el siguiente teorema de convergencia.

**Teorema 4.30.** *Si ``A`` es estrictamente diagonalmente dominante (edd), entonces el método de Gauss-Seidel converge a partir de cualquier iteración inicial ``x_0`` a la única solución ``x=A^{-1}b`` de ``Ax=b``.*

*Además, existe ``\mu_{GS}<1`` tal que*

```math
\|x-x_{k+1}\|_\infty \leq \mu_{GS}\|x-x_k\|_\infty, \qquad k=0,1,2,\dots.
```

*Más aún, ``\mu_{GS}=\max_{1\leq i\leq n}\frac{\beta_i}{1-\alpha_i}``, con ``\alpha_1=0``, ``\alpha_i=\sum_{j<i}\frac{|a_{ij}|}{|a_{ii}|}``, ``i=2,3,\dots,n``, ``\beta_i=\sum_{j>i}\frac{|a_{ij}|}{|a_{ii}|}``, ``i=1,2,\dots,n-1``, ``\beta_n=0``.*

*Demostración.* Recordemos que de la definición del método de Gauss-Seidel
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d6
md"""
```math
x_{k+1} = (L+D)^{-1}[b-Ux_k] = \underbrace{-(L+D)^{-1}U}_{M_{GS}}x_k+\underbrace{(L+D)^{-1}b}_{c},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d7
md"""
o bien
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d8
md"""
```math
(L+D)x_{k+1} = b-Ux_k.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000d9
md"""
Además la solución exacta satisface ``(L+D)x=b-Ux``. Luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000da
md"""
```math
D(x_{k+1}-x) = -L(x_{k+1}-x)-U(x_k-x),
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000db
md"""
y definiendo ``e_k=x_k-x`` tenemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000dc
md"""
```math
e_{k+1} = -D^{-1}Le_{k+1}-D^{-1}Ue_k.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000dd
md"""
Para ``i=1,2,\dots,n``
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000de
md"""
```math
e_{k+1}^i = -\text{fila}_i(D^{-1}L)e_{k+1}-\text{fila}_i(D^{-1}U)e_k,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000df
md"""
y entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e0
md"""
```math
|e_{k+1}^i| \leq \|\text{fila}_i(D^{-1}L)\|_1\|e_{k+1}\|_\infty+\|\text{fila}_i(D^{-1}U)\|_1\|e_k\|_\infty = \alpha_i\|e_{k+1}\|_\infty+\beta_i\|e_k\|_\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e1
md"""
Puesto que ``A`` es edd, resulta que ``\alpha_i+\beta_i<1`` por lo que ``\alpha_i<1`` y ``\beta_i<1`` (ya que son ambos no negativos). Dado ``k``, sea ``i_0`` tal que ``\|e_{k+1}\|_\infty=|e_{k+1}^{i_0}|``, entonces
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e2
md"""
```math
\|e_{k+1}\|_\infty(1-\alpha_{i_0}) \leq \beta_{i_0}\|e_k\|_\infty \implies \|e_{k+1}\|_\infty \leq \frac{\beta_{i_0}}{1-\alpha_{i_0}}\|e_k\|_\infty \leq \underbrace{\max_{1\leq i\leq n}\frac{\beta_i}{1-\alpha_i}}_{\mu_{GS}}\|e_k\|_\infty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e3
md"""
Para asegurar que el algoritmo converge, y finalizar la demostración del teorema, falta ver que ``\frac{\beta_i}{1-\alpha_i}<1``, ``i=1,2,\dots,n``. Veremos que ``\frac{\beta_i}{1-\alpha_i}<\alpha_i+\beta_i<1``. La última desigualdad vale por el teorema anterior de convergencia del método de Jacobi. Para ver la primera, observemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e4
md"""
```math
\alpha_i+\beta_i-\frac{\beta_i}{1-\alpha_i} = \frac{\alpha_i-\alpha_i^2+\beta_i-\alpha_i\beta_i-\beta_i}{1-\alpha_i} = \frac{\alpha_i}{1-\alpha_i}\underbrace{(1-\alpha_i-\beta_i)}_{>0} \geq 0. 
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e5
md"""
**Observación 4.31.** Del último paso de la demostración anterior concluimos que para matrices edd,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e6
md"""
```math
\mu_J = \max_{1\leq i\leq n}(\alpha_i+\beta_i) \geq \max_{1\leq i\leq n}\frac{\beta_i}{1-\alpha_i} = \mu_{GS},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e7
md"""
y por lo tanto, para este tipo de matrices, el método de Gauss-Seidel funciona mejor o igual al de Jacobi.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e8
md"""
### Convergencia del método SOR
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000e9
md"""
Recordemos que el método de sobrerrelajación sucesiva (SOR) es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ea
md"""
```math
x_{k+1} = -\underbrace{\Big(\frac{1}{\omega}D+L\Big)^{-1}\Big(\frac{\omega-1}{\omega}D+U\Big)}_{M_{SOR}^\omega}x_k+\Big(\frac{1}{\omega}D+L\Big)^{-1}b,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000eb
md"""
con ``x_0`` dado.

Para estudiar la convergencia de este método, vemos primero el siguiente teorema:

**Teorema 4.32.** *Sea ``A\in\mathbb{R}^{n\times n}`` sdp. Si ``A=A_1+A_2`` con ``A_1`` invertible y tal que ``A_1^T-A_2`` es sdp, entonces ``\rho(-A_1^{-1}A_2)<1``.*

*Demostración.* Puesto que ``A`` es sdp, sabemos que ``\|x\|_A=(x^TAx)^{1/2}`` es una norma en ``\mathbb{R}^n``, y denotemos también con ``\|\cdot\|_A`` a la norma inducida por ésta en ``\mathbb{R}^{n\times n}``. Veremos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ec
md"""
```math
\|A_1^{-1}A_2\|_A = \max_{\|x\|_A=1}\|A_1^{-1}A_2x\|_A < 1, \qquad (4.4)
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ed
md"""
y esto implicará la afirmación del teorema.

Sea ``x\in\mathbb{R}^n`` que realiza el máximo de (4.4). Entonces ``\|x\|_A=1`` y ``\|A_1^{-1}A_2\|_A=\|A_1^{-1}A_2x\|_A``.

Notemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ee
md"""
```math
A_1^{-1}A_2x = A_1^{-1}(A-A_1)x = \underbrace{A_1^{-1}Ax}_{y}-x,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ef
md"""
y por lo tanto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f0
md"""
```math
\|A_1^{-1}A_2x\|_A^2 = \|y-x\|_A^2 = (y-x)^TA(y-x) = y^TAy-x^TAy-y^TAx+x^TAx. \qquad (4.5)
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f1
md"""
Como ``A`` es simétrica, ``x^TAy=y^TAx`` y como ``y=A_1^{-1}Ax``, resulta que ``x=A^{-1}A_1y``, por lo que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f2
md"""
```math
x^TAy = y^TAx = y^TA(A^{-1}A_1y) = y^TA_1y = y^TA_1^Ty.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f3
md"""
Reemplazando esto último en (4.5), resulta
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f4
md"""
```math
\|A_1^{-1}A_2x\|_A^2 = y^TAy-y^T(A_1+A_1^T)y+x^TAx
```
```math
= x^TAx+y^T(\underbrace{A-A_1}_{A_2}-A_1^T)y = x^TAx-y^T\underbrace{(A_1^T-A_2)}_{=K\text{ sdp}}y
```
```math
= x^TAx-\underbrace{y^TKy}_{>0\text{ pues }y\neq0} < x^TAx = \|x\|_A^2 = 1.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f5
md"""
Finalmente, ``\rho(A_1^{-1}A_2)\leq\|A_1^{-1}A_2x\|_A<1``. 

**Teorema 4.33** (Convergencia del método SOR). *Si ``A`` es sdp, entonces el método SOR converge a partir de cualquier ``x_0\in\mathbb{R}^n``, si ``\omega\in(0,2)``.*

*Demostración.* Se deja como ejercicio. Se logra aplicando el teorema anterior. 

Si bien el teorema anterior da una condición suficiente para que SOR converja, el siguiente teorema dice que no es posible obtener un algoritmo convergente si se toma ``\omega\notin(0,2)``.

**Teorema 4.34.** *Dada una matriz invertible ``A\in\mathbb{R}^{n\times n}``, y ``b\in\mathbb{R}^n``, si ``\omega\notin(0,2)`` entonces la matriz de iteración del método SOR tiene radio espectral mayor o igual a uno.*

*Demostración.* Recordemos que la matriz de iteración de SOR es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f6
md"""
```math
M_{SOR}^\omega = -\Big(\frac{1}{\omega}D+L\Big)^{-1}\Big(\frac{\omega-1}{\omega}D+U\Big),
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f7
md"""
y su determinante es
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f8
md"""
```math
\det(M_{SOR}^\omega) = (-1)^n\frac{\det(\frac{\omega-1}{\omega}D+U)}{\det(\frac{1}{\omega}D+L)} = (-1)^n\frac{\det(\frac{\omega-1}{\omega}D)}{\det(\frac{1}{\omega}D)} = (-1)^n\Big(\frac{\omega-1}{\omega}\Big)^n\omega^n\frac{\det(D)}{\det(D)} = (1-\omega)^n.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000f9
md"""
Si ``\omega\notin(0,2)`` entonces ``|\omega-1|\geq1``. Pero además, como el determinante es el producto de los autovalores,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000fa
md"""
```math
\prod_{\lambda\in\sigma(M_{SOR}^\omega)}|\lambda| = |\det(M_{SOR}^\omega)| = |\omega-1|^n \geq 1,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000fb
md"""
por lo que algún autovalor de ``M_{SOR}^\omega`` debe tener módulo mayor o igual a uno, y luego
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000fc
md"""
```math
\rho(M_{SOR}^\omega) \geq 1. 
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000fd
md"""
### Convergencia del método de Richardson
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000fe
md"""
Recordemos que el método de Richardson se basa en la iteración
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-0000000000ff
md"""
```math
x_{k+1} = (I-\alpha A)x_k+\alpha b, \qquad k=0,1,2,\dots,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000100
md"""
con ``\alpha\neq0`` y ``x_0`` dados.

Consideremos ``A`` sdp, y estudiemos ``\rho(I-\alpha A)``. Como ``A`` es sdp, sus autovalores son todos positivos. Denotamos con ``\lambda_{\min}`` y ``\lambda_{\max}`` al menor y al mayor autovalor de ``A``, respectivamente. Además, los autovalores de ``I-\alpha A`` son ``1-\alpha\lambda`` para ``\lambda\in\sigma(A)``, más precisamente
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000101
md"""
```math
\sigma(I-\alpha A) = \{1-\alpha\lambda : \lambda\in\sigma(A)\} \subset [1-\alpha\lambda_{\max},1-\alpha\lambda_{\min}].
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000102
md"""
Por lo tanto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000103
md"""
```math
\rho(I-\alpha A) = \max\{|1-\alpha\lambda_{\max}|,|1-\alpha\lambda_{\min}|\}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000104
md"""
Así,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000105
md"""
```math
\begin{aligned}
\rho(I-\alpha A)<1 &\iff |1-\alpha\lambda_{\max}|<1 \text{ y } |1-\alpha\lambda_{\min}|<1\\
&\iff -1<1-\alpha\lambda_{\max}<1 \text{ y } -1<1-\alpha\lambda_{\min}<1\\
&\iff 0<\alpha<\frac{2}{\lambda_{\max}} \text{ y } 0<\alpha<\frac{2}{\lambda_{\min}}\\
&\iff 0<\alpha<\frac{2}{\lambda_{\max}},
\end{aligned}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000106
md"""
pues ``\lambda_{\min}\leq\lambda_{\max}``.

Concluimos entonces que

**Lema 4.35.** *El método de Richardson converge para ``A`` sdp, a partir de cualquier ``x_0\in\mathbb{R}^n`` si y solo si*

```math
0<\alpha<\frac{2}{\lambda_{\max}}.
```

Queremos ahora averiguar cuál es el valor de ``\alpha`` que minimiza ``\rho(I-\alpha A)``. Por lo visto anteriormente,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000107
md"""
```math
\alpha_{\text{opt}} = \text{argmin}\max\{|1-\alpha\lambda_{\min}|,|1-\alpha\lambda_{\max}|\},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000108
md"""
es decir, ``\alpha_{\text{opt}}`` es el valor de ``\alpha`` que minimiza la función ``\alpha\to\max\{|1-\alpha\lambda_{\min}|,|1-\alpha\lambda_{\max}|\}``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000109
md"""
**Figura 4.1** (descripción, no reproducida gráficamente). El gráfico muestra, en función de ``\alpha\geq0``, las dos rectas quebradas ``|1-\alpha\lambda_{\max}|`` (que vale 1 en ``\alpha=0``, decrece hasta anularse en ``\alpha=1/\lambda_{\max}``, y luego vuelve a crecer) y ``|1-\alpha\lambda_{\min}|`` (que vale 1 en ``\alpha=0``, decrece más lentamente hasta anularse en ``\alpha=1/\lambda_{\min}``, y luego vuelve a crecer, con ``1/\lambda_{\min}\geq1/\lambda_{\max}``). El máximo de ambas (la curva marcada en el gráfico) decrece con la primera recta y luego crece con la segunda; ``\alpha_{\text{opt}}`` es exactamente el punto donde ambas rectas se cruzan, es decir donde ``|1-\alpha\lambda_{\max}|=|1-\alpha\lambda_{\min}|`` con ambas expresiones de signos opuestos.

A partir del gráfico de la Figura 4.1 vemos que el mismo se alcanza cuando
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000010a
md"""
```math
\alpha\lambda_{\max}-1 = 1-\alpha\lambda_{\min} \iff \alpha = \frac{2}{\lambda_{\min}+\lambda_{\max}}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000010b
md"""
En ese caso, además
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000010c
md"""
```math
\rho(I-\alpha_{\text{opt}}A) = 1-\alpha_{\text{opt}}\lambda_{\min} = 1-\frac{2\lambda_{\min}}{\lambda_{\min}+\lambda_{\max}} = \frac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\max}+\lambda_{\min}} = \frac{\frac{\lambda_{\max}}{\lambda_{\min}}-1}{\frac{\lambda_{\max}}{\lambda_{\min}}+1},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000010d
md"""
y como ``A`` es sdp, ``\kappa_2(A)=\frac{\lambda_{\max}}{\lambda_{\min}}``. Hemos probado entonces el siguiente teorema:

**Teorema 4.36.** *Si ``A`` es sdp, y aplicamos el método de Richardson con ``\alpha_{opt}=\frac{2}{\lambda_{\min}+\lambda_{\max}}``, a partir de cualquier ``x_0`` obtenemos*

```math
\|x-x_{k+1}\|_2 \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x-x_k\|_2, \qquad k=0,1,2,\dots,
```

*donde ``\kappa_2(A)`` es el número de condición de ``A`` con respecto a la norma euclídea de vectores, que satisface*

```math
\kappa_2(A) = \|A\|_2\|A^{-1}\|_2 = \frac{\lambda_{\max}}{\lambda_{\min}}.
```

Es importante observar que el decrecimiento del error (en norma euclídea) es más lento cuando el número de condición de la matriz es más grande. El mismo decrecimiento del error ocurre cuando se mide el mismo en la norma generada por ``A``, más precisamente, se tiene el siguiente teorema.

**Teorema 4.37.** *Si ``A`` es sdp, y aplicamos el método de Richardson con ``\alpha_{opt}=\frac{2}{\lambda_{\min}+\lambda_{\max}}``, a partir de cualquier ``x_0`` obtenemos*

```math
\|x-x_{k+1}\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x-x_k\|_A, \qquad k=0,1,2,\dots,
```

*donde ``\|y\|_A=(y^TAy)^{1/2}``.*

*Demostración.* Para probar la afirmación bastaría ver que ``\|I-\alpha_{\text{opt}}A\|_A\leq\frac{\kappa_2(A)-1}{\kappa_2(A)+1}``, es decir
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000010e
md"""
```math
\max_{\|v\|_A=1}\|(I-\alpha_{\text{opt}}A)v\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000010f
md"""
Sea entonces ``v\in\mathbb{R}^n`` tal que ``\|v\|_A=1``. Observemos que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000110
md"""
```math
\|(I-\alpha_{\text{opt}}A)v\|_A^2 = v^T(I-\alpha_{\text{opt}}A)^TA(I-\alpha_{\text{opt}}A)v = v^T(I-\alpha_{\text{opt}}A)^TA^{1/2}A^{1/2}(I-\alpha_{\text{opt}}A)v,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000111
md"""
donde ``A^{1/2}`` es la *raíz cuadrada* de ``A``. Como ``A`` es sdp, existen ``Q\in\mathbb{R}^{n\times n}`` ortogonal y ``\Lambda\in\mathbb{R}^{n\times n}`` diagonal (con los autovalores de ``A`` en la diagonal) tales que ``A=Q\Lambda Q^T``. Como ``\Lambda`` es diagonal y sus entradas son no-negativas, se define ``\Lambda^{1/2}`` por ``(\Lambda^{1/2})_{ij}=(\Lambda_{ij})^{1/2}``, y cumple ``\Lambda^{1/2}\Lambda^{1/2}=\Lambda``. Entonces se define ``A^{1/2}=Q\Lambda^{1/2}Q^T`` que cumple
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000112
md"""
```math
A^{1/2}A^{1/2} = Q\Lambda^{1/2}Q^TQ\Lambda^{1/2}Q^T = Q\Lambda^{1/2}\Lambda^{1/2}Q^T = Q\Lambda Q^T = A.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000113
md"""
Además
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000114
md"""
```math
AA^{1/2} = Q\Lambda Q^TQ\Lambda^{1/2}Q^T = Q\Lambda\Lambda^{1/2}Q^T = Q\Lambda^{1/2}\Lambda Q^T = Q\Lambda^{1/2}Q^TQ\Lambda Q^T = A^{1/2}A,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000115
md"""
es decir que ``A`` conmuta con ``A^{1/2}``. Por lo tanto
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000116
md"""
```math
\|(I-\alpha_{\text{opt}}A)v\|_A^2 = \underbrace{v^TA^{1/2}}_{w^T}(I-\alpha_{\text{opt}}A)^T(I-\alpha_{\text{opt}}A)\underbrace{A^{1/2}v}_{w} = \|(I-\alpha_{\text{opt}}A)w\|_2^2.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000117
md"""
Ahora bien ``\|w\|_2^2=w^Tw=v^TA^{1/2}A^{1/2}v=v^TAv=\|v\|_A^2=1``. Finalmente, resulta que cualquiera sea ``v\in\mathbb{R}^n`` con ``\|v\|_A=1``,
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000118
md"""
```math
\|(I-\alpha_{\text{opt}}A)v\|_A \leq \max_{\|w\|_2=1}\|(I-\alpha_{\text{opt}}A)w\|_2 = \rho(I-\alpha_{\text{opt}}A) = \frac{\kappa_2(A)-1}{\kappa_2(A)+1}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000119
md"""
Por lo tanto, ``\|I-\alpha_{\text{opt}}A\|_A\leq\frac{\kappa_2(A)-1}{\kappa_2(A)+1}`` y la afirmación se sigue de la consistencia de la norma inducida.
"""

# ╔═╡ d4e5f6a7-9999-4a1b-8c2d-000000000001
md"""
**Resumen (no está en el apunte, viene de las diapositivas de la clase 7).** Condiciones de convergencia para los cuatro métodos vistos:

| Método | Condición suficiente | Condición necesaria y suficiente |
|---|---|---|
| Richardson | ``A`` sdp, ``0<\alpha<2/\lambda_{\max}`` | ``\rho(I-\alpha A)<1`` |
| Jacobi | ``A`` edd | ``\rho(M_J)<1`` |
| Gauss-Seidel | ``A`` edd (o ``A`` sdp) | ``\rho(M_{GS})<1`` |
| SOR | ``A`` sdp, ``\omega\in(0,2)`` | ``\rho(M_{SOR}^\omega)<1`` |

En todos los casos, la condición necesaria y suficiente es ``\rho(M)<1``; las condiciones suficientes (edd o sdp) son las que se usan en la práctica porque son más fáciles de verificar. Para matrices edd, además, ``\mu_{GS}\leq\mu_J<1`` (Observación 4.31): Gauss-Seidel no es peor que Jacobi. Que Gauss-Seidel también converja cuando ``A`` es sdp se sigue de tomar ``\omega=1`` en el Teorema 4.33 de convergencia de SOR (dado que ``\omega=1\in(0,2)``), recordando que el caso ``\omega=1`` de SOR es exactamente Gauss-Seidel.
"""

# ╔═╡ d4e5f6a7-9999-4a1b-8c2d-000000000002
md"""
**Código (no está en el apunte, viene de las diapositivas de la clase 7).** Implementación en MATLAB/OCTAVE de los métodos de Jacobi y de SOR (que incluye Gauss-Seidel tomando ``\omega=1``):

```octave
function x = jacobi(A, b, x0, tol, maxit)
  D  = diag(diag(A));
  LU = A - D;
  x  = x0;
  for k = 1:maxit
    x_new = D \ (b - LU*x);
    if norm(x_new - x) < tol, break; end
    x = x_new;
  end
end
```

```octave
function x = sor(A, b, x0, omega, tol, maxit)
  n = length(b);
  x = x0;
  for k = 1:maxit
    x_old = x;
    for i = 1:n
      s = b(i) - A(i,:)*x + A(i,i)*x(i);
      x(i) = (1-omega)*x(i) + omega*s/A(i,i);
    end
    if norm(x-x_old) < tol, break; end
  end
end
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000120
md"""
## Problemas
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000121
md"""
### 4.1)

Sea ``(\cdot,\cdot)`` un producto escalar en un espacio vectorial ``X``. Demostrar que la función ``\|x\|=(x,x)^{1/2}`` define una norma en ``X``. *Sugerencia:* la desigualdad triangular se sigue de la desigualdad de Cauchy-Schwarz.
"""

# ╔═╡ d992510f-3a20-4b23-9de9-3a93fd89f7d3
md"""
**Demostración**

1. Por Ax. III de producto escalar, ``(x, x)\geq 0 \Longrightarrow (x, x)^{1/2}\geq 0``.


2. ``\|\alpha x\| = (\alpha x, \alpha x)^{1/2} = \big[\alpha (x, \alpha x)\big]^{1/2} = \big[\alpha^2 (x, x)\big]^{1/2} = |\alpha|(x, x)^{1/2} = |\alpha|\|x\|``


3. Veamos que ``\|x+y\|\leq \|x\| + \|y\|``.
```math
\begin{align*}
	\|x + y\|^2 &= (x+y, x+y)\\
		&= (x, x+y) + (y, x+y)\\
		&= (x, x) + (x, y) + (y, x) + (y,y)\\
		&= (x, x) + 2(x, y) + (y,y)\\
		&\leq (x, x) + 2(x, x)^{1/2}(y, y)^{1/2} + (y,y) &\text{(Por desigualdad C-S)}\\
		&= \big[(x,x)^{1/2} + (y,y)^{1/2}\big]^{2}\\
		&= \big[\|x\| + \|y\|\big]^{2}\\
\end{align*}
```

Luego, como ``\|x\|\geq 0`` y ``\|y\|\geq 0`` se tiene ``\|x+y\|\leq \|x\| + \|y\|``.

------------------
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000122
md"""
### 4.2)
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000140
md"""
#### 4.2.a)

Demostrar que ``\|\cdot\|_1`` y ``\|\cdot\|_\infty`` definen normas en ``\mathbb{R}^n``.
"""

# ╔═╡ c2c61466-811a-4cdd-aa57-a4999a9475be
md"""
-----------------
**Proposición:** ``\|\cdot\|_1 = \sum_{i=1}^{n}|x_i|`` satisface los axiomas de norma.

**Demostración**
1. Como ``|x_i|\geq 0``, se tiene que ``\sum_{i=1}^{n} |x_i| \geq 0``.


2. ``\|\alpha x\| = \sum_{i=1}^{n}|\alpha x_i| = |\alpha| \sum_{i=1}^{n}|x_i|=|\alpha| \|x\|``


3. ``\|x + y\| = \sum_{i=1}^n |x_i + y_i| \leq \sum_{i=1}^n |x| + \sum_{i=1}^n |y| = \|x\| + \|y\|``
"""

# ╔═╡ ccb62783-92d4-4a6a-ba28-95fe20a384cd
md"""
----------------
**Proposición:** ``\|\cdot\|_\infty = \max_{1\leq i\leq n} |x_i|`` satisface los axiomas de norma.

**Demostración**
1. Como cada ``|x_i|\geq 0``, se tiene que ``\|x\|_\infty = \max_{1\leq i\leq n} |x_i|\geq0``.


2. ``\|\alpha x\|_\infty = \max_{1\leq i\leq n} |\alpha x_i| = |\alpha| \max_{1\leq i\leq n} |x_i| = |\alpha|\|x\|_\infty``


3. ``\|x + y\| = \max_{1\leq i\leq n} |x_i + y_i| \leq \max_{1\leq i \leq n}|x_i| + \max_{1\leq i \leq n}|y_i| = \|x\| + \|y\|``
-------------------
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000141
md"""
#### 4.2.b)

Demostrar que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000123
md"""
```math
\|x\|_\infty \leq \|x\|_2 \leq \|x\|_1 \qquad\text{y}\qquad \|x\|_1 \leq \sqrt n\|x\|_2 \leq n\|x\|_\infty,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000124
md"""
para todo ``x\in\mathbb{R}^n``. Demostrar que estas estimaciones son *sharp*, mostrando que existe algún ``x`` para el que se cumple la igualdad en cada una de estas estimaciones.
"""

# ╔═╡ 3640fdd2-ee7b-40f2-967d-d37c00a68115
md"""
**Demostración**

1. Veamos que ``\|x\|_\infty \leq \|x\|_2`` o que ``\|x\|_\infty^2 \leq \|x\|_2^2``

```math
\begin{align*}
	\|x\|_2^2 &= \sum_{i=1}^{n}x_{i}^2 \geq \max_{1\leq i\leq n}|x_i|^2 = \|x\|_\infty^2.
\end{align*}
```

2. Veamos que ``\|x\|_2 \leq \|x\|_1`` o que ``\|x\|_2^2 \leq \|x\|_1^2``
```math
\begin{align*}
	\|x\|_2^2 &= \sum_{i=1}^{n}x_{i}^2 \stackrel{🍎}{\leq} \bigg(\sum_{i=1}^{n}|x_i|\bigg)^2 = \|x\|_1^2
\end{align*}
```

Para demostrar 🍎 supongamos que ``x = (x_i)`` son números positivos y hagamos la demostración por inducción sobre la cantidad de términos de una suma.

- Si ``m=2``, entonces ``(\sum_{i=1}^2 x_i)^2=x_1^2 + 2x_1x_2 + x_2^2 = \sum_{i=1}^2x_i^2 + 2x_1x_2\geq \sum_{i=1}^2x_i^2``.


- Supongamos que se cumple para ``m=n``, es decir, ``(\sum_{i=1}^n x_i)^2 \geq \sum_{i=1}^n x_i^2``. Veamos si se cumple para ``m=n+1``.

```math
\begin{align*}
	\bigg(\sum_{i=1}^{n+1}x_i\bigg)^2 &= \bigg(\sum_{i=1}^{n}x_i + x_{n+1}\bigg)^2\\
	&= \bigg(\sum_{i=1}^{n}x_i\bigg)^2 + 2\sum_{i=1}^{n}x_ix_{n+1} + x_{n+1}^2\\
	&\leq  \sum_{i=1}^{n}x_i^2 + 2\sum_{i=1}^{n}x_ix_{n+1} + x_{n+1}^2 &\text{(Hipótesis inductiva)}\\
	&= \sum_{i=1}^{n+1}x_i^2 + 2\sum_{i=1}^{n}x_ix_{n+1}\\
	&\leq \sum_{i=1}^{n+1}x_i^2 &\text{(Pues $x_i\geq0$)}
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000142
md"""
#### 4.2.c)

Dibujar el conjunto ``\{x\in\mathbb{R}^2 : \|x\|=1\}`` para las normas ``\|\cdot\|_1``, ``\|\cdot\|_2``, ``\|\cdot\|_\infty`` y ``\|\cdot\|_A``. Para la norma ``\|\cdot\|_A`` tomar ``A=\begin{bmatrix}2&0\\0&1\end{bmatrix}``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000125
md"""
### 4.3)

Si ``\|\cdot\|`` denota una norma en ``\mathbb{R}^n``, entonces en el espacio de las matrices cuadradas de orden ``n`` se define la norma inducida que denotamos con el mismo símbolo de la siguiente manera:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000126
md"""
```math
\|A\| := \max_{x\neq0}\frac{\|Ax\|}{\|x\|}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000127
md"""
Demostrar que:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000143
md"""
#### 4.3.a)

``\|\cdot\|`` es efectivamente una norma en ``\mathbb{R}^{n\times n}``.
"""

# ╔═╡ b5ba5312-5d6f-4a12-b5e9-e236bf0bf5bd
md"""
**Demostración**

1. ``\|A\| \geq 0``, pues ``\|Ax\|\geq0`` y ``\|x\|\geq 0`` por definición de norma de vectores.


2.
```math
\begin{align*}
\|\alpha A\| &= \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|\alpha Ax\|}{\|x\|}\\ 
	&= \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{|\alpha|\|Ax\|}{\|x\|}\\
	&= |\alpha|\max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|Ax\|}{\|x\|}\\
	&= |\alpha|\|A\|.
\end{align*}
```


3.
```math
\begin{align*}
\|A + B\| &= \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|(A+B)x\|}{\|x\|}\\
	&= \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|Ax+Bx\|}{\|x\|}\\
	&\leq \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|Ax\|+\|Bx\|}{\|x\|}\\
	&= \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\bigg[\frac{\|Ax\|}{\|x\|}+\frac{\|Bx\|}{\|x\|}\bigg]\\
	&\leq \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|Ax\|}{\|x\|} + \max_{\substack{x\in\mathbb{R}^n\\x\neq 0}}\frac{\|Bx\|}{\|x\|}\\
	&= \|A\| + \|B\|
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000144
md"""
#### 4.3.b)

``\displaystyle\|A\| \stackrel{🍉}{=} \max_{\|x\|=1}\|Ax\|\stackrel{🍌}{=} \max_{\|x\|\leq1}\|Ax\| \stackrel{🍎}{=} \max_{0<\|x\|\leq1}\frac{\|Ax\|}{\|x\|}``.
"""

# ╔═╡ b1690f45-8f68-42d9-973d-96842ac94b67
md"""
**Demostración**

Por definición,

```math
\begin{align*}
	\|A\| &= \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}\frac{\|Ax\|}{\|x\|}\\
		&= \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}\frac{1}{\|x\|}{\|Ax\|}\\
		&= \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}{\Bigg\|\frac{1}{\|x\|}Ax\Bigg\|}\\
		&= \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}{\Bigg\|A\frac{x}{\|x\|}\Bigg\|} &👉👽\\
		&= \max_{\substack{x\in\mathbb{R}^n\\x\neq0}}{\|Au\|} &\text{donde $u=x/\|x\|$}\\
\end{align*}
```

Como ``u`` es cualquier vector unitario con dirección y sentido de ``x``, tenemos demostrado 🍉. (Para demostrar las otras, usar 👽)

"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000128
md"""
### 4.4)

Sea ``\|\cdot\|`` una norma en ``\mathbb{R}^n`` y denotemos también por ``\|\cdot\|`` a la norma inducida en el espacio de matrices de ``n\times n``. Demostrar que:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000145
md"""
#### 4.4.a)

``\|Ax\|\leq\|A\|\|x\|``, para toda matriz ``A`` y para todo vector ``x``.
"""

# ╔═╡ aff826c8-7bab-4794-911e-0d4507ae4034
md"""
**Demostración**

Dado que 
```math
\|A\| = \max_{x\neq0}\frac{\|Ax\|}{\|x\|} \geq\frac{\|Ax\|}{\|x\|}
```
para cada ``x\neq 0``. Luego, ``\|A\|\|x\|\geq \|Ax\|``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000146
md"""
#### 4.4.b)

``\|AB\|\leq\|A\|\|B\|``, para todo par de matrices ``A`` y ``B``.
"""

# ╔═╡ 2666d455-03fa-43df-be85-94e89aabc31c
md"""
**Demostración**

```math
\begin{align*}
	\|AB\| &= \max_{\|u\|=1}\|(AB)u\|\\ 
		&= \max_{\|u\|=1}\|A(Bu)\|\\
		&\leq \max_{\|u\|=1}\|A\|\|(Bu)\|\\
		&\leq \|A\|\max_{\|u\|=1}\|(Bu)\|\\
		&\leq \|A\|\|B\|\\
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000129
md"""
### 4.5)

Consideremos las normas matriciales inducidas por las normas ``\|\cdot\|_1``, ``\|\cdot\|_\infty``, y ``\|\cdot\|_2`` en ``\mathbb{R}^n``. Demostrar que:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000147
md"""
#### 4.5.a)

``\displaystyle\|A\|_1 = \max_{1\leq j\leq n}\sum_{i=1}^n|a_{ij}|``
"""

# ╔═╡ b8550327-a521-46c4-bdc1-dddc83f1da7d
md"""
**Demostración**

```math
\begin{align*}
	\|A\|_1 &= \max_{x \neq 0} \frac{\|Ax\|_1}{\|x\|_1} \\
		&= \max_{x\neq 0} \frac{\sum_{i=1}^{n}|(Ax)_i|}{\sum_{i=1}^{n}|x_i|} &👉🍎
\end{align*}
```

La expresión ``|(Ax)_i|`` de 🍎 se puede acotar por

```math
	|(Ax)_i| = \left|\sum_{j=1}^{n} a_{ij}x_{j}\right|\leq \sum_{j=1}^{n} |a_{ij}x_{j}| = \sum_{j=1}^{n} |a_{ij}||x_{j}|.
```

Reemplazando, tenemos
```math
\begin{align*}
	\|A\|_1 &= \max_{x \neq 0} \frac{\|Ax\|_1}{\|x\|_1} \\
		&\leq \max_{x\neq 0} \frac{\sum_{i=1}^{n}\sum_{j=1}^{n} |a_{ij}||x_{j}|}{\sum_{i=1}^{n}|x_i|}\\
		&= \max_{x\neq 0} \frac{\sum_{j=1}^{n}|x_{j}|\sum_{i=1}^{n} |a_{ij}|}{\sum_{k=1}^{n}|x_k|}\\
		&= \max_{x\neq 0} \sum_{j=1}^{n}\frac{|x_{j}|}{\sum_{k=1}^{n}|x_k|}\sum_{i=1}^{n} |a_{ij}|\\
		&= \max_{x\neq 0} \sum_{j=1}^{n}\omega_j\sum_{i=1}^{n} |a_{ij}|\\
		&\leq \max_{x\neq 0} \sum_{j=1}^{n}\omega_j\left(\max_{1\leq k \leq n}\sum_{i=1}^{n} |a_{ik}|\right)\\
		&\leq \max_{x\neq 0}\left(\max_{1\leq k \leq n}\sum_{i=1}^{n} |a_{ik}|\right)\sum_{j=1}^{n}\omega_j\\
		&\leq\max_{1\leq k \leq n}\sum_{i=1}^{n} |a_{ik}|\\
\end{align*}
```
"""

# ╔═╡ 5a5ec9ff-26b1-4166-8c57-f95372ac6d26
md"""
Veamos que esa cota es alcanzable y que, por lo tanto, ese valor es un máximo.

Supongamos que la columna que produce una suma máxima es la columna ``j``. Sea ``x = e_j`` el vector canónico.

```math
\begin{align*}
	\|A\|_1 &= \max_{x \neq 0} \frac{\|Ae_j\|_1}{\|e_j\|_1}\\
			&\geq \|Ae_j\|_1\\
			&= \|\operatorname{col}_j(A)\|_1\\
			&= \sum_{i=1}^{n} |a_{ij}|
\end{align*}
```
es decir, la suma de los elementos de la columna con suma máxima.

Luego
```math
\|A\|_1 = \max_{x \neq 0} \frac{\|Ax\|_1}{\|x\|_1} = \max_{1\leq j \leq n}\sum\operatorname{col}_j(A)
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000148
md"""
#### 4.5.b)

``\displaystyle\|A\|_\infty = \max_{1\leq i\leq n}\sum_{j=1}^n|a_{ij}|``
"""

# ╔═╡ 8f09c463-94c8-4bcf-a001-87d363543817
md"""
**Demostración**

```math
\begin{align*}
	\|A\|_\infty &= \max_{x\neq 0}\frac{\|Ax\|_\infty}{\|x\|_\infty}\\
		&= \max_{x\neq 0}\frac{\displaystyle\max_{1\leq i \leq n} |(Ax)_i|}{\displaystyle\max_{1\leq i \leq n} |x_i|}\\
\end{align*}
```

Supongamos que ``\mu = \max_{1\leq i\leq n} |x_i|``. La componente ``(Ax)_i``, es

```math
	|(Ax)_i| = \sum_{j=1}^{n}|a_{ij}x_j| = \sum_{j=1}^{n}|a_{ij}||x_j| \leq \mu\sum_{j=1}^{n}|a_{ij}|
```

Luego

```math
\begin{align*}
	\|A\|_\infty &= \max_{x\neq 0}\frac{\displaystyle\max_{1\leq i \leq n} |(Ax)_i|}{\displaystyle\max_{1\leq i \leq n} |x_i|}\\
	&\leq \max_{x\neq 0}\frac{\displaystyle\max_{1\leq i \leq n}  \mu\sum_{j=1}^{n}|a_{ij}|}{\mu}\\
	&= \max_{1\leq i \leq n}\sum_{j=1}^{n}|a_{ij}|\\
\end{align*}
```
"""

# ╔═╡ d43d48e0-e0d5-4444-8bcc-0dbd94f05899
md"""
Veamos que existe un vector que produce el valor de la cota. Sea ``f`` la fila que tiene mayor suma de valores absolutos en ``A``, y ``x`` el vector definido por

```math
x_j := \operatorname{sign}(a_{fj}) = \begin{cases}+1 & \text{si } a_{fj}\geq0\\-1 & \text{si } a_{fj}<0\end{cases}
```

es decir, ``1`` en las componentes donde ``f`` tiene componentes no negativas, y ``-1`` en las componentes negativas. Como ``x_j=\pm1`` para todo ``j``, resulta ``x\neq0`` y ``\|x\|_\infty=1``.

Veamos que ``\|Ax\|_\infty`` alcanza exactamente el valor de la cota, distinguiendo la fila ``f`` del resto.

**Para cualquier fila ``i``** (en particular también sirve para ``i=f``), por la desigualdad triangular y ``|x_j|=1``:

```math
|(Ax)_i| = \left|\sum_{j=1}^{n}a_{ij}x_j\right| \leq \sum_{j=1}^{n}|a_{ij}||x_j| = \sum_{j=1}^{n}|a_{ij}| \leq \max_{1\leq k\leq n}\sum_{j=1}^{n}|a_{kj}|
```

(la última desigualdad porque ``f`` es, por definición, la fila de suma máxima). Ninguna fila puede superar ese valor.

Combinando ambos casos: la fila ``f`` alcanza exactamente ``\max_k\sum_j|a_{kj}|``, y ninguna otra fila la supera, así que

```math
\|Ax\|_\infty = \max_{1\leq i\leq n}|(Ax)_i| = \max_{1\leq k\leq n}\sum_{j=1}^{n}|a_{kj}|
```

Como ``\|x\|_\infty=1``:

```math
\|A\|_\infty = \max_{x\neq 0}\frac{\|Ax\|_\infty}{\|x\|_\infty} \geq \frac{\|Ax\|_\infty}{\|x\|_\infty} = \max_{1\leq k\leq n}\sum_{j=1}^{n}|a_{kj}|
```

Junto con la cota superior ya probada, concluimos

```math
\|A\|_\infty = \max_{1\leq i \leq n}\sum_{j=1}^{n}|a_{ij}|
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000149
md"""
#### 4.5.c)

``\|A\|_2 = \max\{\sqrt\lambda : \lambda \text{ es autovalor de } A^TA\}``
"""

# ╔═╡ 799e8331-9b36-4739-85ad-296ccae4a79f
md"""
**Demostración**
"""

# ╔═╡ 0fd65f61-4e9c-4234-834f-f0ccb27e5efe
md"""
```math
\begin{align*}
	\|A\|_2^2 &= \max_{x\neq 0}\frac{\|Ax\|^2_2}{\|x\|_2^2}\\
		&= \max_{x\neq 0} \frac{(Ax)^T(Ax)}{\|x\|_2^2}\\
		&= \max_{x\neq 0} \frac{x^TA^TAx}{\|x\|_2^2}\\
\end{align*}
```
"""

# ╔═╡ 25683840-7722-457d-aeac-18c4e9bdeb94
md"""
Como ``A^TA`` es sdp, es diagonalizable, existe ``Q`` ortogonal tal que 

```math
A^TA = Q^T\Lambda Q
```
donde ``\Lambda`` es una matriz diagonal con los autovalores ``\lambda_u,\ldots,\lambda_i,\ldots,\lambda_l`` de ``A^TA`` en ella. Por lo que 
```math
x^TA^TAx = x^TQ^T\Lambda Qx = (Qx)^T\Lambda Qx
```

Sea ``y=Qx``. Sabemos que ``\|y\|_2^2 = \|Qx\|_2^2 = \|x\|_2^2``.
"""

# ╔═╡ 3d45eabe-2e0f-4bea-b739-0f2c1da5aacd
md"""
```math
\begin{align*}
	\|A\|_2^2 &= \max_{x\neq 0} \frac{x^TA^TAx}{\|x\|_2^2}\\
		&= \max_{x\neq 0} \frac{(Qx)^T\Lambda Qx}{\|x\|_2^2}\\
		&\stackrel{?}{=} \max_{x\neq 0} \frac{\|(Qx)^T\Lambda Qx\|_\infty}{\|x\|_2^2}\\
		&\leq \max_{x\neq 0} \frac{\|(Qx)^T\Lambda Qx\|_\infty}{\|x\|_\infty^2} &&\text{(pues $\|x\|_2 \geq \|x\|_\infty$)}\\
		&\leq \max_{x\neq 0} \frac{\|(Qx)^T\|_\infty\|\Lambda\|_\infty \|Qx\|_\infty}{\|x\|_\infty^2}\\
		&= \max_{x\neq 0} \frac{\|x\|_\infty\|\Lambda\|_\infty \|x\|_\infty}{\|x\|_\infty^2}\\
		&= \max_{x\neq 0} \frac{\|x\|_\infty^2\|\Lambda\|_\infty}{\|x\|_\infty^2}\\
		&= \|\Lambda\|_\infty\\
		&= \lambda_u &&\text{(cada $\lambda_i>0$)}\\
		&= \max\{\lambda:\text{$\lambda$ autovalor de $A^TA$}\}
\end{align*}
```
por lo que ``\|A\|_2=\max\{\sqrt\lambda:\text{$\lambda$ autovalor de $A^TA$}\}``.
"""

# ╔═╡ 07a5ad49-b79c-4778-9d28-6c4572fd34f6
md"""
Veamos que existe ``x`` que permite garantizar la igualdad. Sea ``x = v_{\lambda_u}`` el autovector correspondiente al máximo autovalor ``\lambda_u``.

```math
\begin{align*}
	\|A\|_2^2 &= \max_{x\neq 0} \frac{x^TA^TAx}{\|x\|_2^2}\\
			&\geq \frac{v_{\lambda_u}^TA^TAv_{\lambda_u}}{\|v_{\lambda_u}\|_2^2}\\
			&= \frac{v_{\lambda_u}^T\lambda_{u} v_{\lambda_u}}{\|v_{\lambda_u}\|_2^2}\\
			&= \frac{\lambda_{u} \|v_{\lambda_u}\|_2^2}{\|v_{\lambda_u}\|_2^2}\\
			&= \lambda_u\\
		&= \max\{\lambda:\text{$\lambda$ autovalor de $A^TA$}\}
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000014a
md"""
#### 4.5.d)

Si ``A`` es simétrica, ``\|A\|_2=\max\{|\lambda| : \lambda \text{ es autovalor de } A\}``
"""

# ╔═╡ ab9318d4-58fc-405d-8ccf-e6c5f0b31427
md"""
Antes de comenzar con la demostración, necesitaremos un par de resultados que aun no fueron probados.
"""

# ╔═╡ 14b0bc67-1b97-4080-881e-c7f3b624f101
md"""
**Lema:** Si ``Q`` es ortogonal, entonces ``\|Q\|_2 = 1`` 👉📐

**Demostración:** Por el ejercicio anterior y dado que ``Q^TQ=I``, se tiene directamente que ``\|Q\|_2 = 1``
"""

# ╔═╡ 0de3ff9f-0cff-4849-b709-41cc1bc3d083
md"""
**Lema:** Si ``D`` es una matriz diagonal con elementos ``d_i``, entonces ``\|D\|_2 = \max |d_i|``. 👉🧰

**Demostración:** 

Sea ``\Delta = D^TD``. Sabemos que ``\Delta`` es diagonal, y que sus componentes ``\delta_i`` satisfacen ``\delta_i=d_i^2``. Luego, por ser diagonal, sabemos que sus autovalores ``\lambda_i`` son los elementos de la diagonal, es decir, ``\lambda_i=d_i^2``

Por el ejercicio anterior, tenemos entonces que ``\|D\|_2 = \max\left\{\sqrt{d_i^2}\right\} = \max |d_i|``.
"""

# ╔═╡ feeeb893-bc14-4929-858b-0a2dbd093fa0
md"""
Ahora sí, continuamos con la demostración del ejercicio.
"""

# ╔═╡ 4e4030c2-2b0d-4116-b2f7-f974034dd035
md"""
**Demostración**
"""

# ╔═╡ 0ecd695b-c019-4d90-9248-c7dc107e2341
md"""
```math
\begin{align*}
\|A\|_2 &= \max_{x\neq0}\frac{\|Ax\|_2}{\|x\|_2}\\
	&= \max_{x\neq0}\frac{\|Q^T\Lambda Qx\|_2}{\|x\|_2}\\
	&\leq \max_{x\neq0}\frac{\|Q^T\Lambda\|_2 \|Qx\|_2}{\|x\|_2}\\
	&= \max_{x\neq0}\frac{\|Q^T\Lambda\|_2 \|x\|_2}{\|x\|_2}\\
	&= \|Q^T\Lambda\|_2\\
	&\leq \|Q^T\|_2\|\Lambda\|_2\\
	&= \|\Lambda\|_2 && 📐\\
	&= \max\{|\lambda_i| : \lambda_i \text{ autovalor de $A$}\} &&🧰
\end{align*}
```
"""

# ╔═╡ 7d42abf0-b539-465a-bda9-94a2477df4ed
md"""
Sea ahora ``x=v_\omega``, donde ``v`` es el autovector correspondiente al autovalor ``\omega`` que presenta mayor valor absoluto.
"""

# ╔═╡ b95f1143-162c-4e4b-be4f-35a5170dccc9
md"""
```math
\begin{align*}
\|A\|_2 &= \max_{x\neq0}\frac{\|Ax\|_2}{\|x\|_2}\\
	&\geq \frac{\|Av_\omega\|_2}{\|v_\omega\|_2}\\
	&= \frac{\|\omega v_\omega\|_2}{\|v_\omega\|_2}\\
	&= \frac{|\omega|\|v_\omega\|_2}{\|v_\omega\|_2}\\
	&= |\omega|\\
	&= \max\{|\lambda_i| : \lambda_i \text{ autovalor de $A$}\}
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000012a
md"""
### 4.6)

Sea ``A\in\mathbb{R}^{n\times n}`` una matriz triangular. Probar que sus autovalores son los valores de la diagonal.
"""

# ╔═╡ 6703414f-d602-4cbb-bb49-ab47ff49f18f
md"""
**Demostración** Sabemos que el determinante de una matriz triangular es el producto de los elementos de la diagonal. Sea ``A`` una matriz triangular. Entonces, sus autovalores vienen dados por la solución a

```math
\begin{align*}
	\det(A - \lambda I) = \det \tilde A= 0. 
\end{align*}
```

Ahora bien, ``\tilde A`` es triangular, pues es como ``A``, pero con las componentes de la diagonal modificadas. Más aun, ``\tilde a_{ii} = a_{ii} - \lambda``. 

Por lo tanto, 

```math
\det(A-\lambda I) = \det\tilde A = \prod_{i=1}^{n} \tilde a_{ii} = \prod_{i=1}^{n} (a_{ii} - \lambda) \stackrel{🍎}{=} 0.
```

Ahora bien, 🍎 se cumple solo si ``a_{ii}-\lambda = 0`` para cada ``i`` o, lo que es lo mismo, ``a_{ii} = \lambda``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000012b
md"""
### 4.7)

Sea ``\rho(A):=\max\{|\lambda| : \lambda \text{ autovalor de } A\}`` el radio espectral de ``A``. Demostrar que:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000014b
md"""
#### 4.7.a)

``\rho(A)\leq\|A\|`` para cualquier norma inducida y para toda matriz ``A``.
"""

# ╔═╡ b58bb658-effc-49e9-b7e7-2bdc3c61c9fa
md"""
**Demostración**

Por definición, savemos que 

```math
\begin{align*}
	\lambda v &= A v\\
	\|\lambda v\| &= \|Av\|\\
	|\lambda|\|v\| &\leq \|A\| \|v\| \\
	|\lambda| &\leq \|A\| &&\text{(pues $\|v\|\neq 0$)}\\
	\rho(A) &\leq \|A\|   &&\text{(pues se cumple para cualquier autovalor)}
\end{align*}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000014c
md"""
#### 4.7.b)

``\rho(A)`` puede ser cero aunque ``A`` no lo sea.
"""

# ╔═╡ c5783e27-39d5-4415-b1eb-eaa96c75ba38
md"""
**Demostración**

Sea ``A`` una matriz triangular con su diagonal nula y al menos un elemento no nulo del resto de las componentes. Entonces

1. ``A`` no es nula.
2. Por el ejercicio anterior, ``\sigma(A)`` es un conjunto de ceros. Así, ``\rho(A)`` es cero, cuando ``A`` no lo es.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000012c
md"""
### 4.8)

Sea ``\|\cdot\|`` una norma en ``\mathbb{R}^n`` y sea ``G\in\mathbb{R}^{n\times n}`` una matriz invertible. Demostrar que
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000012d
md"""
```math
\|x\|_G := \|Gx\|, \qquad \forall x\in\mathbb{R}^n,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000012e
md"""
es una norma en ``\mathbb{R}^n`` y que la norma matricial inducida está dada por
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000012f
md"""
```math
\|A\|_G := \|GAG^{-1}\|, \qquad \forall A\in\mathbb{R}^{n\times n}.
```
"""

# ╔═╡ 11db1c2a-7452-49af-809e-47ca58ad3839
md"""
**Demostración**


"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000130
md"""
### 4.9)

Sea ``M\in\mathbb{R}^{n\times n}``. Demostrar que si ``\rho(M)\geq1`` entonces existen ``x_0,c\in\mathbb{R}^n`` tales que la iteración
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000131
md"""
```math
x_{k+1} = Mx_k+c, \qquad k\geq0,
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000132
md"""
no converge.
"""

# ╔═╡ bd108856-bc61-41f6-b2b1-350c8d7dd551
md"""
**Demostración**


"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000133
md"""
### 4.10)

Demostrar que si una matriz ``A\in\mathbb{R}^{n\times n}`` es *estrictamente diagonalmente dominante* (e.d.d), entonces resulta ser no singular.

*(Sugerencia: Probar que ``Ax=0\implies x=0``)*
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000134
md"""
### 4.11)

Estudiar analíticamente la convergencia de los métodos de Jacobi y Gauss-Seidel, para cualesquiera valores iniciales, para el sistema
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000135
md"""
```math
\begin{aligned}
x &&&+&& z &= 2\\
-x &+& y &&&&= 0\\
x &+& 2y &-& 3z &&= 0
\end{aligned}
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000136
md"""
**Observación:** La condición de edd es suficiente pero no necesaria para la convergencia de estos métodos.

*Ayuda:* ¿Conoces alguna condición necesaria y suficiente para la convergencia de métodos de este tipo?
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000137
md"""
### 4.12)

Demostrar que si ``A`` es sdp, entonces todas las entradas de su diagonal son positivas.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000138
md"""
### 4.13)

Demostrar que si ``A`` es sdp, entonces el método SOR converge globalmente (i.e., para cualquier aproximación inicial ``x_0``) para cualquier elección de ``\omega\in(0,2)``.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-000000000139
md"""
### 4.14)

Crear las funciones `x = triinf(L,b)` y `x = trisup(U,b)` que resuelvan sistemas con matrices cuadradas triangulares inferiores/superiores con elementos no nulos en la diagonal. Para cada función, utilizar a lo sumo un lazo `for`.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000013a
md"""
### 4.15)

Programar los métodos de Jacobi, Gauss-Seidel, y SOR y utilizarlos para resolver los siguientes sistemas.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000013b
md"""
```math
(i)\quad \begin{bmatrix}2&1&-1\\-2&-10&0\\-1&-1&4\end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix} = \begin{bmatrix}1\\-12\\2\end{bmatrix},
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000013c
md"""
```math
(ii)\quad \begin{bmatrix}10&1&2&3&4\\1&9&-1&2&-3\\2&-1&7&3&-5\\3&2&3&12&-1\\4&-3&-5&-1&15\end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\\x_4\\x_5\end{bmatrix} = \begin{bmatrix}12\\-27\\14\\-17\\12\end{bmatrix}.
```
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000013d
md"""
Para cada uno de los sistemas:
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000014d
md"""
#### 4.15.a)

Realizar un gráfico del decrecimiento del error (en escala logarítmica) en términos del número de iteraciones para cada uno de los métodos utilizados. Extraer conclusiones.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000014e
md"""
#### 4.15.b)

Eligiendo alguna norma conveniente, realizar un gráfico de ``\frac{\|x_{k+1}-x\|}{\|x_k-x\|}`` en término del número de iteraciones. ¿Cuál es asintóticamente el factor de reducción del error para cada método? Comparar con el radio espectral de cada matriz de iteración.
"""

# ╔═╡ d4e5f6a7-8888-4a1b-8c2d-00000000013e
md"""
### 4.16)

Encontrar una matriz ``A`` para la cual sea imposible encontrar ``\alpha`` de manera que la iteración de Richardson ``x_{k+1}=x_k-\alpha(Ax_k-b)`` converja para todo ``x_0``.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Luxor = "~4.5.0"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "51783892003ed8b9dcc6bed99e443755cff7f722"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1fa950ebc3e37eccd51c6a8fe1f92f7d86263522"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.7+0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.1+2"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f4d39eee89f1e58c26bf447f1d4156c0125d6838"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "7a58e45171b63ed4782f2d36fdee8713a469e6e0"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.2+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "090526e65de8f69648ac156daae153de8b56df62"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.88.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "c3f99f8e7c98031b8845ece9db86dca713ab9bf8"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.14003.0+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b7970cef8ae1c990ba0c09cd8bdc1145e006632f"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "22.1.7+0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "Cairo_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pango_jll", "XML2_jll", "gdk_pixbuf_jll"]
git-tree-sha1 = "e6ab5dda9916d7041356371c53cdc00b39841c31"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.7+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "PolygonAlgorithms", "PrecompileTools", "Random", "Rsvg"]
git-tree-sha1 = "fe8060b3d693f682e14f1019b058c64effb62b43"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "4.5.0"

    [deps.Luxor.extensions]
    LuxorExtLatex = ["LaTeXStrings", "MathTeXEngine"]
    LuxorExtTypstry = ["Typstry"]

    [deps.Luxor.weakdeps]
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    MathTeXEngine = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolygonAlgorithms]]
git-tree-sha1 = "c1092ada65e6d59d6361d5086ddb0a5ea63ae204"
uuid = "32a0d02f-32d9-4438-b5ed-3a2932b48f96"
version = "0.4.0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "e53dad0507631c0b8d5d946d93458cbabd0f05d7"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.1.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "895f21b699121d1a57ecac57e65a852caf569254"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.13+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ef17c47d22224aaecc76e597ab21a072e025cf7b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.14.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cb007192783c56d8249db4cf0e3495001edfe414"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.5+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "28e57478e8a160d346a19c28b3fffb9273bcc9c2"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.134+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╠═d4e5f6a7-8888-4a1b-8c2d-000000000001
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000002
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000003
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000004
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000005
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000006
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000007
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000008
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000009
# ╟─2f9d20ff-9e4b-439a-9c00-fc8565d4bf61
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000000a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000000b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000000c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000000d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000000e
# ╟─f0437696-4357-42f8-89ea-6e3a254a0c2a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000000f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000010
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000011
# ╟─cf19f3da-e87a-4624-a896-4b6276409318
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000012
# ╟─75fb8468-25e0-4506-bf65-1c2d93d58f84
# ╟─28d73e21-98c3-40c1-a2c2-e4f78a605f5a
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000013
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000014
# ╟─dfcb451e-7d09-4a8a-8f60-a44569572824
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000015
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000016
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000017
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000018
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000019
# ╟─5466d779-1be8-40ac-9126-60fc3161784d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000001a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000001b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000001c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000001d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000001e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000001f
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000001f1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000001f2
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000020
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000021
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000022
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000023
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000024
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000025
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000251
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000252
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000026
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000027
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000028
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000029
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000002a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000002b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000002c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000002d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000002e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000002f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000030
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000031
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000032
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000033
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000034
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000035
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000036
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000037
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000038
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000039
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000003a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000003b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000003c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000003d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000003e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000003f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000040
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000041
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000042
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000043
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000044
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000045
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000046
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000047
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000048
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000049
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000004a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000004b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000004c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000004d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000004e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000004f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000050
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000051
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000052
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000053
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000054
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000055
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000056
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000057
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000058
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000059
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000005a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000005b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000005c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000005d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000005e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000005f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000060
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000061
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000062
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000063
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000064
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000065
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000066
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000067
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000068
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000069
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000006a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000006b
# ╟─f6ca739b-2bd9-46f1-88c0-4f277da6f81b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000006c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000006d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000006e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000006f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000070
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000071
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000072
# ╟─7e76323b-8d03-4454-ada1-7e175d930a43
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000073
# ╟─afa824e3-52c4-41b0-a94a-53a793ccaeaf
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000074
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000075
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000076
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000077
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000078
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000079
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000007a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000007b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000007c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000007d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000007e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000007f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000080
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000081
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000082
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000083
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000084
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000085
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000086
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000087
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000088
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000089
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000008a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000008b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000008c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000008d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000008e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000008f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000090
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000091
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000092
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000093
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000094
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000095
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000096
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000097
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000098
# ╟─e9ed9b75-1f05-4c90-bb5c-94c5e9325820
# ╟─18d8ba73-e4aa-4a16-a83e-60c4243466a2
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000099
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000009a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000009b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000009c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000009d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000009e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000009f
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a0
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a2
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a3
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a4
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a5
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a6
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a7
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a8
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000a9
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000aa
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ab
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ac
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ad
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ae
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000af
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b0
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b2
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b3
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b4
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b5
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b6
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b7
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b8
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000b9
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ba
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000bb
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000bc
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000bd
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000be
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000bf
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c0
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c2
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c3
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c4
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c5
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c6
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c7
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c8
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000c9
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ca
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000cb
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000cc
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000cd
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ce
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000cf
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d0
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d2
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d3
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d4
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d5
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d6
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d7
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d8
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000d9
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000da
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000db
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000dc
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000dd
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000de
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000df
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e0
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e2
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e3
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e4
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e5
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e6
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e7
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e8
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000e9
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ea
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000eb
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ec
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ed
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ee
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ef
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f0
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f1
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f2
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f3
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f4
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f5
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f6
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f7
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f8
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000f9
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000fa
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000fb
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000fc
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000fd
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000fe
# ╟─d4e5f6a7-8888-4a1b-8c2d-0000000000ff
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000100
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000101
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000102
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000103
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000104
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000105
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000106
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000107
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000108
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000109
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000010a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000010b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000010c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000010d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000010e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000010f
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000110
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000111
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000112
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000113
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000114
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000115
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000116
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000117
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000118
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000119
# ╟─d4e5f6a7-9999-4a1b-8c2d-000000000001
# ╟─d4e5f6a7-9999-4a1b-8c2d-000000000002
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000120
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000121
# ╟─d992510f-3a20-4b23-9de9-3a93fd89f7d3
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000122
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000140
# ╟─c2c61466-811a-4cdd-aa57-a4999a9475be
# ╟─ccb62783-92d4-4a6a-ba28-95fe20a384cd
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000141
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000123
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000124
# ╟─3640fdd2-ee7b-40f2-967d-d37c00a68115
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000142
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000125
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000126
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000127
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000143
# ╟─b5ba5312-5d6f-4a12-b5e9-e236bf0bf5bd
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000144
# ╟─b1690f45-8f68-42d9-973d-96842ac94b67
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000128
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000145
# ╟─aff826c8-7bab-4794-911e-0d4507ae4034
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000146
# ╟─2666d455-03fa-43df-be85-94e89aabc31c
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000129
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000147
# ╟─b8550327-a521-46c4-bdc1-dddc83f1da7d
# ╟─5a5ec9ff-26b1-4166-8c57-f95372ac6d26
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000148
# ╟─8f09c463-94c8-4bcf-a001-87d363543817
# ╟─d43d48e0-e0d5-4444-8bcc-0dbd94f05899
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000149
# ╟─799e8331-9b36-4739-85ad-296ccae4a79f
# ╟─0fd65f61-4e9c-4234-834f-f0ccb27e5efe
# ╟─25683840-7722-457d-aeac-18c4e9bdeb94
# ╟─3d45eabe-2e0f-4bea-b739-0f2c1da5aacd
# ╟─07a5ad49-b79c-4778-9d28-6c4572fd34f6
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000014a
# ╟─ab9318d4-58fc-405d-8ccf-e6c5f0b31427
# ╟─14b0bc67-1b97-4080-881e-c7f3b624f101
# ╟─0de3ff9f-0cff-4849-b709-41cc1bc3d083
# ╟─feeeb893-bc14-4929-858b-0a2dbd093fa0
# ╟─4e4030c2-2b0d-4116-b2f7-f974034dd035
# ╟─0ecd695b-c019-4d90-9248-c7dc107e2341
# ╟─7d42abf0-b539-465a-bda9-94a2477df4ed
# ╟─b95f1143-162c-4e4b-be4f-35a5170dccc9
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000012a
# ╟─6703414f-d602-4cbb-bb49-ab47ff49f18f
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000012b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000014b
# ╟─b58bb658-effc-49e9-b7e7-2bdc3c61c9fa
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000014c
# ╟─c5783e27-39d5-4415-b1eb-eaa96c75ba38
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000012c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000012d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000012e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000012f
# ╟─11db1c2a-7452-49af-809e-47ca58ad3839
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000130
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000131
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000132
# ╟─bd108856-bc61-41f6-b2b1-350c8d7dd551
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000133
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000134
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000135
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000136
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000137
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000138
# ╟─d4e5f6a7-8888-4a1b-8c2d-000000000139
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000013a
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000013b
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000013c
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000013d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000014d
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000014e
# ╟─d4e5f6a7-8888-4a1b-8c2d-00000000013e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
