### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ b2c3d4e5-3333-4a1b-8c2d-000000000003
md"""
# Capítulo 2: Matrices estrictamente diagonalmente dominantes (2.11 y 2.12)
"""


# ╔═╡ 0ac08613-e31a-45b7-8e5d-9b3deb334da2
md"""
### 2.11
**Proposición:** Si ``A`` es *estrictamente diagonalmente dominante* (e.d.d.), entonces es no singular.


**Demostración**

Veamos, por induccion, que reducir por Gauss produce siempre pivotes no nulos. Para ello, veremos que si un pivote es no nulo y la submatriz restante es e.d.d., entonces la submatriz que se produce al eliminar debajo de ese pivote también es e.d.d.

Sea
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000021
md"""
```math
	A = \begin{bmatrix}
		a_{11} & a_{12} & a_{13} & \cdots & a_{1n}\\
		a_{21} & a_{22} & a_{23} & \cdots & a_{2n}\\
		\vdots & \vdots & \ddots & \ddots & \vdots\\
		a_{n1} & a_{n2} & a_{n3} & \cdots & a_{nn}
	\end{bmatrix}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000022
md"""
**Caso inicial:** Como ``A`` es e.e.d., ``|a_{11}| > \sum_{j\neq 1}|a_{1j}| \geq 0``, en particular ``a_{11}\neq 0``. Al eliminar la primera columna:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000023
md"""
```math
\def\celda#1#2{a_{#1#2} - a_{#11}\dfrac{a_{1#2}}{a_{11}}}

\begin{align*}
	P_1 &= \begin{bmatrix}
		a_{11} & a_{12} & a_{13} & \cdots & a_{1n}\\
		0 & \celda{2}{2} & \celda{2}{3} & \cdots & \celda{2}{n}\\
		\vdots & \vdots & \ddots & \ddots & \vdots\\
		0 & \celda{n}{2} & \celda{n}{3} & \cdots & \celda{n}{n}\\
	\end{bmatrix}\\
		&= \begin{bmatrix}
			a_{11} & a_{12} & a_{13} & \cdots & a_{1n}\\
			0 & b_{11} & b_{12} & \cdots & b_{1,n-1}\\
			\vdots & \vdots & \ddots & \ddots & \vdots\\
			0 & b_{n-1;1} & b_{n-1,2} & \cdots & b_{n-1;n-1}
		\end{bmatrix}
\end{align*}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000024
md"""
Sea ``m=n-1``. Veamos que
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000025
md"""
```math
B = \begin{bmatrix}
	b_{11} & b_{12} & \cdots & b_{1m}\\
	b_{21} & b_{22} & \cdots & b_{2m}\\
	\vdots & \vdots & \ddots & \vdots\\
	b_{m1} & b_{m2} & \cdots & b_{mm}
\end{bmatrix}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000026
md"""
es tambien e.d.d. 

Para cada índice ``2\leq i\leq n``, sea
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005c
md"""
```math
	D=\{2,\dots,n\}\setminus\{i\}\quad (👉👽)\qquad\text{y}\qquad \lambda_i = |a_{i1}|/|a_{11}|\quad (👉🐔)
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005d
md"""
Como ``A`` es e.d.d. en la fila 1, ``|a_{i1}| < |a_{11}|``, así que ``\lambda_i \in [0,1)``.

Veamos, ahora sí, que ``|b_{ii}| > \sum_{j\neq i}|b_{ij}| \geq 0`` acotando por separado la diagonal y la suma de ``b``. 

Por la desigualdad triangular
"""

# ╔═╡ 7b9157c7-b2b6-4f03-bcfc-745af8843d25
md"""
```math
|a| - |b| \stackrel{👉🎧}{\leq} |a-b| \stackrel{👉🎤}{\leq} |a| + |b|
```
"""

# ╔═╡ 470bde64-8a4b-459e-86fc-9cc4292dfbb4
md"""
que, aplicada a cada uno de los términos de interés, nos da
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005e
md"""
```math
\begin{align*}
	👉🐱&& |b_{ii}| &= \bigg|a_{ii} - a_{i1}\dfrac{a_{1i}}{a_{11}}\bigg| \stackrel{\substack{🎧\\🐔}}{\geq} |a_{ii}| - \lambda_i|a_{1i}|\\
	👉🐶&& \sum_{j\in D}|b_{ij}| &= \sum_{j\in D}\bigg|a_{ij} - a_{i1}\dfrac{a_{1j}}{a_{11}}\bigg| \stackrel{\substack{🎤\\🐔}}{\leq} \sum_{j\in D}|a_{ij}| + \lambda_i\sum_{j\in D}|a_{1j}|
\end{align*}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005f
md"""
Para ver que la componente de la diagonal es mayor que la suma, veamos que la diferencia entre la primera y la segunda es mayor que cero. Restando ambas cotas:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000060
md"""
```math
\begin{align}
	👉🍌\qquad |b_{ii}| - \sum_{j\in D}|b_{ij}| & \stackrel{\substack{🐱\\🐶}}{\geq}
		|a_{ii}| - \lambda_i|a_{1i}| - \bigg[\sum_{j\in D}|a_{ij}| + \lambda_i\sum_{j\in D}|a_{1j}|\bigg]\\
	&\geq \bigg[|a_{ii}| - \sum_{j\in D}|a_{ij}|\bigg] - \lambda_i\bigg[|a_{1i}| + \sum_{j\in D}|a_{1j}|\bigg]
\end{align}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000061
md"""
El término ``|a_{1i}| + \sum_{j\in D}|a_{1j}|`` es la suma completa de la fila 1 sin el pivote (porque ``i\notin D``), es decir,
"""

# ╔═╡ eff4a025-7755-4bc5-8fc7-5aba85fa41f4
md"""
```math
👉🍉\qquad|a_{1i}| + \sum_{j\in D}|a_{1j}| = \sum_{j\neq 1}|a_{1j}|
```
"""

# ╔═╡ cb2e75a2-1fb2-4b3a-8cd0-e97244b055fd
md"""
Usamos las hopótesis sobre la fila 1 de ``A`` para acotar el segundo corchete, y la de de la fila ``i`` de ``A`` (sumando todos los ``j\neq i``, incluyendo ``j=1``) para acotar el primero:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000062
md"""
```math
\begin{align*}
	👉🍎 &&\lambda_i\sum_{j\neq 1}|a_{1j}| &< \lambda_i|a_{11}| =\frac{|a_{i1}|}{|a_{11}|} |a_{11}|  = |a_{i1}|\qquad\text{y}\\
	👉🍐 &&|a_{ii}| - \sum_{j\in D}|a_{ij}| &> |a_{i1}| \qquad\text{pues}\qquad |a_{ii}| > \sum_{j\neq i}|a_{ij}| = |a_{i1}| + \sum_{j\in D}|a_{ij}|
\end{align*}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000063
md"""
Encadenando:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000064
md"""
```math
	👉🍇\qquad|a_{ii}| - \sum_{j\in D}|a_{ij}| \stackrel{🍐}{>} |a_{i1}| \stackrel{🍎}{>} \lambda_i\sum_{j\neq 1}|a_{1j}|
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000065
md"""
y por lo tanto
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000066
md"""
```math
\begin{align}
	|b_{ii}| - \sum_{j\in D}|b_{ij}| & \stackrel{🍌}{\geq}
		\bigg[|a_{ii}| - \sum_{j\in D}|a_{ij}|\bigg] - \lambda_i\bigg[|a_{1i}| + \sum_{j\in D}|a_{1j}|\bigg]\\
	&\stackrel{\genfrac{}{}{0pt}{1}{🍇}{🍉}}{>} \lambda_i\sum_{j\neq 1}|a_{1j}| - \lambda_i \sum_{j\neq 1}|a_{1j}|\\
	& = 0
\end{align}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000067
md"""
por lo que
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000068
md"""
```math
	|b_{ii}| > \sum_{j\in D}|b_{ij}|
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000069
md"""
En particular, como ``\sum_{j\in D}|b_{ij}| \geq 0``, también se obtiene ``|b_{ii}|>0``. Así, ``b`` es estrictamente diagonalmente dominante.

**Paso inductivo:** Habiendo mostrado que la submatriz ``B`` (de tamaño ``m=n-1``) hereda e.e.d., el mismo argumento aplicado a ``B`` muestra que la siguiente submatriz de tamaño ``m-1`` (por ejemplo ``C``) también es e.e.d., y así sucesivamente hasta agotar las ``n-1`` columnas. Por inducción (rigor????), todos los pivotes de la eliminación gaussiana sin pivoteo sobre ``A`` son no nulos. Entonces la matriz ``U`` resultante tiene diagonal no nula, por lo que ``\det(U)\neq 0`` y es no singular.
"""

# ╔═╡ c6f5fd8d-505b-4b55-bd73-e5c028b771e1
md"""
### 2.12

Demostración directa del ejercicio anterior.
"""


# ╔═╡ Cell order:
# ╟─b2c3d4e5-3333-4a1b-8c2d-000000000003
# ╟─0ac08613-e31a-45b7-8e5d-9b3deb334da2
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000021
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000022
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000023
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000024
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000025
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000026
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005d
# ╟─7b9157c7-b2b6-4f03-bcfc-745af8843d25
# ╟─470bde64-8a4b-459e-86fc-9cc4292dfbb4
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005e
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000060
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000061
# ╟─eff4a025-7755-4bc5-8fc7-5aba85fa41f4
# ╟─cb2e75a2-1fb2-4b3a-8cd0-e97244b055fd
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000062
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000063
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000064
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000065
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000066
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000067
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000068
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000069
# ╟─c6f5fd8d-505b-4b55-bd73-e5c028b771e1
