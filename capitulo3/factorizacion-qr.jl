### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000001
using LinearAlgebra, PlutoUI, PlutoTeachingTools, BenchmarkTools, GLM

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000002
TableOfContents()

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000003
md"""
# Capítulo 3: Factorización QR
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000004
md"""
## Mínimos cuadrados
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000005
md"""
Consideremos la situación en que se quiere *ajustar* (*fitear*) una lista dada de puntos ``(X_i,Y_i)``, ``i=1,2,\dots,n``, por una función cuadrática. Lo más usual es hacer lo que se conoce por *mínimos cuadrados*, que consiste en hallar los parámetros ``\alpha``, ``\beta``, ``\gamma`` que minimizan el funcional
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000006
md"""
```math
J(\alpha,\beta,\gamma) = \sum_{i=1}^{n}[(\alpha X_i^2+\beta X_i+\gamma)-Y_i]^2.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000007
md"""
Si definimos
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000008
md"""
```math
A = \begin{bmatrix}X_1^2&X_1&1\\X_2^2&X_2&1\\\vdots&\vdots&\vdots\\X_n^2&X_n&1\end{bmatrix}, \qquad b = \begin{bmatrix}Y_1\\Y_2\\\vdots\\Y_n\end{bmatrix}, \qquad x = \begin{bmatrix}\alpha\\\beta\\\gamma\end{bmatrix},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000009
md"""
entonces ``J(\alpha,\beta,\gamma) = \|Ax-b\|_2^2``, y el problema consiste en minimizar ``\|Ax-b\|_2``, el cuadrado de la norma euclídea del vector ``Ax-b``.

Consideremos entonces el problema de minimizar ``\|Ax-b\|_2`` para ``A\in\mathbb{R}^{n\times m}``, ``b\in\mathbb{R}^n``, con ``n\geq m``.

Si recordamos que ``Ax`` es una combinación lineal de las columnas de ``A``, el problema equivale a hallar el elemento del espacio columna de ``A`` que está más cerca de ``b``. Dado que usamos la norma euclídea, esto equivale a hallar ``x\in\mathbb{R}^m`` tal que ``Ax-b`` es ortogonal a las columnas de ``A``. Es decir
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000000a
md"""
```math
\begin{cases}\text{col}_1(A)^T(Ax-b)=0\\\text{col}_2(A)^T(Ax-b)=0\\\quad\vdots\\\text{col}_m(A)^T(Ax-b)=0\end{cases} \iff \begin{cases}\text{fila}_1(A^T)(Ax-b)=0\\\text{fila}_2(A^T)(Ax-b)=0\\\quad\vdots\\\text{fila}_m(A^T)(Ax-b)=0\end{cases} \iff A^T(Ax-b)=0.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000000b
md"""
Llegamos así a que para minimizar ``\|Ax-b\|_2`` debemos resolver las *ecuaciones normales*
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000000c
md"""
```math
A^TAx = A^Tb.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000000d
md"""
Si las columnas de ``A`` son linealmente independientes, entonces ``A^TA\in\mathbb{R}^{m\times m}`` es simétrica y definida positiva (sdp) y el sistema tiene solución única ``x=(A^TA)^{-1}A^Tb``. Por este motivo, ``(A^TA)^{-1}A^T`` se conoce usualmente como *pseudoinversa* de ``A``. En Julia se calcula con `pinv`.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000000e
md"""
## Matrices ortogonales
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000000f
md"""
**Definición 3.1.** Se dice que una matriz ``Q\in\mathbb{R}^{n\times m}`` es *ortogonal* si ``Q^TQ=I``.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000010
md"""
**Observación 3.2.** Vale la pena observar lo siguiente:

- Si ``Q`` es ortogonal, la entrada ``ij`` de ``Q^TQ`` es ``\delta_{ij}``, por lo que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000011
md"""
```math
\delta_{ij} = (Q^TQ)_{ij} = \text{fila}_i(Q^T)\,\text{col}_j(Q) = \text{col}_i(Q)^T\text{col}_j(Q) = \text{col}_i(Q)\cdot\text{col}_j(Q).
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000012
md"""
Por lo tanto, las columnas de ``Q`` forman un conjunto *ortonormal*, y en consecuencia las columnas de ``Q`` son linealmente independientes. Como son ``m`` columnas li. en ``\mathbb{R}^n`` necesariamente ``n\geq m``.

- Si ``P`` y ``Q`` son matrices ortogonales y el número de columnas de ``P`` coincide con el número de filas de ``Q``, entonces ``PQ`` es ortogonal, pues
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000013
md"""
```math
(PQ)^TPQ = Q^TP^TPQ = Q^TIQ = Q^TQ = I.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000014
md"""
- Si ``Q`` es ortogonal, entonces ``\|Qx\|_2=\|x\|_2`` y ``(Qx)\cdot(Qy)=x\cdot y``. En efecto
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000015
md"""
```math
(Qx)\cdot(Qy) = (Qx)^T(Qy) = x^TQ^TQy = x^TIy = x^Ty.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000016
md"""
Recordemos que queremos encontrar una manera de minimizar ``\|Ax-b\|_2`` para ``A\in\mathbb{R}^{n\times m}`` y ``b\in\mathbb{R}^n`` con ``n\geq m``. La idea es hallar una manera mejor que la de recurrir a las ecuaciones normales.

Observemos el siguiente razonamiento. Por un lado, si ``Q\in\mathbb{R}^{n\times n}`` es ortogonal, el problema equivale a minimizar ``\|Q^T(Ax-b)\|_2^2`` (como ``Q`` es cuadrada, ``Q^T`` también es ortogonal).

Supongamos que es posible encontrar ``Q\in\mathbb{R}^{n\times n}`` ortogonal tal que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000017
md"""
```math
Q^TA = R = \begin{bmatrix}R_1\\0\end{bmatrix}\begin{matrix}m\times m\\(n-m)\times m\end{matrix}
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000018
md"""
con ``R_1\in\mathbb{R}^{m\times m}`` triangular superior. Luego,
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000019
md"""
```math
\|Q^T(Ax-b)\|_2^2 = \|Q^TAx-Q^Tb\|_2^2 = \|Rx-c\|_2^2,
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000001a
md"""
donde ``c=Q^Tb=\begin{pmatrix}c_1\\c_2\end{pmatrix}`` con ``c_1\in\mathbb{R}^m`` y ``c_2\in\mathbb{R}^{n-m}``. Además, ``Rx=\begin{pmatrix}R_1x\\0\end{pmatrix}``. Por lo tanto, queremos minimizar
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000001b
md"""
```math
\|Q^T(Ax-b)\|_2^2 = \|Rx-c\|_2^2 = \|R_1x-c_1\|_2^2+\|c_2\|_2^2.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000001c
md"""
El mínimo se logra tomando ``x`` la solución del sistema triangular ``R_1x=c_1`` y el mínimo da ``\|c_2\|_2^2``.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000001d
md"""
**Conclusión:**

Dada ``A\in\mathbb{R}^{n\times m}``, (``n\geq m``), ``b\in\mathbb{R}^n``, si existe ``Q\in\mathbb{R}^{n\times n}`` ortogonal tal que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000001e
md"""
```math
Q^TA = R = \begin{bmatrix}R_1\\0\end{bmatrix}
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000001f
md"""
con ``R_1`` triangular superior no singular, entonces, el ``x`` que minimiza ``\|Ax-b\|_2^2`` es la solución de
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000020
md"""
```math
R_1x = (Q^Tb)_{1:m}
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000021
md"""
y además,
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000022
md"""
```math
\|Ax-b\|_2^2 = \|(Q^Tb)_{m+1:n}\|_2^2.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000023
md"""
**Observación 3.3.** Observemos que a los efectos de hallar la solución ``x`` del problema de mínimos cuadrados, basta con obtener **las primeras ``m`` columnas de ``Q``** y **las primeras ``m`` filas de ``R``** (💡). Si llamamos ``Q`` y ``R`` a estas sub-matrices de las ``Q`` y ``R`` anteriores, el problema consiste en hallar ``Q\in\mathbb{R}^{n\times m}`` ortogonal y ``R\in\mathbb{R}^{m\times m}`` triangular superior tales que ``A=QR``. A esto nos abocaremos en lo que sigue.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000024
md"""
**Observación 3.4.** Planteando el sistema normal hace falta primero efectuar el producto ``A^TA``, luego factorizar (vía Cholesky, por ejemplo; ver la Sección 2.7) y luego resolver dos sistemas triangulares. Utilizando la factorización QR de la Observación 3.3 hay que obtener la factorización, y luego resolver ``Rx=Q^Tb``. ¿Qué es más costoso computacionalmente?

Para poder responder esta pregunta debemos primero saber el costo de hallar la factorización QR y saber cuándo es posible.

Observemos que si ``Q\in\mathbb{R}^{n\times m}`` es ortogonal y ``R\in\mathbb{R}^{m\times m}`` triangular superior tales que ``A=QR``, entonces, teniendo en cuenta que ``R`` es triangular superior:
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000025
md"""
```math
\text{col}_1(A) = Q\,\text{col}_1(R) = Q\begin{bmatrix}r_{11}\\0\\\vdots\\0\end{bmatrix} = r_{11}\,\text{col}_1(Q),
```
```math
\text{col}_2(A) = Q\,\text{col}_2(R) = Q\begin{bmatrix}r_{12}\\r_{22}\\0\\\vdots\\0\end{bmatrix} = r_{12}\,\text{col}_1(Q)+r_{22}\,\text{col}_2(Q),
```
``\vdots``
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000026
md"""
Es decir, la columna ``i``-ésima de ``Q`` es ortogonal a las columnas ``1,\dots,i-1`` de ``Q``, y es combinación lineal de las primeras ``i`` columnas de ``A``. Esto tiene una gran reminiscencia de la ortogonalización de Gram-Schmidt. Veamos que este procedimiento nos conduce a una factorización QR.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000027
md"""
## Ortogonalización de Gram-Schmidt
"""

# ╔═╡ 058d902c-5c99-42ab-bde6-f6d569b22794
"""
	gram_schmidth_legible(M::AbstractMatrix)

Implementación legible.
```
function gram_schmidth_legible(M::AbstractMatrix)
	A = similar(M) #buffer de retorno
	nrow, ncol = size(M)
	V = Vector{eltype(M)}(undef, nrow) #buffer para operar por columnas.

	A[:, 1] = M[:, 1] / norm(M[:, 1]) # normaliza primer elemento

	for col in 2:ncol
		V .= M[:, col] - A[:, 1:col-1] * (A[:, 1:col-1]' * M[:, col])
		A[:, col] = normalize(V)
	end

	return A
end
```
"""
function gram_schmidth_legible(M::AbstractMatrix)
	A = similar(M) #buffer de retorno
	nrow, ncol = size(M)
	V = Vector{eltype(M)}(undef, nrow) #buffer para operar por columnas.

	A[:, 1] = M[:, 1] / norm(M[:, 1]) # normaliza primer elemento

	for col in 2:ncol
		V .= M[:, col] - A[:, 1:col-1] * (A[:, 1:col-1]' * M[:, col])
		A[:, col] = normalize(V)
	end

	return A
end

# ╔═╡ 9927ae15-19e3-45b6-aa0f-ac00b6f866ae
"""
	gram_schmidth(M::AbstractMatrix)

Implementación más eficiente.
```
function gram_schmidth!(M::AbstractMatrix)
	rows_qty, cols_qty = size(M)
	
	M[:, 1] ./= BLAS.nrm2(rows_qty, M, 1)

	@inbounds @views for nc in 2:cols_qty
		for ncc in 1:nc-1
			BLAS.axpy!( #αX + Y !-> Y. (a, X, Y)
				-BLAS.dot( # Z ⋅ W. (nrows, Z, zstep, W, wstep)
					rows_qty, M[:,nc], 1, M[:,ncc], 1), #-> α
					M[:, ncc], #-> X
					M[:, nc]) #-> Y
		end
		
		normalize!(M[:, nc])
	end
end
```
"""
function gram_schmidth!(M::AbstractMatrix)
	rows_qty, cols_qty = size(M)

	M[:, 1] ./= BLAS.nrm2(rows_qty, M, 1)

	@inbounds @views for nc in 2:cols_qty
		for ncc in 1:nc-1
			BLAS.axpy!( #αX + Y !-> Y. (a, X, Y)
				-BLAS.dot( # Z ⋅ W. (nrows, Z, zstep, W, wstep)
					rows_qty, M[:,nc], 1, M[:,ncc], 1), #-> α
					M[:, ncc], #-> X
					M[:, nc]) #-> Y
		end

		normalize!(M[:, nc])
	end
end

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000028
md"""
Repasemos el procedimiento de ortogonalización de Gram-Schmidt. Dada una familia de ``m`` vectores ``v_1,v_2,\dots,v_m\in\mathbb{R}^n``, linealmente independientes, el procedimiento es

- ``u_1 = v_1/\|v_1\|``, normalizamos
- ``w_2 = v_2-(v_2\cdot u_1)u_1``, restamos a ``v_2`` la componente en la dirección de ``u_1`` (proyectamos en el subespacio ortogonal a ``u_1``)
- ``u_2 = w_2/\|w_2\|``, normalizamos
- ``w_3 = v_3-(v_3\cdot u_1)u_1-(v_3\cdot u_2)u_2``, restamos a ``v_3`` las componentes en las direcciones de ``u_1`` y ``u_2`` (proyectamos en el subespacio ortogonal a ``u_1`` y ``u_2``)
- ``u_3 = w_3/\|w_3\|``, normalizamos
- ``\vdots``

De esta manera, resulta que ``u_i`` es una combinación lineal de ``v_1,\dots,v_i``, y el factor de ``v_i`` es no nulo. Así, resulta que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000029
md"""
```math
\underbrace{\begin{bmatrix}v_1&v_2&\dots&v_m\end{bmatrix}}_{A} = \underbrace{\begin{bmatrix}u_1&u_2&\dots&u_m\end{bmatrix}}_{Q_1}\underbrace{\begin{bmatrix}\times&\times&\times&\dots&\times\\0&\times&\times&\dots&\times\\0&0&\times&\dots&\times\\\vdots&\vdots&\ddots&\ddots&\vdots\\0&0&\dots&0&\times\end{bmatrix}}_{R_1},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000002a
md"""
y los vectores ``u_i`` son ortonormales.

Hemos demostrado el siguiente teorema:
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000002b
md"""
**Teorema 3.5.** *Si las columnas de ``A\in\mathbb{R}^{n\times m}`` son linealmente independientes, entonces existe una descomposición QR del siguiente tipo:*

``Q_1\in\mathbb{R}^{n\times m}`` *es ortogonal*,

``R_1\in\mathbb{R}^{m\times m}`` *es triangular superior*,

``A=Q_1R_1``.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000002c
md"""
**Observación 3.6.** Bajo las mismas condiciones sobre ``A`` existe una descomposición QR del siguiente tipo:

``Q\in\mathbb{R}^{n\times n}`` *es ortogonal*,

``R = \begin{bmatrix}R_1\\0\end{bmatrix} \in\mathbb{R}^{n\times m}`` *con* ``R_1\in\mathbb{R}^{m\times m}`` *triangular superior*,

``A=QR``.

La matriz ``R_1`` de esta descomposición es la del teorema, y la ``Q`` se obtiene *completando* una base ortonormal de ``\mathbb{R}^n`` a partir de las columnas de ``Q_1``.

Pero existe otra manera de obtener la factorización QR, que es muy interesante y más económica, utilizando las *matrices de Householder*.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000002d
md"""
## Matrices de Householder
"""

# ╔═╡ c3d4e5f6-9999-4a1b-8c2d-000000000005
md"""
**Motivación (no está en el apunte, viene de las diapositivas de la clase).** Antes de definir la matriz de Householder, conviene ver primero una matriz de *proyección*: dado ``w\in\mathbb{R}^n`` con ``\|w\|_2=1``, la matriz
```math
P = I - ww^T
```
proyecta cualquier vector ``x`` sobre el hiperplano ``\{v : w^Tv=0\}`` (ortogonal a ``w``).

Si en cambio queremos *reflejar* ``x`` respecto de ese hiperplano (en vez de proyectarlo), hay que restar el doble de la componente en la dirección de ``w``:
```math
H = I - 2ww^T.
```

Esta ``H`` es exactamente la matriz de Householder que define el apunte, tomando ``u=\sqrt2\,w`` (así ``\|u\|_2=\sqrt2`` y ``H=I-uu^T``).
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000002e
md"""
El material de esta sección fue extraído de [St, Lecture 8].

**Definición 3.7.** Una transformación o matriz de Householder es una matriz de la forma
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000002f
md"""
```math
H = I - uu^T,
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000030
md"""
con ``u\in\mathbb{R}^n`` tal que ``\|u\|_2=\sqrt2``.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000031
md"""
**Figura 3.1** (descripción, no reproducida gráficamente). La matriz de Householder ``H=I-uu^T`` (con ``\|u\|_2=\sqrt2``) refleja cada vector ``x`` respecto del hiperplano ``\{v:u^Tv=0\}``, ortogonal a ``u``: la componente de ``x`` en la dirección de ``u`` cambia de signo y la componente sobre el hiperplano permanece igual.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000032
md"""
**Observación 3.8.** Vale la pena observar lo siguiente, para una matriz ``H=I-uu^T`` de Householder:

- ``H`` es simétrica: ``H^T=H`` (inmediato)

- ``H`` es ortogonal. En efecto
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000033
md"""
```math
H^TH = HH = (I-uu^T)(I-uu^T) = I-2uu^T+u\underbrace{u^Tu}_{\|u\|_2^2}u^T = I-2uu^T+2uu^T = I.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000034
md"""
- Si ``X\in\mathbb{R}^{n\times k}``, hay dos maneras (al menos) de calcular ``HX``:

  - Ensamblar ``H`` (``\frac{n(n+1)}{2}+n\cong\frac{n^2}{2}`` operaciones), y luego hacer el producto ``HX`` que lleva ``2n^2k`` operaciones. Totalizando ``\cong2n^2k`` operaciones.

  - Observar lo siguiente:
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000035
md"""
```math
HX = (I-uu^T)X = X-uu^TX = X-u(u^TX).
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000036
md"""
Si hacemos primero ``u^TX`` nos lleva ``2nk`` operaciones, luego multiplicar ``u`` por lo obtenido ``u(u^TX)`` nos lleva ``nk`` operaciones más. Finalmente, hacer ``X`` menos ese resultado nos lleva ``nk`` operaciones, totalizando ¡solamente ``4nk`` operaciones!

La idea es que la matriz ``Q`` sea una composición o producto de matrices de Householder.
Veamos cómo puede construirse una transformación de Householder con ciertas propiedades.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000037
md"""
**Teorema 3.9.** *Sea ``x`` un vector unitario (``\|x\|_2=1``) y definamos*
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000038
md"""
```math
u = \frac{1}{\sqrt{1\pm x_1}}(x\pm e_1),
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000039
md"""
*donde el signo en ``\pm`` se elige de manera que ``1\pm x_1>0``. Entonces ``H=I-uu^T`` es de Householder y ``Hx=\mp e_1``.*
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000003a
md"""
**Figura 3.2** (descripción, no reproducida gráficamente). Ilustración del Teorema 3.9 en ``\mathbb{R}^2``: dado ``x`` con ``\|x\|_2=1`` y ``x_1>0`` (tomando el signo ``+``), la reflexión ``H=I-uu^T`` con ``u=\frac{1}{\sqrt{1+x_1}}(x+e_1)`` lleva ``x`` al eje generado por ``e_1``, resultando ``Hx=-e_1``. El hiperplano de reflexión biseca el ángulo entre ``x`` y ``-e_1``.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000003b
md"""
**Observación 3.10.** Antes de demostrarlo observemos para qué puede servir. Si ``x`` es la primera columna de una matriz ``A``, y ``H`` se construye como en el teorema, entonces
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000003c
md"""
```math
HA = \begin{bmatrix}\mp1&\times&\dots&\times\\0&\times&\dots&\times\\\vdots&\vdots& &\vdots\\0&\times&\dots&\times\end{bmatrix}, \quad\text{ó}\quad \text{col}_1(HA)=\mp e_1,
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000003d
md"""
donde ``\times`` representa números no necesariamente nulos.

**Demostración.** Veamos primero que ``\|u\|_2=\sqrt2``:
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000003e
md"""
```math
\|u\|_2^2 = u^Tu = \frac{1}{1\pm x_1}\big((x\pm e_1)^T(x\pm e_1)\big) = \frac{1}{1\pm x_1}\big(x^Tx\pm2e_1^Tx+e_1^Te_1\big) = \frac{1}{1\pm x_1}(1\pm2x_1+1) = 2\frac{1\pm x_1}{1\pm x_1} = 2.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000003f
md"""
Veamos ahora cuál es el efecto de ``H`` sobre ``x``: Por un lado,
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000040
md"""
```math
Hx = (I-uu^T)x = x-uu^Tx = x-u(u^Tx),
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000041
md"""
y por otro lado
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000042
md"""
```math
u^Tx = \frac{1}{\sqrt{1\pm x_1}}(x\pm e_1)^Tx = \frac{1}{\sqrt{1\pm x_1}}(x^Tx\pm e_1^Tx) = \frac{1\pm x_1}{\sqrt{1\pm x_1}} = \sqrt{1\pm x_1}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000043
md"""
Por lo tanto
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000044
md"""
```math
Hx = x-u(u^Tx) = x-\sqrt{1\pm x_1}\,u = x-\sqrt{1\pm x_1}\frac{1}{\sqrt{1\pm x_1}}(x\pm e_1) = x-(x\pm e_1) = \mp e_1.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000045
md"""
**Observación 3.11.** Conviene elegir el signo de manera que ``\sqrt{1\pm x_1}`` sea lo más grande posible, pues dividirá otros números. Observemos que ``\|x\|_2=1`` implica ``-1\leq x_1\leq1``, y que:

si ``x_1`` está cerca de ``1`` ``\implies`` ``1-x_1`` estará cerca de cero, y

si ``x_1`` está cerca de ``-1`` ``\implies`` ``1+x_1`` estará cerca de cero.

Para evitar inconvenientes, y para que los divisores que aparezcan sean lo más grandes posibles, elegiremos el signo positivo cuando ``x_1>0`` y negativo cuando ``x_1<0``.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000046
md"""
**Observación 3.12.** Si ``\|x\|_2\neq1``, entonces tomamos
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000047
md"""
```math
u = \frac{1}{\sqrt{1\pm\frac{x_1}{\|x\|_2}}}\Big(\frac{x}{\|x\|_2}\pm e_1\Big),
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000048
md"""
y con ``H=I-uu^T`` resulta
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000049
md"""
```math
Hx = \mp\|x\|_2e_1.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000004a
md"""
A partir del teorema y las observaciones llegamos al siguiente teorema.
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000004b
md"""
**Teorema 3.13.** *Sea ``A\in\mathbb{R}^{n\times m}`` con columnas linealmente independientes. Entonces existen ``Q\in\mathbb{R}^{n\times n}`` ortogonal y ``R\in\mathbb{R}^{n\times m}``, con ``R=\begin{bmatrix}R_1\\0\end{bmatrix}`` y ``R_1\in\mathbb{R}^{m\times m}`` triangular superior no singular, tales que ``A=QR``.*

**Demostración.** Haremos la demostración por inducción sobre ``m``. Si ``m=1``, ``A=x\in\mathbb{R}^n\setminus\{0\}`` (para cualquier ``n\geq1``). Tomamos ``\sigma=1`` si ``x_1\geq0``, o ``\sigma=-1`` si ``x_1<0``, y
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000004c
md"""
```math
u = \frac{1}{\sqrt{1+\sigma\frac{x_1}{\|x\|_2}}}\Big(\frac{x}{\|x\|_2}+\sigma e_1\Big), \qquad H = I-uu^T.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000004d
md"""
Luego, por el teorema anterior y las observaciones
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000004e
md"""
```math
HA = \begin{bmatrix}-\sigma\|x\|_2\\0\\\vdots\\0\end{bmatrix} =: R, \qquad\implies\qquad A = \underbrace{H}_{Q}R.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000004f
md"""
Supongamos válido el teorema para ``m`` (y todo ``n\geq m``) y demostrémoslo para ``m+1``.
Sea ``A\in\mathbb{R}^{n\times(m+1)}`` con columnas linealmente independientes. Sea ``x=\text{col}_1(A)\neq0``, tomemos ``\sigma=1`` si ``x_1\geq0``, o ``\sigma=-1`` si ``x_1<0`` y
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000050
md"""
```math
u = \frac{1}{\sqrt{1+\sigma\frac{x_1}{\|x\|_2}}}\Big(\frac{x}{\|x\|_2}+\sigma e_1\Big), \qquad H = I-uu^T.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000051
md"""
Luego
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000052
md"""
```math
HA = \begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat A\end{bmatrix}, \qquad\text{con } r\in\mathbb{R}^m \text{ y } \hat A\in\mathbb{R}^{(n-1)\times m}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000053
md"""
Las columnas de ``\hat A`` son linealmente independientes, pues si no ``A`` no tendría ``m+1`` columnas linealmente independientes (demostrado abajo 👉🍎). Por la hipótesis inductiva, existen ``\hat Q\in\mathbb{R}^{(n-1)\times(n-1)}`` ortogonal y ``\hat R=\begin{bmatrix}\hat R_1\\0\end{bmatrix}\in\mathbb{R}^{(n-1)\times m}`` con ``\hat R_1\in\mathbb{R}^{m\times m}`` triangular superior tales que ``\hat A=\hat Q\hat R``. Luego
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000054
md"""
```math
HA = \begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat Q\hat R\end{bmatrix} = \begin{bmatrix}1&0^T\\0&\hat Q\end{bmatrix}\begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat R\end{bmatrix},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000055
md"""
y por lo tanto
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000056
md"""
```math
A = H\underbrace{\begin{bmatrix}1&0^T\\0&\hat Q\end{bmatrix}}_{Q}\underbrace{\begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat R\end{bmatrix}}_{R}.
```
"""

# ╔═╡ 007f43a6-15ec-45fa-be4d-3a25261abf3f
md"""
**Proposición: 🍎** Si las columnas de ``A`` son *li* entonces las columnas de ``\hat{A}`` son *li*.
"""

# ╔═╡ 95aa7cc0-4534-46d0-91c3-33a721d94435
md"""
**Demostración:**

Supongamos que las columnas de ``\hat A`` no son *li*. Entonces, existe ``\alpha\in\mathbb{R}^m`` no nulo tal que ``\hat A\alpha = 0``. Tomemos
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000008a
md"""
```math
\beta_1 := \frac{r^T\alpha}{\sigma\|x\|_2}
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000008b
md"""
(bien definido pues ``\sigma\|x\|_2\neq0``, ya que ``x=\text{col}_1(A)\neq0``) y sea ``\beta=(\beta_1,\alpha)\in\mathbb{R}^{m+1}``, que es no nulo pues ``\alpha\neq0``. Entonces
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000008c
md"""
```math
HA\beta = \begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat A\end{bmatrix}\begin{bmatrix}\beta_1\\\alpha\end{bmatrix} = \begin{bmatrix}-\sigma\|x\|_2\beta_1+r^T\alpha\\\hat A\alpha\end{bmatrix} = \begin{bmatrix}0\\0\end{bmatrix},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000008d
md"""
donde la primera fila se anula por la elección de ``\beta_1``, y la segunda por ser ``\hat A\alpha=0``.

Como ``H`` es una transformación de Householder, es ortogonal (``H=H^{-1}``), y por lo tanto invertible. Luego ``HA\beta=0 \iff A\beta=0``. Así, ``\beta\neq0`` satisface ``A\beta=0``, es decir, las columnas de ``A`` no serían *li*, contradiciendo la hipótesis.
"""

# ╔═╡ a5c99dc0-b35a-4fcf-9547-32a7c7fb40cc
md"""
-----------------------
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000057
md"""
**Implementación.**

Para implementar el algoritmo seguiremos los pasos de la demostración del teorema.

Dada ``A\in\mathbb{R}^{n\times m}`` (con sus columnas linealmente independientes), hacemos (aquí y en lo que sigue usamos la convención ``\text{sgn}(x_1)=1`` si ``x_1\geq0`` y ``\text{sgn}(x_1)=-1`` si ``x_1<0``)
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000058
md"""
```math
x=\text{col}_1(A), \quad \sigma=\text{sgn}(x_1), \quad u_1=\frac{1}{\sqrt{1+\sigma\frac{x_1}{\|x\|_2}}}\Big(\frac{x}{\|x\|_2}+\sigma e_1\Big), \quad H_1=I-u_1u_1^T.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000059
md"""
Luego, sabemos que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000005a
md"""
```math
H_1A = \begin{bmatrix}\rho_1&r_1^T\\0&A_1\end{bmatrix}, \quad\text{con } \rho_1=-\sigma\|x\|_2, \text{ y } \begin{bmatrix}r_1^T\\A_1\end{bmatrix}=H_1A_{:,2:m}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000005b
md"""
Y también
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000005c
md"""
```math
A = H_1\begin{bmatrix}\rho_1&r_1^T\\0&A_1\end{bmatrix}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000005d
md"""
Ahora hacemos
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000005e
md"""
```math
x = \text{col}_1(A_1) \in \mathbb{R}^{n-1}, \quad \sigma=\text{sgn}(x_1), \quad u_2 = \frac{1}{\sqrt{1+\sigma\frac{x_1}{\|x\|_2}}}\Big(\frac{x}{\|x\|_2}+\sigma e_1\Big), \quad H_2 = I-u_2u_2^T \in \mathbb{R}^{(n-1)\times(n-1)}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000005f
md"""
Entonces
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000060
md"""
```math
H_2A_1 = \begin{bmatrix}\rho_2&r_2^T\\0&A_2\end{bmatrix}, \quad\text{con } \rho_2=-\sigma\|x\|_2, \text{ y } \begin{bmatrix}r_2^T\\A_2\end{bmatrix}=H_2(A_1)_{:,2:m-1}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000061
md"""
Y también
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000062
md"""
```math
A_1 = H_2\begin{bmatrix}\rho_2&r_2^T\\0&A_2\end{bmatrix},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000063
md"""
por lo que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000064
md"""
```math
A = H_1\begin{bmatrix}\rho_1&r_1^T\\0&A_1\end{bmatrix} = H_1\begin{bmatrix}\rho_1&r_1^T\\0&H_2\begin{bmatrix}\rho_2&r_2^T\\0&A_2\end{bmatrix}\end{bmatrix} = H_1\begin{bmatrix}1&0^T\\0&H_2\end{bmatrix}\begin{bmatrix}\rho_1&&r_1^T\\0&\rho_2&r_2^T\\0&0&A_2\end{bmatrix}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000065
md"""
Ahora hacemos
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000066
md"""
```math
x=\text{col}_1(A_2)\in\mathbb{R}^{n-2}, \quad \sigma=\text{sgn}(x_1), \quad u_3=\frac{1}{\sqrt{1+\sigma\frac{x_1}{\|x\|_2}}}\Big(\frac{x}{\|x\|_2}+\sigma e_1\Big), \quad H_3=I-u_3u_3^T \in \mathbb{R}^{(n-2)\times(n-2)}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000067
md"""
Entonces
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000068
md"""
```math
H_3A_2 = \begin{bmatrix}\rho_3&r_3^T\\0&A_3\end{bmatrix}, \quad\text{con } \rho_3=-\sigma\|x\|_2, \text{ y } \begin{bmatrix}r_3^T\\A_3\end{bmatrix}=H_3(A_2)_{:,2:m-2}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000069
md"""
Y también
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000006a
md"""
```math
A_2 = H_3\begin{bmatrix}\rho_3&r_3^T\\0&A_3\end{bmatrix},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000006b
md"""
por lo que
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000006c
md"""
```math
A = H_1\begin{bmatrix}1&0^T\\0&H_2\end{bmatrix}\begin{bmatrix}1&0&0^T\\0&1&0^T\\0&0&H_3\end{bmatrix}\begin{bmatrix}\rho_1&&&r_1^T\\0&\rho_2&&r_2^T\\0&0&\rho_3&r_3^T\\0&0&0&A_3\end{bmatrix}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000006d
md"""
Así sucesivamente ....

El algoritmo para MATLAB/OCTAVE (tal como aparece en el apunte) es el siguiente; a continuación se lo traduce a Julia.

```octave
function [U,R] = qrhouse(A)
% function [U,R] = qrhouse(A)
%
% Descomposicion QR de A
% La matriz R tiene el mismo tamaño que A pero es triangular superior
% La matriz Q es cuadrada y ortogonal
%   A = Q R
%
% Pero no se calcula Q, sino que se devuelven los vectores u que
% generan las transformaciones de Householder que conforman Q
%
%   Q' = Qk ... Q2 Q1
%
% con Q1 = H1 = eye(n) - u1 u1'
%     Q2 = [1 zeros(1,n-1); zeros(n-1,1) H2]   y   H2 = eye(n-1) - u2 u2'
%     Q3 = [eye(2) zeros(2,n-2); zeros(n-2,2) H3]
%           y H3  = eye(n-2) - u3 u3'
%       ...
%   uk = U(1:n-k+1, k)
%
% ver qrsolve

[n, m] = size(A);
R = A;
U = zeros(n, m);

for i = 1:m
    normx = norm(R(i:n, i));
    u = R(i:n, i)/normx;
    if (u(1) >= 0)
        sg = 1;
    else
        sg = -1;
    end
    u(1) = u(1) + sg;
    u = u / sqrt(1+sg*R(i,i)/normx);

    U(1:n-i+1, i) = u;
    R(i,i) = -sg*normx;
    R(i+1:n,i) = 0;
    R(i:n,i+1:m) = R(i:n,i+1:m) - u*(u'*R(i:n,i+1:m));
end
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000006e
"""
    qrhouse_libro(A)

Descomposición QR de ``A`` mediante transformaciones de Householder.
La matriz ``R`` tiene el mismo tamaño que ``A`` pero es triangular superior.
No se calcula ``Q`` explícitamente: se devuelven en las columnas de ``U`` los
vectores ``u`` que generan las transformaciones de Householder que la
conforman (``Q^T = H_m ⋯ H_2H_1``, con ``H_k`` construida a partir de
``U[1:n-k+1, k]``, ver `qrsolve_libro`).

```julia
function qrhouse_libro(A)
    n, m = size(A)
    R = copy(A)
    U = zeros(n, m)

    for i in 1:m
        normx = norm(R[i:n, i])
        u = R[i:n, i] / normx
        sg = u[1] >= 0 ? 1 : -1
        u[1] = u[1] + sg
        u = u / sqrt(1 + sg * R[i, i] / normx)

        U[1:n-i+1, i] = u
        R[i, i] = -sg * normx
        R[i+1:n, i] .= 0
        R[i:n, i+1:m] .-= u * (u' * R[i:n, i+1:m])
    end

    return U, R
end
```
"""
function qrhouse_libro(A::AbstractMatrix{T}) where T
    n, m = size(A)
    R = copy(A)
    U = zeros(T, n, m)

    for i in 1:m
        normx = norm(R[i:n, i])
        u = R[i:n, i] / normx
        sg = u[1] >= 0 ? 1 : -1
        u[1] = u[1] + sg
        u = u / sqrt(1 + sg * R[i, i] / normx)

        U[1:n-i+1, i] = u
        R[i, i] = -sg * normx
        R[i+1:n, i] .= 0
        R[i:n, i+1:m] .-= u * (u' * R[i:n, i+1:m])
    end

    return U, R
end

# ╔═╡ d37ade2c-3528-4aa1-91bb-3294a9a4c849
"""
    qrhouse_con_cami(A)
```julia
function qrhouse_con_cami(A)
    R = copy(A)
    Q = similar(A)

    nrow, ncol = size(A)

    function vectoru(x, nx, σ)
        a = 1 / sqrt(1 + σ * x[1] / nx)
        v = x ./ nx
        v[1] += σ

        return a * v
    end

    for nc in 1:ncol
        x = R[nc:nrow, nc]
        norma_x = norm(x)

        σ = x[1] >= 0 ? 1.0 : -1.0
        ρ = -σ * norma_x
        u = vectoru(x, norma_x, σ)
        R[nc, nc] = ρ
        R[nc+1:nrow,nc] .= 0
        R[nc:nrow,nc+1:ncol] .-= u * (u' * R[nc:nrow,nc+1:ncol])


        if nc == 1
            Q .= I - u * u'
        else
            Q[1:nrow, nc:ncol] .-= (Q[1:nrow, nc:ncol] * u) * u'
        end
    end

    return Q, R
end
```
"""
function qrhouse_con_cami(A)
    R = copy(A)
    Q = similar(A)

    nrow, ncol = size(A)

    function vectoru(x, nx, σ)
        a = 1 / sqrt(1 + σ * x[1] / nx)
        v = x ./ nx
        v[1] += σ

        return a * v
    end

    for nc in 1:ncol
        x = R[nc:nrow, nc]
        norma_x = norm(x)

        σ = x[1] >= 0 ? 1.0 : -1.0
        ρ = -σ * norma_x
        u = vectoru(x, norma_x, σ)
        R[nc, nc] = ρ
        R[nc+1:nrow,nc] .= 0
        R[nc:nrow,nc+1:ncol] .-= u * (u' * R[nc:nrow,nc+1:ncol])


        if nc == 1
            Q .= I - u * u'
        else
            Q[1:nrow, nc:ncol] .-= (Q[1:nrow, nc:ncol] * u) * u'
        end
    end

    return Q, R
end

# ╔═╡ 22ea694a-3fc1-4752-bac3-77b3b420191d
"""
    qrhouse(A)

Implmementacion optimizada
```julia
function qrhouse(A)
    nrow, ncol = size(A)
    R = copy(A)
    Q = Matrix{eltype(A)}(I, nrow, nrow)

    vi = Vector{Float64}(undef, nrow)
    wr = Vector{Float64}(undef, ncol) #para actualizar R
    wq = Vector{Float64}(undef, nrow) #para actualizar Q

    function vectoru!(u, x, nx, σ)
        @. u = x / nx
        u[1] += σ # u + σe₁
        a = 1 / sqrt(1 + σ * x[1] / nx)
        u .*= a
        return u
    end
    
    @inbounds @views for nc in 1:ncol
        x = R[nc:nrow, nc]
        norma_x = norm(x)

        σ = x[1] >= 0 ? 1.0 : -1.0
        ρ = -σ * norma_x
        len = nrow - nc + 1
        u = vi[1:len]
        vectoru!(u, x, norma_x, σ)

        R[nc, nc] = ρ
        R[nc+1:nrow,nc] .= 0

        mcols = ncol - nc
        if mcols > 0
            Rsub = R[nc:nrow, nc+1:ncol]
            w = wr[1:mcols]
            mul!(w, Rsub', u)
            mul!(Rsub, u, w', -1.0, 1.0)
        end

        Qsub = Q[1:nrow, nc:ncol]
        mul!(wq, Qsub, u)
        mul!(Qsub, wq, u', -1.0, 1.0)
    end

    return Q, R
end
```
"""
function qrhouse(A)
    nrow, ncol = size(A)
    R = copy(A)
    Q = Matrix{eltype(A)}(I, nrow, nrow)

    vi = Vector{Float64}(undef, nrow)
    wr = Vector{Float64}(undef, ncol) #para actualizar R
    wq = Vector{Float64}(undef, nrow) #para actualizar Q

    function vectoru!(u, x, nx, σ)
        @. u = x / nx
        u[1] += σ #en lugar de multiplicar por canonico y sumar.
        a = 1 / sqrt(1 + σ * x[1] / nx)
        u .*= a
        return u
    end
    
    @inbounds @views for nc in 1:ncol
        x = R[nc:nrow, nc]
        norma_x = norm(x)

        σ = x[1] >= 0 ? 1.0 : -1.0
        ρ = -σ * norma_x
        len = nrow - nc + 1
        u = vi[1:len]
        vectoru!(u, x, norma_x, σ)

        R[nc, nc] = ρ
        R[nc+1:nrow,nc] .= 0

        mcols = ncol - nc
        if mcols > 0
            Rsub = R[nc:nrow, nc+1:ncol]
            w = wr[1:mcols]
            mul!(w, Rsub', u)
            mul!(Rsub, u, w', -1.0, true)
        end

        Qsub = Q[1:nrow, nc:ncol]
        mul!(wq, Qsub, u)
        mul!(Qsub, wq, u', -1.0, true)
    end

    return Q, R
end

# ╔═╡ b91efe6d-df9d-43c1-a33e-b944eb1436e1
"""
    urhouse(A)

Implmementacion optimizada
```julia
function urhouse(A)
    nrow, ncol = size(A)
    R = copy(A)
    U = Matrix{eltype(A)}(undef, nrow, nrow)

    vi = Vector{Float64}(undef, nrow)
    wr = Vector{Float64}(undef, ncol) #para actualizar R
    wq = Vector{Float64}(undef, nrow) #para actualizar Q

    function vectoru!(u, x, nx, σ)
        @. u = x / nx
        u[1] += σ # u + σe₁
        a = 1 / sqrt(1 + σ * x[1] / nx)
        u .*= a
        return u
    end
    
    @inbounds @views for nc in 1:ncol
        x = R[nc:nrow, nc]
        norma_x = norm(x)

        σ = x[1] >= 0 ? 1.0 : -1.0
        ρ = -σ * norma_x
        len = nrow - nc + 1
        u = vi[1:len]
        vectoru!(u, x, norma_x, σ)

        R[nc, nc] = ρ
        R[nc+1:nrow,nc] .= 0

        mcols = ncol - nc
        if mcols > 0
            Rsub = R[nc:nrow, nc+1:ncol]
            w = wr[1:mcols]
            mul!(w, Rsub', u)
            mul!(Rsub, u, w', -1.0, 1.0)
        end

        U[1:len,nc] = u
    end

    return U, R
end
```
"""
function urhouse(A)
    nrow, ncol = size(A)
    R = copy(A)
    U = Matrix{eltype(A)}(undef, nrow, nrow)

    vi = Vector{Float64}(undef, nrow)
    wr = Vector{Float64}(undef, ncol)

    function vectoru!(u, x, nx, σ)
        @. u = x / nx
        u[1] += σ # u + σe₁
        a = 1 / sqrt(1 + σ * x[1] / nx)
        u .*= a
        return u
    end

    sgn(x) = x >= 0 ? 1.0 : -1.0
    
    @inbounds @views for nc in 1:ncol
        x = R[nc:nrow, nc]
        norma_x = norm(x)
        σ = sgn(x[1])
        ρ = -σ * norma_x

        len = nrow - nc + 1
        u = vi[1:len]
        vectoru!(u, x, norma_x, σ)
        U[1:len,nc] = u

        R[nc, nc] = ρ
        R[nc+1:nrow,nc] .= 0

        mcols = ncol - nc
        if mcols > 0
            Rsub = R[nc:nrow, nc+1:ncol]
            w = wr[1:mcols]
            mul!(w, Rsub', u)
            mul!(Rsub, u, w', -1.0, 1.0)
        end
    end

    return U, R
end

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000006f
md"""
Esta función se complementa con `qrsolve_libro`, que utiliza la descomposición QR dada por `qrhouse_libro` para resolver un problema de mínimos cuadrados. Utiliza la función `trisup_libro(T,c)`, que resuelve el sistema triangular superior ``Tx=c`` mediante sustitución hacia atrás (retrosustitución), tal como la función `trisup` de los ejercicios del Capítulo 2 a la que hace referencia el apunte.
"""

# ╔═╡ 4a096e7e-9f35-4df1-a85b-2192f61fc965
"""
    ursolve(U, R, b)

```julia
function ursolve(U, R, b)
    nrow, ncol = size(U)
    bb = copy(b)

    @inbounds @views for nc in 1:ncol
        u = U[1:(nrow-nc+1), nc]
        bv = bb[nc:nrow]
        escalar = u' * bv
        axpy!(-escalar, u, bv)
    end

    @views(UpperTriangular(R[1:size(R,2),:]) \\ bb[1:size(R,2)])
end
```
"""
function ursolve(U, R, b)
    nrow, ncol = size(U)
    bb = copy(b)

    @inbounds @views for nc in 1:ncol
        u = U[1:(nrow-nc+1), nc]
        bv = bb[nc:nrow]
        escalar = u' * bv
        axpy!(-escalar, u, bv)
    end

    @views(UpperTriangular(R[1:size(R,2),:]) \ bb[1:size(R,2)])
end

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000072
md"""
## Conteo de operaciones y estabilidad
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000073
md"""
**Conteo de operaciones.**

Contamos las operaciones a orden dominante, en el caso ``n\gg m`` (muchas más filas que columnas), que es el habitual en mínimos cuadrados.

- *Ecuaciones normales* (vía Cholesky). Formar ``A^TA`` cuesta ``\cong nm^2`` operaciones: la matriz es simétrica, y cada una de las ``\cong m^2/2`` entradas de su triángulo superior es un producto interno de dos columnas de ``A`` (``\cong2n`` operaciones). Formar ``A^Tb`` cuesta ``\cong2nm``, la factorización de Cholesky de ``A^TA`` cuesta ``\cong\frac{1}{3}m^3``, y resolver los dos sistemas triangulares ``\cong2m^2``. En total
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000074
md"""
```math
\cong nm^2+\frac{1}{3}m^3\cong nm^2 \qquad(\text{si } n\gg m).
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000075
md"""
- *QR por Gram-Schmidt.* Ortonormalizar las ``m`` columnas de ``A`` cuesta ``\cong2nm^2`` operaciones.

- *QR por Householder.* Aplicar las ``m`` transformaciones ``H_i`` sin ensamblar ``Q`` (usando ``HX=X-u(u^TX)``) cuesta ``\cong2nm^2-\frac{2}{3}m^3\cong2nm^2``; resolver ``R_1x=(Q^Tb)_{1:m}`` y aplicar los ``H_i`` a ``b`` agregan sólo términos de orden inferior.

Resumimos los costos a orden dominante (``n\gg m``):

| Método | Operaciones |
|---|---|
| Ecuaciones normales (Cholesky) | ``nm^2+\frac{1}{3}m^3\cong nm^2`` |
| QR por Gram–Schmidt | ``\cong2nm^2`` |
| QR por Householder | ``2nm^2-\frac{2}{3}m^3\cong2nm^2`` |

Así, las ecuaciones normales requieren aproximadamente *la mitad* de las operaciones que la vía QR. Sin embargo, como veremos enseguida, son menos *estables*. El conteo detallado de cada método queda propuesto como ejercicio al final del capítulo (no se resuelve en este cuaderno; ver Problema 3.5 del apunte).

**Estabilidad y número de condición.**

Para comparar la *precisión* de los métodos conviene introducir el *número de condición* de una matriz. Para una matriz cuadrada e invertible ``B``, el número de condición en norma 2 es ``\kappa_2(B)=\|B\|_2\|B^{-1}\|_2\geq1``; mide cuánto puede amplificarse, en el peor caso, el error relativo de la solución de ``Bx=c`` frente a perturbaciones en los datos (lo estudiaremos en detalle más adelante, en el capítulo de métodos iterativos). Si ``B`` es sdp, resulta ``\kappa_2(B)=\lambda_{\max}/\lambda_{\min}``, el cociente entre su mayor y su menor autovalor.

Como ``A^TA`` es sdp (lo vimos al plantear las ecuaciones normales), podemos medir el condicionamiento del problema de mínimos cuadrados a través de él. Definimos el número de condición de ``A\in\mathbb{R}^{n\times m}`` (de rango ``m``) en norma 2 como
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000076
md"""
```math
\kappa_2(A) = \sqrt{\kappa_2(A^TA)} = \sqrt{\frac{\lambda_{\max}(A^TA)}{\lambda_{\min}(A^TA)}}, \qquad\text{de modo que}\qquad \kappa_2(A^TA) = \kappa_2(A)^2.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000077
md"""
*(Equivalentemente, en términos de los valores singulares ``\sigma_1\geq\dots\geq\sigma_m>0`` de ``A`` (las raíces cuadradas de los autovalores de ``A^TA``), es ``\kappa_2(A)=\sigma_1/\sigma_m=\|A\|_2\|A^+\|_2``, donde ``A^+=(A^TA)^{-1}A^T`` es la pseudoinversa de ``A``. Para ``A`` cuadrada e invertible, ``A^+=A^{-1}`` y se recupera ``\kappa_2(A)=\|A\|_2\|A^{-1}\|_2``.)*

Decimos que un algoritmo es *estable hacia atrás* (*backward stable*) si la solución ``\hat x`` que produce en aritmética de punto flotante es la solución *exacta* de un problema *levemente perturbado*. Para ``\min\|Ax-b\|_2``, esto significa que existen perturbaciones ``\Delta A``, ``\Delta b`` con
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000078
md"""
```math
\frac{\|\Delta A\|}{\|A\|}, \quad \frac{\|\Delta b\|}{\|b\|} \cong \varepsilon_{\text{maq}}
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000079
md"""
(del orden del *épsilon de máquina*, ``\varepsilon_{\text{maq}}\approx2.2\times10^{-16}`` en doble precisión) tales que ``\hat x`` resuelve exactamente ``\min\|(A+\Delta A)x-(b+\Delta b)\|_2``. En palabras: el algoritmo no introduce más error que el inevitable al representar los datos en la máquina.

Conviene distinguir el *error hacia adelante* ``\|\hat x-x\|`` (cuán lejos está la solución calculada de la exacta) del *error hacia atrás* (cuánto hay que perturbar los datos para que ``\hat x`` sea exacta). Ambos se relacionan, a grandes rasgos, por
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000007a
md"""
```math
\text{error hacia adelante} \cong \kappa_2(A) \times (\text{error hacia atrás}).
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000007b
md"""
Por lo tanto, un algoritmo estable hacia atrás produce un error hacia adelante ``\cong\kappa_2(A)\varepsilon_{\text{maq}}``: lo mejor que puede esperarse dado el condicionamiento del problema.

**Observación 3.14.** La factorización QR por Householder es estable hacia atrás: ``Q`` es un producto de reflexiones, que son isometrías (``\|Qx\|_2=\|x\|_2``) y por lo tanto no amplifican los errores de redondeo. En cambio, el método de las ecuaciones normales forma ``A^TA``, cuyo número de condición es ``\kappa_2(A^TA)=\kappa_2(A)^2``. Aun resolviendo ``A^TAx=A^Tb`` de manera perfecta, el error hacia adelante puede ser ``\cong\kappa_2(A)^2\varepsilon_{\text{maq}}``: se pierde *el doble* de dígitos significativos que con la vía QR. Por eso, cuando ``A`` está mal condicionada, se prefiere resolver el problema de mínimos cuadrados mediante la factorización QR (Householder).
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000007c
md"""
## Rotaciones de Givens
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000007d
md"""
Consideraremos en esta sección el caso particular en que queremos hallar la factorización QR de una matriz *upper Hessenberg*. Esta clase de matrices aparece en métodos iterativos para sistemas lineales como el GMRes y también en el cálculo de autovalores.

**Definición 3.15.** Se dice que una matriz ``H`` es *upper Hessenberg* si sus componentes ``h_{ij}=0`` para todo ``i>j+1``. Es decir, las entradas por debajo de la primera subdiagonal son nulas.

Para matrices de tipo *upper Hessenberg*, o simplemente *Hessenberg*, se puede obtener una factorización QR mucho más eficientemente que usando las transformaciones de Householder, utilizando las *rotaciones de Givens*.
"""

# ╔═╡ 33f6e3d7-50a6-4ad2-a272-e74cbbd59ed8
md"""
**Definición 3.16.** Una *rotación de Givens* de ``2\times2`` es una matriz de la forma
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000007e
md"""
```math
G = \begin{bmatrix}\cos(\theta)&-\text{sen}(\theta)\\\text{sen}(\theta)&\cos(\theta)\end{bmatrix}
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-00000000007f
md"""
considerando ``\theta\in[-\pi,\pi]``. Esta matriz es ortogonal y *rota* el vector ``\begin{bmatrix}c\\-s\end{bmatrix}=\begin{bmatrix}\cos(\theta)\\-\text{sen}(\theta)\end{bmatrix}`` hacia el vector ``\begin{bmatrix}1\\0\end{bmatrix}``. Es decir,
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000080
md"""
```math
\text{si } c^2+s^2=1, \qquad G=\begin{bmatrix}c&-s\\s&c\end{bmatrix} \qquad\text{resulta}\qquad G\begin{bmatrix}c\\-s\end{bmatrix} = \begin{bmatrix}1\\0\end{bmatrix}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000081
md"""
Una *rotación de Givens* de ``n\times n`` es una matriz identidad de ``n\times n`` donde se ha reemplazado un bloque de tamaño ``2\times2`` de la diagonal por una matriz de rotación de Givens de ``2\times2``. Es decir, una matriz de la forma:
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000082
md"""
```math
G = \begin{bmatrix}1&0&\dots& & &0\\0&\ddots&\ddots& & &\\ &\ddots&c&-s& &\vdots\\\vdots& &s&c&0& \\ & & &0&1&\ddots\\ & & &\ddots&\ddots&0\\0& & & &0&1\end{bmatrix}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000083
md"""
Usaremos la notación ``G_j(c,s)`` para indicar una matriz de rotación de Givens de ``n\times n`` con la matriz de Givens de ``2\times2`` en las filas y columnas ``j`` y ``j+1``.

La idea es utilizar rotaciones de Givens para reducir matrices de Hessenberg a matrices triangulares superiores, de la siguiente manera: Sea ``H`` una matriz de Hessenberg de ``n\times m`` (``n\geq m``) de rango ``m``. El objetivo es reducirla a triangular superior multiplicándola a izquierda por rotaciones de Givens. Si queremos producir un cero en la posición ``(2,1)`` de la matriz ``H`` tenemos que pensar que estamos rotando el vector ``\begin{bmatrix}h_{11}\\h_{21}\end{bmatrix}``. Por ello, definimos ``G_1=G_1(c_1,s_1)`` con
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000084
md"""
```math
c_1 = \frac{h_{11}}{\sqrt{h_{11}^2+h_{21}^2}} \qquad\text{y}\qquad s_1 = \frac{-h_{21}}{\sqrt{h_{11}^2+h_{21}^2}}.
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000085
md"""
Si multiplicamos ``G_1`` por ``H``, en la primera columna de ``H^1=G_1H`` queda un único elemento no nulo en la posición ``(1,1)``. Vale la pena observar que el producto ``G_1H`` sólo afecta las dos primeras filas de ``H``.
Con la nueva matriz ``H^1``, definimos
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000086
md"""
```math
c_2 = \frac{h_{22}}{\sqrt{h_{22}^2+h_{32}^2}} \qquad\text{y}\qquad s_2 = \frac{-h_{32}}{\sqrt{h_{22}^2+h_{32}^2}},
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000087
md"""
donde ahora ``h_{22}`` y ``h_{32}`` denotan las entradas correspondientes de ``H^1``, y ``G_2=G_2(c_2,s_2)``, y así sucesivamente, anulando una por una las entradas de la primera subdiagonal. Si ``n=m`` la última entrada subdiagonal es la ``(m,m-1)`` y bastan ``p=m-1`` pasos; si ``n>m`` hay que anular también la entrada ``(m+1,m)`` y hacen falta ``p=m`` pasos. Llegamos así a
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000088
md"""
```math
\underbrace{G_pG_{p-1}\dots G_1}_{Q^T} H = R,
```
"""

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000089
md"""
con ``R`` triangular superior y ``Q^T=G_pG_{p-1}\dots G_1`` ortogonal de ``n\times n``.

Observamos también que al aplicar ``G_2`` la primera columna de ``H^1`` no se modifica, al aplicar ``G_3`` no se modifican la primera y la segunda columna, y así sucesivamente.

*Nota:* en este punto, el apunte deja pendiente el cálculo del número de operaciones necesario para obtener ``Q`` y ``R`` con este método ("Queda pendiente hacer un cálculo del número de operaciones necesario para obtener ``Q`` y ``R``."). No se resuelve aquí, tal como está en la fuente.
"""

# ╔═╡ 617150c9-0cb5-4dc8-9b96-e780b5680e62
"""
	rotacion_de_givens(M)
```
function rotacion_de_givens(M)
	nrow, ncol = size(M)
	A = copy(M)

	@views @inbounds @simd for icol in 1:ncol-1
		col = A[icol:icol+1, icol]
		nrm = hypot(col[1], col[2])
		c =  col[1] / nrm
		s = -col[2] / nrm

		for j in icol:ncol #En lugar del prod. matricial para evitar overhead
			a = A[icol, j]
			b = A[icol+1, j]
			A[icol, j] = c*a - s*b
			A[icol+1, j] = s*a + c*b
		end
	end
	
	return A
end
```
"""
function rotacion_de_givens(M)
	nrow, ncol = size(M)
	A = copy(M)

	@views @inbounds @simd for icol in 1:ncol-1
		col = A[icol:icol+1, icol]
		nrm = hypot(col[1], col[2])
		c =  col[1] / nrm
		s = -col[2] / nrm

		for j in icol:ncol #En lugar del producto matricial para evitar overhead
			a = A[icol, j]
			b = A[icol+1, j]
			A[icol, j] = c*a - s*b
			A[icol+1, j] = s*a + c*b
		end
	end
	
	return A
end

# ╔═╡ c32c7715-5bd6-4682-8f14-ed52a5e2b9dd
md"""
## 3.7 Problemas
"""

# ╔═╡ db376a73-bfa2-40d0-a282-074bb4ba8bae
md"""
### 3.1 Ortogonalización de Gram-Schmidt.

#### 3.1.a) 

(*Está en el libro*). Lo importante es que es iterativo, y que no es paralelizable.
"""

# ╔═╡ bd039dc5-e6ce-4a2b-be77-4481ad11a2b5
md"""
#### 3.1.b)

(*También está probado en el libro*)
"""

# ╔═╡ 59b39939-05d9-488b-a5a6-2df6581d05d8
md"""
### 3.2)

Realizar una función `mincuad(A,b)` que resuelva las ecuaciones normales correspondientes al sistema ``Ax=b``. ¿Qué pasa si para un sistema rectangular con más filas que columnas escribimos `A\b`? Verificar y comparar el resultado arrojado por esta última fórmula y el obtenido por la función `mincuad(A,b)`. Utilizando la documentación, indagar el funcionamiento del símbolo `\` y de la función `pinv` e informar brevemente.
"""

# ╔═╡ a3126c2a-b66a-453a-b208-5bb2988a6d8d
"""
	mincuad1(A, b)

```julia
function mincuad1(A, b)
	inv(transpose(A) * A) * transpose(A) * b #<- implementación directa.
end
```
"""
function mincuad1(A, b)
	inv(transpose(A) * A) * transpose(A) * b #<- implementación directa.
end

# ╔═╡ fe793033-46cd-4c77-87e0-69cd82afdaa1
let
	A = rand(1000,3)
	b = rand(1000)

	@benchmark mincuad1($A, $b)
end

# ╔═╡ 84bfa63d-cfe8-4c51-8017-4dee5f2cc906
"""
	mincuad2(A, b)
```julia
function mincuad2(A, b)
	(transpose(A) * A) \\ (transpose(A) * b) # sin invertir
end
```
"""
function mincuad2(A, b)
	(transpose(A) * A) \ (transpose(A) * b)
end

# ╔═╡ f3d59908-9bbe-475b-8c49-38aab198368d
let
	A = rand(1000,3)
	b = rand(1000)

	@benchmark mincuad2($A, $b)
end

# ╔═╡ 04f5f8f5-4f44-45d5-abee-7adb1c1e4ae1
md"""
**Documentación de `\`.**
"""

# ╔═╡ d38c1824-93bc-4ebf-a700-0850d9335df0
"""
"""
\

# ╔═╡ c5d58661-106c-468d-af76-0e35194f35aa
md"""
**Documentación de `pinv`.**
"""

# ╔═╡ 90c6ba41-71d5-4451-8709-e5d2a3453454
"""
"""
pinv

# ╔═╡ 41b40cc5-9d7a-4c78-be4b-666367339dbd
md"""
### 3.3)

Supongamos que el precio de una casa depende principalmente de dos variables: su edad (``x_1``) y el área de su superficie cubierta (``x_2``) en metros cuadrados. Se propone el modelo lineal ``y=\beta_0+\beta_1x_1+\beta_2x_2``. Se observaron cinco casas seleccionadas aleatoriamente del mercado y se obtuvieron los siguientes datos:

| Precio en miles (``y``) | Edad en años (``x_1``) | Superficie cubierta en m² ×100 (``x_2``) |
|---|---|---|
| 50 | 1 | 1 |
| 40 | 5 | 1 |
| 52 | 5 | 2 |
| 47 | 10 | 2 |
| 65 | 20 | 3 |

Encontrar los valores de ``\beta_0``, ``\beta_1`` y ``\beta_2`` que minimizan el error por cuadrados mínimos. ¿Cuál es el error? Estimar el costo de una vivienda de 25 años y 250 m² de superficie cubierta.
"""

# ╔═╡ 5fd755e7-7a79-42f1-a1d2-0bd3c3bd8826
begin
	precio_miles = [50, 40, 52, 47, 65]
	edad_anios = [1 5 5 10 20]'
	superficie = [1 1 2 2 3]'

	C = [ones(5) edad_anios superficie]
	S = C \ precio_miles

	@info "coeficientes:" S
end

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000070
"""
    trisup_libro(T, c)

Resuelve el sistema ``Tx=c`` mediante sustitución hacia atrás (retrosustitución),
suponiendo ``T`` triangular superior cuadrada con todos los elementos de la
diagonal no nulos. Es la función `trisup` a la que hace referencia el apunte
(ver Sección 2.3).

```julia
function trisup_libro(T, c)
    n = length(c)
    x = copy(c)
    for i in n:-1:1
        x[i] = (x[i] - T[i, i+1:n]' * x[i+1:n]) / T[i, i]
    end
    return x
end
```
"""
function trisup_libro(T::AbstractMatrix{S}, c::AbstractVector{S}) where S
    n = length(c)
    x = copy(c)
    for i in n:-1:1
        x[i] = (x[i] - T[i, i+1:n]' * x[i+1:n]) / T[i, i]
    end
    return x
end

# ╔═╡ c3d4e5f6-7777-4a1b-8c2d-000000000071
"""
    qrsolve_libro(U, R, b)

Utiliza la descomposición QR de ``A`` dada por `qrhouse_libro` (en `U` y `R`)
para resolver el problema de cuadrados mínimos ``min|Ax-b|``.
Devuelve la solución `x` y el error `e = |Ax-b|`.

```julia
function qrsolve_libro(U, R, b)
    n, m = size(U)
    c = copy(b)

    for i in 1:m
        u = U[1:n-i+1, i]
        c[i:n] .-= u * (u' * c[i:n])
    end

    x = trisup_libro(R[1:m, :], c[1:m])
    e = norm(c[m+1:n])

    return x, e
end
```
"""
function qrsolve_libro(U::AbstractMatrix{T}, R::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n, m = size(U)
    c = copy(b)

    for i in 1:m
        u = U[1:n-i+1, i]
        c[i:n] .-= u * (u' * c[i:n])
    end

    x = trisup_libro(R[1:m, :], c[1:m])
    e = norm(c[m+1:n])

    return x, e
end

# ╔═╡ 6da13ed3-84fb-41d4-84ba-e91daee28eb8
curva_regresion(S, x1, x2) = S[1] + S[2] * x1 + S[3] * x2

# ╔═╡ bb79b6e6-ba51-46ad-93c9-196eb2f8ae9d
precio_prediccion = curva_regresion(S, 25, 2.5)

# ╔═╡ 1953d197-4f1b-4c50-89c9-c3fd3937673b
begin
	m = lm(C, precio_miles)
	predict(m, [1 25 2.5])
end

# ╔═╡ 25b3cc26-92b9-4be7-ba46-2dde89f8ecdb
error = sum((curva_regresion.(Ref(S), edad_anios, superficie) - precio_miles).^2)

# ╔═╡ c18d028d-43cd-41dd-858a-56522afadae5
residuals(m)

# ╔═╡ c3d4e5f6-9999-4a1b-8c2d-000000000001
md"""
------------------------------
### 3.4) Transformaciones de Householder

#### 3.4.a)
Programar una función `[x,e] = mincuadqr(A,b)` que resuelva el problema de cuadrados mínimos correspondiente al sistema ``Ax=b``, utilizando las dos funciones `[U,R]=qrhouse(A)` y `[x,e]=qrsolve(U,R,b)` vistas en clase.

#### 3.4.b)
Utilizando la función `mincuadqr(A,b)`, crear la función `a = ajustar(x,y,n)` que reciba como argumentos:
- `x`, `y`: coordenadas de puntos a ajustar con un polinomio.
- `n`: grado del polinomio con el que se ajustarán los puntos.

La salida deben ser los coeficientes del polinomio obtenido, donde `a[1]` es el coeficiente de orden `n`, `a[2]` el de orden `n-1`, y así sucesivamente, hasta `a[n+1]` que es el coeficiente de orden cero.
"""

# ╔═╡ c3d4e5f6-9999-4a1b-8c2d-000000000002
md"""
------------------------------
### 3.5)
Realizar el conteo de operaciones necesarias para resolver el problema de mínimos cuadrados asociado al sistema ``Ax=b`` con ``A\in\mathbb{R}^{n\times k}`` con ``n`` mucho mayor que ``k`` mediante:

#### 3.5.a)
El método de las ecuaciones normales, es decir, resolviendo el sistema ``A^TAx=A^Tb`` mediante la descomposición de Cholesky (Sección 2.7).

#### 3.5.b)
La descomposición QR de ``A`` utilizando el método de Gram-Schmidt.

#### 3.5.c)
La descomposición QR de ``A`` utilizando las transformaciones de Householder (función `[x,e] = mincuadqr(A,b)`).
"""

# ╔═╡ c3d4e5f6-9999-4a1b-8c2d-000000000003
md"""
------------------------------
### 3.6)
Analizar experimentalmente la velocidad y la precisión (o estabilidad) de las funciones `mincuad(A,b)`, `mincuadqr(A,b)` y `\\` (de MATLAB/OCTAVE; en Julia, `A\\b`) al resolver el siguiente problema de mínimos cuadrados: ``Ax=b`` donde ``A\in\mathbb{R}^{n\times4}`` es la matriz `A = [a.^3 a.^2 a ones(n,1)]` para `a = range(1,5,length=n)`, y el vector ``b`` elegido de manera que la solución exacta del problema sea ``x=(1,1,1,1)^T``. Considerar distintos valores grandes de ``n``.
"""

# ╔═╡ c3d4e5f6-9999-4a1b-8c2d-000000000004
md"""
------------------------------
### 3.7)
Realizar un conteo de las operaciones necesarias para hallar la factorización QR de una matriz *upper Hessenberg* de ``n\times n``.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
GLM = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.8.0"
GLM = "~1.9.5"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "51a783cb0bdd0d2c5c1c33cdac8bac0a8ac5686d"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
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

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "PrecompileTools", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "9670d3febc2b6da60a0ae57846ba74670290653f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.8.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.1+2"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "d2facc77c08c1c2bfb1a77c148edd05b3db5410b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.130"

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
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bad39456d9f0166184fce2248783dd9862645c1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "LogExpFunctions", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsModels"]
git-tree-sha1 = "c963639ae5b9aab54f543bdc7504f42f59880bec"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.5"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "31bb6c92405c084617facc1d7ed9eb6c402d061e"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.30"

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

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c7345ab1a7ca4dc8a02c9f6510da0d9857bbe513"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.7.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "df7566479bd64f20bd16b09960145e70160ffb3b"
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

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
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
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

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

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "3de8f5e6e90ebfa8d6d1f86997d6cdcd6a912ff3"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.7"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

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

[[deps.Profile]]
deps = ["StyledStrings"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "13a9e0164267bb9ebc55c1c5d72fceea8f09555b"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.7"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "429071b23f4c9a13fb6582f807cc2ef454082408"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "e4d7a1a0edc20af42689ea6f4f3587a2175d50ee"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.12"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "0db41c4e0d9f3fa195395a6401a8290752c9cd3d"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "2d0fc55c61321ba245c47be599570d11bac50303"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.5"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
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
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"
"""

# ╔═╡ Cell order:
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000001
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000002
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000003
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000004
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000005
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000006
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000007
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000008
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000009
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000000a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000000b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000000c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000000d
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000000e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000000f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000010
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000011
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000012
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000013
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000014
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000015
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000016
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000017
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000018
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000019
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000001a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000001b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000001c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000001d
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000001e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000001f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000020
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000021
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000022
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000023
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000024
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000025
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000026
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000027
# ╟─058d902c-5c99-42ab-bde6-f6d569b22794
# ╟─9927ae15-19e3-45b6-aa0f-ac00b6f866ae
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000028
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000029
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000002a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000002b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000002c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000002d
# ╟─c3d4e5f6-9999-4a1b-8c2d-000000000005
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000002e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000002f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000030
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000031
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000032
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000033
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000034
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000035
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000036
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000037
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000038
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000039
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000003a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000003b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000003c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000003d
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000003e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000003f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000040
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000041
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000042
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000043
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000044
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000045
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000046
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000047
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000048
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000049
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000004a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000004b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000004c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000004d
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000004e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000004f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000050
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000051
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000052
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000053
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000054
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000055
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000056
# ╟─007f43a6-15ec-45fa-be4d-3a25261abf3f
# ╟─95aa7cc0-4534-46d0-91c3-33a721d94435
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000008a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000008b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000008c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000008d
# ╟─a5c99dc0-b35a-4fcf-9547-32a7c7fb40cc
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000057
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000058
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000059
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000005a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000005b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000005c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000005d
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000005e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000005f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000060
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000061
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000062
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000063
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000064
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000065
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000066
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000067
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000068
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000069
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000006a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000006b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000006c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000006d
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000006e
# ╟─d37ade2c-3528-4aa1-91bb-3294a9a4c849
# ╟─22ea694a-3fc1-4752-bac3-77b3b420191d
# ╟─b91efe6d-df9d-43c1-a33e-b944eb1436e1
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000006f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000070
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000071
# ╟─4a096e7e-9f35-4df1-a85b-2192f61fc965
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000072
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000073
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000074
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000075
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000076
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000077
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000078
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000079
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000007a
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000007b
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000007c
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000007d
# ╟─33f6e3d7-50a6-4ad2-a272-e74cbbd59ed8
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000007e
# ╟─c3d4e5f6-7777-4a1b-8c2d-00000000007f
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000080
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000081
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000082
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000083
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000084
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000085
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000086
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000087
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000088
# ╟─c3d4e5f6-7777-4a1b-8c2d-000000000089
# ╟─617150c9-0cb5-4dc8-9b96-e780b5680e62
# ╟─c32c7715-5bd6-4682-8f14-ed52a5e2b9dd
# ╟─db376a73-bfa2-40d0-a282-074bb4ba8bae
# ╟─bd039dc5-e6ce-4a2b-be77-4481ad11a2b5
# ╟─59b39939-05d9-488b-a5a6-2df6581d05d8
# ╟─a3126c2a-b66a-453a-b208-5bb2988a6d8d
# ╟─fe793033-46cd-4c77-87e0-69cd82afdaa1
# ╟─84bfa63d-cfe8-4c51-8017-4dee5f2cc906
# ╟─f3d59908-9bbe-475b-8c49-38aab198368d
# ╟─04f5f8f5-4f44-45d5-abee-7adb1c1e4ae1
# ╟─d38c1824-93bc-4ebf-a700-0850d9335df0
# ╟─c5d58661-106c-468d-af76-0e35194f35aa
# ╟─90c6ba41-71d5-4451-8709-e5d2a3453454
# ╟─41b40cc5-9d7a-4c78-be4b-666367339dbd
# ╟─5fd755e7-7a79-42f1-a1d2-0bd3c3bd8826
# ╟─6da13ed3-84fb-41d4-84ba-e91daee28eb8
# ╟─bb79b6e6-ba51-46ad-93c9-196eb2f8ae9d
# ╟─1953d197-4f1b-4c50-89c9-c3fd3937673b
# ╟─25b3cc26-92b9-4be7-ba46-2dde89f8ecdb
# ╟─c18d028d-43cd-41dd-858a-56522afadae5
# ╟─c3d4e5f6-9999-4a1b-8c2d-000000000001
# ╟─c3d4e5f6-9999-4a1b-8c2d-000000000002
# ╟─c3d4e5f6-9999-4a1b-8c2d-000000000003
# ╟─c3d4e5f6-9999-4a1b-8c2d-000000000004
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
