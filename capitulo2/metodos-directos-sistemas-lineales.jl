### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 43dcd902-9b72-11f1-b709-bf3a6b7c9b3e
using LinearAlgebra, BenchmarkTools, PlutoUI, SIMD

# ╔═╡ 30461d29-5bf1-48d6-b7bd-6ec383304b64
TableOfContents()

# ╔═╡ 387e2688-5e2d-4b9b-991a-f99904f71c8c
md"# Capítulo 2: Métodos directos para sistemas lineales"

# ╔═╡ 84ddd9c4-19a8-4dda-9464-764722d8e553
md"## Breve repaso de Álgebra Lineal"

# ╔═╡ 8f839c88-9f4d-4e40-9349-8e7d3fb9b68e
md"Durante este curso, los vectores serán columna. De este modo, si ``x, y\in\mathbb{R}^n``, entonces"

# ╔═╡ a0600810-db07-4e5b-8fc3-c56d0711dc13
md"""
- Producto escalar: ``x\cdot y = x^Ty\in \mathbb{R}``.
- Producto tensorial: ``x\otimes y = xy^T\in \mathcal{M}_{n\times n}``, con
```math
\text{fila}_{i}(xy^T)=x_iy^T,\qquad \text{col}_j(xy^T)=xy_j = y_j x.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000001
md"""
Si ``A\in\mathbb{R}^{n\times m}`` y ``B\in\mathbb{R}^{m\times p}`` entonces ``AB\in\mathbb{R}^{n\times p}`` y

- ``(AB)_{ij} = \text{fila}_i(A)\,\text{col}_j(B)``;
- ``\text{fila}_i(AB) = \text{fila}_i(A)B``;
- ``\text{col}_j(AB) = A\,\text{col}_j(B)``.

Si ``A\in\mathbb{R}^{n\times m}`` y ``x\in\mathbb{R}^m``, entonces ``Ax\in\mathbb{R}^n`` y ``Ax`` es una combinación lineal de las columnas de ``A``,
```math
Ax = x_1\,\text{col}_1(A) + x_2\,\text{col}_2(A) + \cdots + x_m\,\text{col}_m(A),
```
es decir ``Ax\in\text{col}(A)``, el *espacio columna* de ``A``.

Si ``A\in\mathbb{R}^{n\times m}`` y ``x\in\mathbb{R}^n``, entonces ``x^TA\in\mathbb{R}^{1\times m}``, es un vector fila de ``m`` componentes, y es además una combinación lineal de las filas de ``A``:
```math
x^TA = x_1\,\text{fila}_1(A) + x_2\,\text{fila}_2(A) + \cdots + x_n\,\text{fila}_n(A).
```

Si ``A\in\mathbb{R}^{n\times m}`` y ``e_i`` es el vector de ``m`` componentes que cumple ``(e_i)_j=\delta_{ij}``, entonces
```math
Ae_i = \text{col}_i(A).
```

Si ``A\in\mathbb{R}^{n\times m}`` y ``e_i`` es el vector de ``n`` componentes que cumple ``(e_i)_j=\delta_{ij}``, entonces
```math
e_i^TA = \text{fila}_i(A).
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000002
md"""
### Sistemas lineales de ``n\times n``

Si
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000000b
md"""
```math
A = \begin{bmatrix}a_{11}&a_{12}&\dots&a_{1n}\\a_{21}&a_{22}&\dots&a_{2n}\\\vdots&\vdots&\ddots&\vdots\\a_{n1}&a_{n2}&\dots&a_{nn}\end{bmatrix} \qquad\text{y}\qquad b = \begin{bmatrix}b_1\\b_2\\\vdots\\b_n\end{bmatrix}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000000c
md"""
entonces es lo mismo escribir ``Ax=b`` que escribir
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000000d
md"""
```math
\begin{aligned}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n &= b_1\\
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n &= b_2\\
&\ \ \vdots\\
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nn}x_n &= b_n.
\end{aligned}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000000e
md"""
En Julia se resuelve con `x = A\\b`.

Sabemos de Álgebra Lineal que las siguientes afirmaciones son equivalentes para una matriz dada ``A\in\mathbb{R}^{n\times n}``:

1. Para todo ``b\in\mathbb{R}^n``, existe ``x\in\mathbb{R}^n`` tal que ``Ax=b``.
2. Si ``b\in\mathbb{R}^n`` y existe ``x\in\mathbb{R}^n`` tal que ``Ax=b``, entonces tal ``x`` es único.
3. Si ``Ax=0`` entonces ``x=0``.
4. Las columnas de ``A`` son linealmente independientes.
5. Las filas de ``A`` son linealmente independientes.
6. Existe ``A^{-1}\in\mathbb{R}^{n\times n}`` tal que ``AA^{-1}=A^{-1}A=I``.
7. ``\det(A)\neq 0``.

**Observación.** Si ``A\in\mathbb{R}^{30\times 30}`` y ``\det(A)=1``, entonces
```math
\det(0.1\,A) = (0.1)^{30}\det(A) = 10^{-30}\det(A) = 10^{-30}.
```
Este ejemplo muestra que el determinante pequeño no indica si una matriz es *casi singular*. El *número de condición* da información de este tipo (lo veremos al estudiar métodos iterativos).
"""

# ╔═╡ cec62a0c-c997-4196-b067-298c19df73d0
md"""
## Matrices triangulares inferiores
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000003
md"""
**Definición.** Se dice que una matriz ``L\in\mathbb{R}^{n\times n}`` es *triangular inferior* si ``\ell_{ij}=0`` siempre que ``i<j``. Se dice que ``L`` es *triangular inferior unitaria* si además es triangular inferior y ``\ell_{ii}=1``, ``i=1,2,\dots,n``. Esquemáticamente, una matriz triangular inferior luce de la siguiente forma
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000000f
md"""
```math
L = \begin{bmatrix}
\ell_{11}&0&0&\dots&0\\
\ell_{21}&\ell_{22}&0&\dots&0\\
\ell_{31}&\ell_{32}&\ell_{33}&\dots&0\\
\vdots&\vdots&\vdots&\ddots&\vdots\\
\ell_{n1}&\ell_{n2}&\ell_{n3}&\dots&\ell_{nn}
\end{bmatrix}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000010
md"""
y en la triangular inferior unitaria todas las entradas de la diagonal son unos.

Si ``L`` es triangular inferior y quiero resolver el sistema ``Lx=b``, observamos que ``\det(L)=\ell_{11}\ell_{22}\dots\ell_{nn}\neq0`` si y solo si ``\ell_{ii}\neq0`` para todo ``i=1,2,\dots,n``. Por lo tanto, el sistema tendrá solución única si todos los elementos de la diagonal son no nulos. Para resolverlo observamos cómo luce:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000011
md"""
```math
\begin{aligned}
\ell_{11}x_1 &&&&&&&=b_1\\
\ell_{21}x_1&+&\ell_{22}x_2&&&&&=b_2\\
\ell_{31}x_1&+&\ell_{32}x_2&+&\ell_{33}x_3&&&=b_3\\
&&&&&\vdots&&\\
\ell_{n1}x_1&+&\ell_{n2}x_2&+&\cdots&+&\ell_{nn}x_n&=b_n.
\end{aligned}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000012
md"""
Claramente se ve que la solución de este sistema se puede obtener de la siguiente forma, respetando el orden de las operaciones:
```math
\begin{aligned}
x_1 &= b_1/\ell_{11}\\
x_2 &= (b_2-\ell_{21}x_1)/\ell_{22}\\
x_3 &= (b_3-\ell_{31}x_1-\ell_{32}x_2)/\ell_{33}\\
&\ \ \vdots\\
x_n &= (b_n-\ell_{n1}x_1-\ell_{n2}x_2-\cdots-\ell_{n(n-1)}x_{n-1})/\ell_{nn}.
\end{aligned}
```
Este procedimiento se conoce como *sustitución hacia adelante* (*forward substitution*), y es lo que implementan las funciones que siguen.
"""

# ╔═╡ 3eeb588f-a98b-4206-9fb0-e9aab80a58a0
"""
	triinf1_libro(L, b)

Resuelve el sistema ``Lx = b``
suponiendo L triangular inferior cuadrada con todos
los elementos de la diagonal no nulos.

Primera implementacion del libro, pero en Julia (*sin optimizaciones*).
```julia
function triinf1_libro(L, b)
	n = length(b)
	x = zeros(n)

	for i in 1:n
		x[i] = b[i]

		for j in 1:i-1
			x[i] = x[i] - L[i, j] * x[j]
		end

		x[i] = x[i] / L[i, i]
	end

	return x
end
```
"""
function triinf1_libro(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
	n = length(b)
	x = zeros(T, n)

	for i in 1:n
		x[i] = b[i]

		for j in 1:i-1
			x[i] = x[i] - L[i, j] * x[j]
		end

		x[i] = x[i] / L[i, i]
	end

	return x
end

# ╔═╡ f26a9cf4-b22d-488c-b958-7ea388490af7
"""
	triinf_libro(L, b)

Resuelve el sistema ``Lx = b`` suponiendo ``L`` triangular inferior cuadrada con todos los elementos de la diagonal no nulos.

Segunda implementación del libro (*sin optimizaciones*)
```julia
function triinf_libro(L, b)
	n = length(b)
	x = copy(b)

	for i in 1:n
		for j in 1:i-1
			x[i] = x[i] - L[i, j] * x[j]
		end

		x[i] = x[i] / L[i, i]
	end

	return x
end
```
"""
function triinf_libro(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
	n = length(b)
	x = copy(b)

	for i in 1:n
		for j in 1:i-1
			x[i] = x[i] - L[i, j] * x[j]
		end

		x[i] = x[i] / L[i, i]
	end

	return x
end

# ╔═╡ b1b0d35f-ee7f-4e22-9bb8-350d1d00363f
"""
	triinf_filas_libro(L, b)

Resuelve el sistema ``Lx = b`` suponiendo ``L`` triangular inferior cuadrada con todos los elementos de la diagonal no nulos.

Se accede la matriz ``L`` por filas.
Tercera implementación del libro (*sin optimizaciones*)
```julia
function triinf_filas_libro(L, b)
	n = length(b)
	x = copy(b)

	for i in 1:n
		x[i] = x[i] - dot(L[i,1:i-1], x[1:i-1])
		x[i] = x[i] / L[i, i]
	end

	return x
end
```
"""
function triinf_filas_libro(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
	n = length(b)
	x = copy(b)

	for i in 1:n
		x[i] = x[i] - dot(L[i,1:i-1], x[1:i-1])
		x[i] = x[i] / L[i, i]
	end

	return x
end

# ╔═╡ 5d6fa8d5-d0f8-4cfb-b0ab-5b351ab29081
"""
	triinf_filas_opt(L, b)

Resuelve el sistema ``Lx = b``
suponiendo ``L`` triangular inferior cuadrada con todos los elementos de la diagonal no nulos.

Se accede a la matriz ``L`` por filas. Implementación optimizada con macros de Julia.
```julia
function triinf_filas_opt(L, b)
	n = length(b)
	x = copy(b)

	@inbounds @views for i in 1:n
		x[i] -= dot(L[i,1:i-1], x[1:i-1])
		x[i] /= L[i, i]
	end

	return x
end
```
"""
function triinf_filas_opt(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
	n = length(b)
	x = copy(b)

	@inbounds @views for i in 1:n
		x[i] -= dot(L[i,1:i-1], x[1:i-1])
		x[i] /= L[i, i]
	end

	return x
end

# ╔═╡ 885e65ff-16ff-4dd3-a759-2533a6a7a635
"""
	triinf_cols_libro(L, b)

Resuelve el sistema ``Lx = b``
suponiendo L triangular inferior cuadrada con todos
los elementos de la diagonal no nulos.

Se accede la matriz ``L`` por columnas.

Primera implementacion del libro, pero en Julia (*sin optimizaciones*).

```julia
function triinf_cols_libro(L, b)
    n = length(b)
    x = copy(b)

    x[1] = x[1] / L[1,1]

    for i in 2:n
        x[i:n] = x[i:n] - x[i-1] * L[i:n, i-1]
        x[i] = x[i] / L[i,i]
    end

    return x
end
```
"""
function triinf_cols_libro(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = length(b)
    x = copy(b)

    x[1] = x[1] / L[1,1]

    for i in 2:n
        x[i:n] = x[i:n] - x[i-1] * L[i:n, i-1]
        x[i] = x[i] / L[i,i]
    end

    return x
end

# ╔═╡ 83c191fc-a03f-4d35-b7c6-fba72ee03aed
"""
	triinf_cols_opt(L, b)

Resuelve el sistema ``Lx = b``
suponiendo L triangular inferior cuadrada con todos
los elementos de la diagonal no nulos.
```julia
function triinf_cols_opt(L, b)
    n = length(b)
    x = copy(b)

    @inbounds begin
        x[1] = x[1] / L[1,1]

        for i in 2:n
            xi = x[i-1]

            @simd for k in i:n
                x[k] = x[k] - xi * L[k,i-1]
            end

            x[i] = x[i] / L[i,i]
        end
    end

    return x
end
```
"""
function triinf_cols_opt(L::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = length(b)
    x = copy(b)

    @inbounds begin
        x[1] = x[1] / L[1,1]

        for i in 2:n
            xi = x[i-1]

            @simd for k in i:n
                x[k] = x[k] - xi * L[k,i-1]
            end

            x[i] = x[i] / L[i,i]
        end
    end

    return x
end

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000004
md"""
**Conteo de operaciones.**

Para conocer la *complejidad* del algoritmo, y entender acerca del *costo computacional*, hacemos un conteo de operaciones. Contemos, por ejemplo, las de `triinf_cols`:

| Momento | multiplicaciones | sumas | divisiones |
|---|---|---|---|
| antes del for | | | 1 |
| primera iteración | ``n-1`` | ``n-1`` | 1 |
| segunda iteración | ``n-2`` | ``n-2`` | 1 |
| tercera iteración | ``n-3`` | ``n-3`` | 1 |
| ``\vdots`` | ``\vdots`` | ``\vdots`` | ``\vdots`` |
| última iteración | 1 | 1 | 1 |
| **total** | ``(n-1)n/2`` | ``(n-1)n/2`` | ``n`` |

Por lo tanto el total de operaciones es ``(n-1)n + n = n^2``. Las demás variantes (`triinf`, `triinf_filas`, etc.) realizan exactamente el mismo número de operaciones; lo que cambia entre ellas es *cómo* se acceden los datos en memoria.

**Observación.** Vale la pena hacer algunos comentarios a esta altura:

- Cualquiera de las versiones que implementamos acá es, en general, más lenta que escribir directamente `L\\b` en Julia: `\\` despacha a rutinas de LAPACK/BLAS, muy optimizadas (multi-hilo, con bloqueo por caché, etc.).
- Las versiones *vectorizadas* (por ejemplo `triinf_cols`, que opera con rebanadas `L[i:n, i-1]`) suelen ser más rápidas que las que hacen todo con bucles escalares explícitos. En Julia esta diferencia es mucho menor que en MATLAB/Octave, ya que Julia compila cada función la primera vez que se ejecuta (*JIT*): un bucle bien tipado (como en `triinf_filas_opt` o `triinf_cols_opt`, con `@inbounds`/`@simd`/`@views`) puede llegar a competir con las versiones vectorizadas.
- Para lenguajes de programación como Fortran, C, C++, Julia, etc., existen rutinas optimizadas para las operaciones siguientes:

  - **saxpy**: ``z\leftarrow \alpha x + y``, con ``\alpha`` escalar y ``x,y`` vectores.
  - **dot**: ``a\leftarrow x\cdot y``, con ``x,y`` vectores.

  Para vectores largos, estas rutinas reparten el trabajo en varios procesadores (paralelizan). En Julia, `LinearAlgebra.BLAS` expone estas rutinas directamente (`BLAS.axpy!`, `BLAS.dot`, etc.), y son las que terminan usándose "por detrás" en operaciones como `L\\b`.
"""

# ╔═╡ 62816356-cfab-45ff-88d9-ab384e26473f
md"""
## Matrices triangulares superiores
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000005
md"""
**Definición.** Se dice que una matriz ``U\in\mathbb{R}^{n\times n}`` es *triangular superior* si ``u_{ij}=0`` siempre que ``i>j``. Esquemáticamente, una matriz triangular superior luce de la siguiente forma
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000013
md"""
```math
U = \begin{bmatrix}
u_{11}&u_{12}&u_{13}&u_{14}&\dots&u_{1n}\\
0&u_{22}&u_{23}&u_{24}&\dots&u_{2n}\\
0&0&u_{33}&u_{34}&\dots&u_{3n}\\
\vdots&\vdots&0&\ddots&\ddots&\vdots\\
\vdots&\vdots&\vdots&&u_{(n-1)(n-1)}&u_{(n-1)n}\\
0&0&0&\dots&0&u_{nn}
\end{bmatrix}.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000014
md"""
La resolución de sistemas ``Ux=b`` con ``U`` triangular superior es análoga a la de sistemas triangulares inferiores, haciendo lo que se conoce como *sustitución hacia atrás* (*backward substitution*): se despeja primero ``x_n``, y se va sustituyendo hacia arriba. También se puede implementar accediendo a ``U`` por filas o por columnas, como hicimos para el caso triangular inferior.
"""

# ╔═╡ f7b1ae60-b108-4860-abfa-8ffc6438a5d4
"""
	trisup_cols_opt(U, b)

Resuelve el sistema Ux = b
suponiendo U triangular superior cuadrada con todos
los elementos de la diagonal no nulos.

```julia
function trisup_cols_opt(U, b)
    n = length(b)
    x = copy(b)

    @inbounds begin
        x[n] = x[n] / U[n,n]

        for i in n-1:-1:1
            xi = x[i+1]

            @simd for k in 1:i
                x[k] = x[k] - xi * U[k,i+1]
            end

            x[i] = x[i] / U[i,i]
        end
    end
    return x
end
```
"""
function trisup_cols_opt(U::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = length(b)
    x = copy(b)
    @inbounds begin
        x[n] = x[n] / U[n,n]
        for i in n-1:-1:1
            xi = x[i+1]
            @simd for k in 1:i
                x[k] = x[k] - xi * U[k,i+1]
            end
            x[i] = x[i] / U[i,i]
        end
    end
    return x
end

# ╔═╡ 5d7e0747-f9c9-4794-98d0-dbbbb0d05e0a
md"""
## Matrices llenas. Método de Gauss
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000006
md"""
Consideremos el problema de resolver ``Ax=b`` con ``A\in\mathbb{R}^{n\times n}``, ``b\in\mathbb{R}^n`` dados, y ``A`` invertible. El método de eliminación de Gauss consiste en realizar operaciones elementales sobre la matriz ``A`` y el lado derecho ``b`` de manera que el sistema se transforme en un sistema triangular superior.

**Primer paso** ``[A=A^{(0)}\to A^{(1)},\ b=b^{(0)}\to b^{(1)}]``. Restar de las filas 2 a ``n`` un múltiplo de la fila 1, para obtener ceros debajo de ``a_{11}``. Esto se logra de la siguiente manera:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000029
md"""
```math
\text{fila}_1(A^{(1)}) = \text{fila}_1(A^{(0)}), \qquad \text{fila}_k(A^{(1)}) = \text{fila}_k(A^{(0)}) - \underbrace{\frac{a_{k1}^{(0)}}{a_{11}^{(0)}}}_{m_{k1}}\text{fila}_1(A^{(0)}), \quad k=2,3,\dots,n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000002a
md"""
```math
b_1^{(1)} = b_1^{(0)}, \qquad b_k^{(1)} = b_k^{(0)} - m_{k1}\,b_1^{(0)}, \quad k=2,3,\dots,n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000002b
md"""
**Segundo paso** ``[A^{(1)}\to A^{(2)},\ b^{(1)}\to b^{(2)}]``. Restar de las filas 3 a ``n`` un múltiplo de la fila 2, para obtener ceros debajo de ``a_{22}``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000002c
md"""
```math
\text{fila}_1(A^{(2)}) = \text{fila}_1(A^{(1)}), \quad \text{fila}_2(A^{(2)}) = \text{fila}_2(A^{(1)}), \quad \text{fila}_k(A^{(2)}) = \text{fila}_k(A^{(1)}) - \underbrace{\frac{a_{k2}^{(1)}}{a_{22}^{(1)}}}_{m_{k2}}\text{fila}_2(A^{(1)}), \quad k=3,\dots,n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000002d
md"""
```math
b_1^{(2)} = b_1^{(1)}, \quad b_2^{(2)} = b_2^{(1)}, \quad b_k^{(2)} = b_k^{(1)} - m_{k2}\,b_2^{(1)}, \quad k=3,\dots,n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000002e
md"""
**Paso ``j`` (``j\leq n``)** ``[A^{(j-1)}\to A^{(j)},\ b^{(j-1)}\to b^{(j)}]``. Restar de las filas ``j+1`` a ``n`` un múltiplo de la fila ``j``, para obtener ceros debajo de ``a_{jj}``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000002f
md"""
```math
\text{fila}_k(A^{(j)}) = \text{fila}_k(A^{(j-1)}), \quad k=1,\dots,j, \qquad \text{fila}_k(A^{(j)}) = \text{fila}_k(A^{(j-1)}) - \underbrace{\frac{a_{kj}^{(j-1)}}{a_{jj}^{(j-1)}}}_{m_{kj}}\text{fila}_j(A^{(j-1)}), \quad k=j+1,\dots,n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000030
md"""
```math
b_k^{(j)} = b_k^{(j-1)}, \quad k=1,\dots,j, \qquad b_k^{(j)} = b_k^{(j-1)} - m_{kj}\,b_j^{(j-1)}, \quad k=j+1,\dots,n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000031
md"""
Al finalizar el paso ``n-1`` obtenemos el sistema
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000032
md"""
```math
A^{(n-1)}x = b^{(n-1)},
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000033
md"""
equivalente al original, con la ventaja de que ``A^{(n-1)}=U`` es triangular superior, y puede resolverse por sustitución hacia atrás.

**Conteo de operaciones.**

Observemos que donde sabemos que el resultado será cero, no hacemos las operaciones, sino que directamente colocamos un cero. Denotaremos con ``d``, ``m``, ``r`` las divisiones, multiplicaciones y restas, respectivamente.

| Paso | operaciones sobre ``A`` | operaciones sobre ``b`` | cuántas veces |
|---|---|---|---|
| ``1`` | ``[1d + (n-1)m + (n-1)r]`` | ``(1m+1r)`` | ``(n-1)`` |
| ``j,\ 2\leq j\leq n-1`` | ``[1d + (n-j)m + (n-j)r]`` | ``(1m+1r)`` | ``(n-j)`` |

Total de operaciones para la *triangulación* sobre la matriz ``A``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000034
md"""
```math
\sum_{j=1}^{n-1}(1+2(n-j))(n-j) = \sum_{k=1}^{n-1}(1+2k)k = \sum_{k=1}^{n-1}k + 2\sum_{k=1}^{n-1}k^2 = \frac{(n-1)n}{2} + 2\frac{(n-1)n(2n-1)}{6} = \frac{2}{3}n^3-\frac{1}{2}n^2-\frac{1}{6}n \cong \frac{2}{3}n^3.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000035
md"""
Total de operaciones al reducir a forma triangular, realizadas sobre el lado derecho ``b``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000036
md"""
```math
\sum_{j=1}^{n-1}2(n-j) = 2\sum_{j=1}^{n-1}j = (n-1)n = n^2-n.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000037
md"""
Hasta aquí, el total de operaciones para reducir a forma triangular un sistema con el método de Gauss es
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000038
md"""
```math
\frac{2}{3}n^3-\frac{1}{2}n^2-\frac{1}{6}n + n^2-n = \frac{2}{3}n^3+\frac{1}{2}n^2-\frac{7}{6}n \cong \frac{2}{3}n^3.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000039
md"""
Si sumamos el número de operaciones para la eliminación hacia atrás (``n^2``, contadas al analizar `trisup`/`triinf`), obtenemos que para resolver un sistema con el método de Gauss hacen falta
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000003a
md"""
```math
\frac{2}{3}n^3+\frac{1}{2}n^2-\frac{7}{6}n + n^2 = \frac{2}{3}n^3+\frac{3}{2}n^2-\frac{7}{6}n \cong \frac{2}{3}n^3 \quad\text{operaciones.}
```
"""

# ╔═╡ a9f741f4-9e65-4fa6-b286-3a621aab6633
"""
    metodo_de_gauss1(M, b)

Implementación directa del método de Gauss.

```julia
function metodo_de_gauss1(M, b)
    m = [M b]
    n = size(M, 1)

    for c in 1:n-1
        idx = argmax(abs.(m[c+1:n, c])) + c

        m[c, :], m[idx, :] = m[idx, :], m[c, :]

        for r in (c+1):n
            m[r, :] = m[r, :] - m[r, c] / m[c, c] .* m[c, :]
        end
    end

    return m[:, 1:end-1], m[:, end]
end
```
"""
function metodo_de_gauss1(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
    m = [M b]
    n = size(M, 1)

    for c in 1:n-1
        idx = argmax(abs.(m[c+1:n, c])) + c

        m[c, :], m[idx, :] = m[idx, :], m[c, :]

        for r in (c+1):n
            m[r, :] = m[r, :] - m[r, c] / m[c, c] .* m[c, :]
        end
    end

    return m[:, 1:end-1], m[:, end]
end

# ╔═╡ 2d6ac9d5-8423-4e74-9abd-3d192e1668fa
"""
    function metodo_de_gauss_opt(M, b)

Implementación optimizada del método de Gauss. (*column mayor*)

```julia
function metodo_de_gauss_opt(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = size(M, 1)
    m = [M b]

    #para no pedir memoria cada vez (aloja los m)
    factores = Vector{T}(undef, n)

    @inbounds for c in 1:n-1
        # buscar pos de pivote sin alojar nada
        maxval = abs(m[c, c])
        idx = c
        for r in c+1:n
            v = abs(m[r, c])
            if v > maxval
                maxval = v
                idx = r
            end
        end

        # intercambiar sin alojar (solo si idx ≠ c)
        if idx != c
            @simd for j in 1:n+1
                m[c, j], m[idx, j] = m[idx, j], m[c, j]
            end
        end

        invpiv = one(T) / m[c, c]
        for r in c+1:n
            factores[r] = m[r, c] * invpiv
        end

        # loop invertido (col fuera, fila dnetor)
        for j in c:n+1
            mcj = m[c, j]
            @simd for r in c+1:n
                m[r, j] -= factores[r] * mcj
            end
        end
    end

    return m[:, 1:end-1], m[:, end]
end
```
"""
function metodo_de_gauss_opt(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = size(M, 1)
    m = [M b]

    factores = Vector{T}(undef, n)

    @inbounds for c in 1:n-1
        maxval = abs(m[c, c])
        idx = c
        for r in c+1:n
            v = abs(m[r, c])
            if v > maxval
                maxval = v
                idx = r
            end
        end

        if idx != c
            @simd for j in 1:n+1
                m[c, j], m[idx, j] = m[idx, j], m[c, j]
            end
        end

        invpiv = one(T) / m[c, c]
        @simd for r in c+1:n
            factores[r] = m[r, c] * invpiv
        end

        for j in c:n+1
            mcj = m[c, j]
            @simd for r in c+1:n
                m[r, j] -= factores[r] * mcj
            end
        end
    end

    return m[:, 1:end-1], m[:, end]
end

# ╔═╡ 51302fe8-2eb6-4e45-aea9-b12d7a434481
md"""
## Factorización LU
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000007
md"""
Recordando que ``\text{fila}_i(GA)=\text{fila}_i(G)A``, observamos que
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000015
md"""
```math
A^{(1)} = \begin{bmatrix}1&0&0&\cdots&0\\-m_{21}&1&0&\cdots&0\\-m_{31}&0&1&\cdots&0\\\vdots&\vdots&\ddots&\ddots&\vdots\\-m_{n1}&0&\cdots&0&1\end{bmatrix}A^{(0)}.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000016
md"""
Más precisamente, el primer paso de la eliminación de Gauss se escribe
```math
\text{si}\quad g_1=\begin{bmatrix}0\\m_{21}\\\vdots\\m_{n1}\end{bmatrix}, \quad\text{y}\quad G^1=I-g_1e_1^T, \quad\text{entonces}\quad A^{(1)}=G^1A^{(0)}=G^1A.
```
Análogamente, el segundo paso se escribe ``A^{(2)}=G^2A^{(1)}``, con ``g_2=(0,0,m_{32},\dots,m_{n2})^T`` y ``G^2=I-g_2e_2^T``. El ``j``-ésimo paso resulta ``A^{(j)}=G^jA^{(j-1)}``, con ``g_j=(0,\dots,0,m_{(j+1)j},\dots,m_{nj})^T`` y ``G^j=I-g_je_j^T``.

Las matrices ``G^j=I-g_je_j^T`` se llaman *matrices de Gauss*. En el Ejercicio 2.6 (más abajo) se demuestra que ``(G^j)^{-1}=I+g_je_j^T``, y que
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000017
md"""
```math
(G^1)^{-1}(G^2)^{-1}\cdots(G^{n-1})^{-1} = I+\sum_{j=1}^{n-1}g_je_j^T = \begin{bmatrix}1&0&0&\cdots&0\\m_{21}&1&0&\cdots&0\\m_{31}&m_{32}&1&\cdots&0\\\vdots&\vdots&\ddots&\ddots&\vdots\\m_{n1}&m_{n2}&m_{n3}&\cdots&1\end{bmatrix} = L.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000018
md"""
**Conclusión:** El método de eliminación de Gauss produce (si ``a_{jj}^{(j-1)}\neq0``, ``j=1,2,\dots,n-1``) una factorización ``A=LU`` con ``U`` triangular superior, ``L`` triangular inferior unitaria, además ``L_{ij}=m_{ij}=a_{ij}^{(j-1)}/a_{jj}^{(j-1)}``, el *multiplicador*, si ``i>j``.

**Observación.**

- Si ``A=LU``, resolver ``Ax=b`` se desglosa en la resolución de dos sistemas triangulares: ``Ax=b \implies L\underbrace{Ux}_{y}=b``. Resolvemos primero ``Ly=b`` (``n^2-n`` operaciones), y luego ``Ux=y`` (``n^2`` operaciones).
- Una vez obtenida la factorización (``\frac{2}{3}n^3-\frac{1}{2}n^2-\frac{1}{6}n`` operaciones), la resolución del sistema lleva ``\frac{2}{3}n^3-\frac{1}{2}n^2-\frac{1}{6}n+2n^2-n = \frac{2}{3}n^3+\frac{3}{2}n^2-\frac{7}{6}n`` operaciones. Esto es exactamente lo mismo que resolver por el método de Gauss, y no debería sorprendernos.
- El vector ``y`` obtenido al resolver el primer sistema triangular inferior coincide con ``b^{(n-1)}``.
- La utilidad de la factorización se ve al resolver muchos sistemas con igual matriz ``A`` y diferentes lados derechos ``b``. Para resolver ``k`` sistemas utilizando el método de Gauss cada vez se necesitan ``k(\frac{2}{3}n^3+\frac{3}{2}n^2-\frac{7}{6}n)\cong k\frac{2}{3}n^3`` operaciones. Si se factoriza primero y luego se resuelven los dos sistemas triangulares ``k`` veces, se necesitan ``(\frac{2}{3}n^3-\frac{1}{2}n^2-\frac{1}{6}n) + k(2n^2-n) \cong \frac{2}{3}n^3+kn^2`` operaciones.
- Las matrices ``L`` y ``U`` se almacenan (en general) en el mismo lugar en que estaba la matriz original ``A``. Los multiplicadores ``m_{ij}`` se almacenan en la posición ``(i,j)``, que es justamente la posición donde se *produce* un cero.
"""

# ╔═╡ 2b09d162-6b0a-45d5-b966-119c5be4b476
"""
    eleu_libro(A)

Devuelve la factorizacion ``LU`` sin pivoteo.

```julia
function eleu_libro(M)
    n = size(M, 1)
    U = copy(M)
    L = Matrix{Float64}(I, n, n)

    for c in 1:n-1
        @views for r in (c+1):n
            m = U[r, c] / U[c, c]
            U[r, :] .-= m .* U[c, :]
            L[r, c] = m
        end
    end
    return L, U
end
```
"""
function eleu_libro(M::AbstractMatrix)
    n = size(M, 1)
    U = float(copy(M))
    L = Matrix{Float64}(I, n, n)
    
    for c in 1:n-1        
        @views for r in (c+1):n
            m = U[r, c] / U[c, c]
            U[r, :] .-= m .* U[c, :]
            L[r, c] = m
        end
    end
    return L, U
end

# ╔═╡ ab7c1a16-c76e-49d9-be60-feb360451a84
md"""
## Problemas con la factorización ``LU``. Pivoteo
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000008
md"""
Podría suceder que ``a_{11}=0`` o que ``a_{jj}^{(j-1)}=0``, en cuyo caso el algoritmo no funciona. Más aún, puede verse que no es posible hallar ``\ell`` y ``u_{11},u_{12},u_{22}`` tales que
```math
\begin{bmatrix}0&1\\1&1\end{bmatrix} = \begin{bmatrix}1&0\\\ell&1\end{bmatrix}\begin{bmatrix}u_{11}&u_{12}\\0&u_{22}\end{bmatrix}.
```
Este simple ejemplo muestra que no siempre es posible hallar una factorización ``LU`` de ``A``, con ``L`` triangular inferior unitaria y ``U`` triangular superior. (Este es, precisamente, el Ejercicio 2.7 de más abajo.)

El siguiente teorema da una condición suficiente para que el algoritmo visto funcione y logre la factorización ``LU``.

**Teorema.** *Sea ``A`` una matriz de ``n\times n``. Si ``A_{1:k,1:k}`` es no singular, es decir, ``\det(A_{1:k,1:k})\neq0``, para todo ``k=1,2,\dots,n-1``, entonces existe descomposición ``LU`` de ``A``. Más precisamente, existe ``L`` triangular inferior unitaria y ``U`` triangular superior tales que ``A=LU``. Además, si ``A`` es invertible, esta descomposición ``LU`` es única.*

**Demostración.** Basta ver que en el proceso de eliminación de Gauss resulta siempre ``a_{jj}^{(j-1)}\neq0`` para ``j=1,2,\dots,n-1``. Demostremos esto por inducción. Si ``j=1``,
```math
a_{1,1}^{(0)} = a_{11} = \det(A_{1:1,1:1}) \neq 0.
```
Supongamos que ``a_{jj}^{(j-1)}\neq0``, ``j=1,2,\dots,k``, y recordemos que ``A^{(k)}=G^kG^{k-1}\dots G^2G^1A``. Como las matrices ``G^j`` son matrices de Gauss (y ``\det(G^j)=1``), resulta que
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000019
md"""
```math
A^{(k)}_{1:k+1,1:k+1} = G^k_{1:k+1,1:k+1}\dots G^2_{1:k+1,1:k+1}G^1_{1:k+1,1:k+1}A_{1:k+1,1:k+1} = \begin{bmatrix}a_{11}^{(0)}&a_{12}^{(0)}&\cdots&a_{1(k+1)}^{(0)}\\0&a_{22}^{(1)}&\cdots&a_{2(k+1)}^{(1)}\\0&0&\ddots&\vdots\\0&0&\cdots&a_{(k+1)(k+1)}^{(k)}\end{bmatrix}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000001a
md"""
y por lo tanto
```math
\det(A^{(k)}_{1:k+1,1:k+1}) = \det(G^k_{1:k+1,1:k+1})\dots\det(G^1_{1:k+1,1:k+1})\det(A_{1:k+1,1:k+1}) = a_{11}^{(0)}a_{22}^{(1)}\cdots a_{(k+1)(k+1)}^{(k)},
```
como esta cantidad es no nula por hipótesis inductiva, resulta que ``a_{(k+1)(k+1)}^{(k)}\neq0``, que es lo que se quería demostrar.

Supongamos ahora que ``A`` es no singular. Para demostrar la unicidad de la descomposición ``LU`` supongamos que hay matrices ``L_1,L_2,U_1,U_2``, tales que ``L_i`` es triangular inferior unitaria y ``U_i`` es triangular superior, ``i=1,2``, con ``L_1U_1=L_2U_2=A``. Luego ``L_2^{-1}L_1=U_2U_1^{-1}``. Como ``L_2^{-1}L_1`` resulta triangular inferior unitaria, y ``U_2U_1^{-1}`` es triangular superior, tenemos que ``L_2^{-1}L_1=U_2U_1^{-1}=I``, por lo que ``L_1=L_2`` y ``U_1=U_2``.

**Pivoteo.** Para poder resolver sistemas que no cumplen con las hipótesis del teorema anterior se utiliza el *pivoteo*. El pivoteo (parcial) consiste en intercambiar, en el paso ``j``-ésimo de la eliminación de Gauss, la fila ``j`` con la fila ``k``, que contiene al elemento de ``A^{(j-1)}_{j:n,j}`` con máximo valor absoluto, antes de realizar las operaciones por filas. De esta manera, uno se asegura que nunca divide por cero, y además que divide por el elemento de módulo máximo, de los posibles *pivotes*. Esto es beneficioso pues dividir por números pequeños aumenta la propagación de errores de redondeo.

**Observación.** Este pivoteo se conoce como *pivoteo parcial*. El *pivoteo total* consiste en hallar el elemento de ``A^{(j-1)}_{j:n,j:n}`` con módulo máximo e intercambiar filas y columnas. Es más complicado y en la práctica no se observan beneficios significativos con respecto al pivoteo parcial.
"""

# ╔═╡ faa6b8f4-656c-4f7c-a761-e7842c89ad29
"""
    L, U, p = eleupiv_libro(A)

Decomposición ``LU`` con pivoteo.
- ``L`` será triangular inferior unitaria.
- ``U`` será triangular superior
- ``LU = A[p, :]`` permutación ``p`` de las filas de ``A``.
```julia
function eleupiv_libro(A::AbstractMatrix)
    n = size(A, 1)
    L, U, p = zeros(n, n), copy(A), collect(1:n)

    for c in 1:n-1
        idx = argmax(abs.(U[c:n, c])) + c - 1

        p[[idx, c]] = p[[c, idx]]
        L[[idx, c], 1:c-1] = L[[c, idx], 1:c-1]
        U[[idx, c], :] = U[[c, idx], :]

        L[c+1:n, c] = U[c+1:n, c] / U[c, c]
        U[c+1:n, c] .= 0
        U[c+1:n, c+1:n] .-= L[c+1:n, c] * U[c, c+1:n]'
    end

    LowerTriangular(L + I), UpperTriangular(U), p
end
```
"""
function eleupiv_libro(A::AbstractMatrix)
    n = size(A, 1)
    L, U, p = zeros(n, n), copy(A), collect(1:n)

    for c in 1:n-1
        idx = argmax(abs.(U[c:n, c])) + c - 1

        p[[idx, c]] = p[[c, idx]]
        L[[idx, c], 1:c-1] = L[[c, idx], 1:c-1]
        U[[idx, c], :] = U[[c, idx], :]

        L[c+1:n, c] = U[c+1:n, c] / U[c, c]
        U[c+1:n, c] .= 0
        U[c+1:n, c+1:n] .-= L[c+1:n, c] * U[c, c+1:n]'
    end

    LowerTriangular(L + I), UpperTriangular(U), p
end

# ╔═╡ ccebc4a5-b633-46f8-8153-e88cd0dbf717
md"""
### Permutaciones
"""

# ╔═╡ 854bbb92-5282-4f9f-9373-a764aeae135a
md"""
Si ``P`` se obtiene intercambiando la fila ``j`` y ``k`` de una matriz identidad, entonces ``P`` sigue teniendo un uno por cada fila (intercambiar filas no cambia la disposición de los elementos de la fila) (🍎) y un uno por cada columna (dado que cada columna tiene un solo uno y el resto ceros, intercambiar filas solo reordena estos valores) (🍐).



-----------------------------------------------
**Proposición:** ``P^T=P``

**Demostración**

Supongamos que de la matriz identidad ``I`` itercambiamos la fila ``i`` por la fila ``j``.


Para ver que  ``P^T = P`` bastará con probar que

- ``p_{i,k} = p_{k,i}`` para ``1\leq k \leq n`` (🐱) y 
- ``p_{j,k} = p_{k,j}`` para ``1\leq k \leq n`` (🐶)

pues estas dos son las únicas filas que alteramos.

- (🐱) Distingamos entre el caso ``k=j`` y ``k\neq j``. Comencemos con ``k=j``. Dado que el uno de la fila ``i`` está en la componente ``ii`` de ``I``, entonces estará en ``ji`` de ``P`` y, dado que el uno de la fila ``j`` está en la componente ``jj`` de ``I``, estará en ``ij`` de ``P``. Así ``p_{ij} = p_{ji} = 1``. En lo que respecta a ``p_{i,k}`` con ``k\neq j``, estas componentes valen todas cero (🍎), y por (🍐)  todas las compontes de la columna ``i`` (salvo ``p^{ji}`` como ya se vió) también valen  cero. Luego, ``p_{i,k} = p_{k,i}`` para ``1\leq k\leq n``.

- (🐶) La demostración es equivalente a la anterior.

Luego ``P = P^T``.




-----------------------------------------------
**Proposición:** ``PP = I``

**Demostración**

Veamos que ``PP=Y_{y_{ij}}`` solo tiene 1 en ``y_{ii}`` y ceros en ``y_{ij}`` con ``i\neq j``.

- Para el cálculo de ``y_{ii}`` tenemos ``\sum_{k=1}^{n}p_{ik}p_{ki}``. Dado que en el inciso anterior se probó que P es simétrica, ``p_{ik}`` y ``p_{ki}`` deben valer 1 simultaneamente y, por 🍎 y 🍐, lo hacen solo una vez. Entonces ``y_{ii} = \sum_{k=1}^{n}p_{ik}p_{ki} = 0 + 0 + \cdots + 0 + 1 + 0 + \cdots + 0 + 0 = 1``

- Para el cálculo de  ``y_{ij}``, con ``j\neq i`` tenemos ``\sum_{k=1}^{n}p_{ik}p_{kj}``. Supongamos que las componente ``p_{im}`` de la fila ``i`` vale uno. Por 🍎 solo está componente vale uno y el resto de componentes vale cero. Así, 
```math
y_{ij} = \sum_{k=1}^{n}p_{ik}p_{kj} = p_{im}p_{mj}
```

Como ``i\neq j``, ``p_{mj}`` debe valer cero, pues al ser ``P`` simétrica y por 🍐, sabemos que, en la columna ``i``, solo ``p_{mi}`` es distinta de cero. Luego, ``y_{ij}=0``.





-----------------------------------------------
**Proposición:** ``PA`` se obtiene intercambiando las fijas ``j`` y ``k`` de ``A``

**Demostración**



-----------------------------------------------
**Proposición:** ``AP`` se obtiene intercambiando las columnas ``j`` y ``k`` de ``A``

**Demostración**
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000009
md"""
Con las permutaciones ya definidas, los pasos del método de Gauss con pivoteo se pueden escribir de la siguiente manera:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000003b
md"""
```math
A^{(1)} = G^1P^1A^{(0)} = G^1P^1A, \qquad A^{(2)} = G^2P^2A^{(1)} = G^2P^2G^1P^1A, \qquad \dots
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000003c
md"""
```math
U = A^{(n-1)} = G^{n-1}P^{n-1}G^{n-2}P^{n-2}\dots G^2P^2G^1P^1A,
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000003d
md"""
con ``P^i`` una matriz de permutación que se obtiene intercambiando la fila ``i`` con una fila ``k\geq i`` y ``G^i=I-g_ie_i^T``,
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000003e
md"""
```math
g_i = (0,\dots,0,m_{(i+1),i},\dots,m_{ni})^T, \qquad\text{y}\qquad m_{\ell i} = \frac{(P^iA^{(i-1)})_{\ell i}}{(P^iA^{(i-1)})_{ii}}.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000003f
md"""
**Teorema.** *Si ``U=A^{(n-1)}`` de la eliminación de Gauss con pivoteo, entonces*
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000040
md"""
```math
G^{n-1}P^{n-1}G^{n-2}P^{n-2}\dots G^2P^2G^1P^1A = U,
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000041
md"""
*con ``U`` triangular superior, y además*
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000042
md"""
```math
PA=LU
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000043
md"""
*con ``P=P^{n-1}P^{n-2}\dots P^2P^1`` y ``L`` triangular inferior unitaria con ``|\ell_{ij}|\leq1``. Además, la ``k``-ésima columna de ``L``, debajo de la diagonal, es una permutación de los multiplicadores del paso ``k``, más precisamente*
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000044
md"""
```math
L_{k+1:n,k} = (P^{n-1}\dots P^{k+2}P^{k+1}g_k)_{k+1:n}.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000045
md"""
**Demostración.** Por la construcción se ve que ``U`` es triangular superior. Veamos ahora que ``PA=LU`` y que ``L`` tiene las características mencionadas. Como para cualquier permutación ocurre que ``PP=I``, resulta
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000046
md"""
```math
\begin{aligned}
U &= G^{n-1}P^{n-1}G^{n-2}P^{n-2}\dots G^2P^2G^1P^1A\\
&= G^{n-1}(P^{n-1}G^{n-2}P^{n-1})(P^{n-1}P^{n-2}G^{n-3}P^{n-2}P^{n-1})\\
&\quad\times(P^{n-1}P^{n-2}P^{n-3}G^{n-4}P^{n-3}P^{n-2}P^{n-1})\dots\\
&\quad\times P^{n-1}\dots P^3P^2G^1P^2P^3\dots P^{n-1}\\
&\quad\times \underbrace{P^{n-1}P^{n-2}\dots P^3P^2P^1}_{P}A\\
&= \tilde G^{n-1}\tilde G^{n-2}\dots\tilde G^2\tilde G^1 PA,
\end{aligned}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000047
md"""
donde
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000048
md"""
```math
\tilde G^k = P^{n-1}\dots P^{k+1}G^kP^{k+1}\dots P^{n-1} = I - \underbrace{P^{n-1}\dots P^{k+1}g_k}_{\tilde g_k}e_k^T = I-\tilde g_ke_k^T.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000049
md"""
Cada ``\tilde G^k`` es entonces una matriz de Gauss con vector ``\tilde g_k``, y por la Conclusión de la sección de Factorización LU, ``(\tilde G^{n-1})^{-1}\dots(\tilde G^1)^{-1} = I+\sum_{k=1}^{n-1}\tilde g_ke_k^T = L``, triangular inferior unitaria con ``|\ell_{ij}|\leq1``, y ``L_{k+1:n,k}`` es exactamente ``(\tilde g_k)_{k+1:n} = (P^{n-1}\dots P^{k+1}g_k)_{k+1:n}``. Luego ``U=(\tilde G^{n-1}\dots\tilde G^1)PA = L^{-1}PA``, es decir ``PA=LU``.
"""

# ╔═╡ 208197f7-5691-480c-bad4-d4676d381a96
md"""
## Factorización  de Cholesky. Matrices simétricas y definidas positivas
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000000a
md"""
**Definición.** Se dice que una matriz ``A\in\mathbb{R}^{n\times n}`` es *simétrica y definida positiva*, o brevemente *sdp*, si
```math
A^T=A \qquad\text{y}\qquad x^TAx>0,\quad\forall x\in\mathbb{R}^n\setminus\{0\}.
```
Se dice *semi-definida positiva* si ``x^TAx\geq0``, ``\forall x\in\mathbb{R}^n``.

**Observación.** Vale la pena observar lo siguiente:

- ``x^TAx=(Ax)\cdot x``. Luego, si ``A`` es sdp el coseno del ángulo entre ``Ax`` y ``x`` es ``>0``, por lo que el ángulo entre ``x`` y ``Ax`` es agudo. Si pensamos en la transformación ``x\to Ax``, vemos que transforma ``x`` en otro vector ``Ax`` que apunta en una dirección parecida a ``x``.
- Si ``A`` es sdp, entonces ``A`` es invertible. En efecto, ``Ax=0`` implica ``x^TAx=0``, lo que implica ``x=0``, por lo que la única solución del problema homogéneo ``Ax=0`` es la trivial.
- Si ``A`` es sdp, entonces es diagonalizable, con todos sus autovalores reales y positivos. Es decir, existe ``U\in\mathbb{R}^{n\times n}`` ortogonal ``U^TU=I`` con ``A=U\Lambda U^T`` y ``\Lambda`` diagonal con todos los elementos de la diagonal positivos.
- Si ``R\in\mathbb{R}^{m\times n}``, entonces ``R^TR`` es simétrica y semi-definida positiva, pues ``x^TR^TRx=(Rx)^T(Rx)=\|Rx\|_2^2\geq0`` para todo ``x\in\mathbb{R}^n``.
- Si las ``n`` columnas de ``R`` son linealmente independientes, entonces ``R^TR`` es sdp (será necesario que ``n\leq m``).
- Si ``A`` es sdp y escribimos ``A=\begin{bmatrix}\alpha&a^T\\a&A_*\end{bmatrix}``, con ``\alpha\in\mathbb{R}``, ``a\in\mathbb{R}^{n-1}`` y ``A_*\in\mathbb{R}^{(n-1)\times(n-1)}``, entonces ``\alpha>0`` y ``A_*`` es sdp. En otras palabras, si ``A\in\mathbb{R}^{n\times n}`` es sdp, entonces ``a_{11}>0`` y ``A_{2:n,2:n}`` es sdp.
- Si ``A`` es sdp, ``a_{ii}=e_i^TAe_i>0``, ``i=1,2,\dots,n``. Es decir, si ``A`` es sdp, entonces todos los elementos de la diagonal son positivos.
- Sabemos que toda matriz se puede escribir como suma de una simétrica más una antisimétrica: ``A = \underbrace{\frac{A+A^T}{2}}_{\text{parte simétrica de }A} + \underbrace{\frac{A-A^T}{2}}_{\text{parte antisimétrica de }A}``. Como ``x^TAx\in\mathbb{R}``, resulta ``x^TAx=(x^TAx)^T=x^TA^Tx``, y por lo tanto ``x^T\big(\frac{A-A^T}{2}\big)x = -x^T\big(\frac{A-A^T}{2}\big)x``, de donde ``x^T\big(\frac{A-A^T}{2}\big)x=0``. En consecuencia
  ```math
  x^TAx = x^T\Big(\frac{A+A^T}{2}\Big)x + x^T\Big(\frac{A-A^T}{2}\Big)x = x^T\Big(\frac{A+A^T}{2}\Big)x.
  ```
  Por lo tanto, si una matriz ``A`` cumple ``x^TAx>0`` para todo ``x\in\mathbb{R}^n\setminus\{0\}``, luego su *parte simétrica* ``\frac{A+A^T}{2}`` es sdp.

**Teorema.** *Si ``A\in\mathbb{R}^{n\times n}`` es sdp, entonces existe ``R`` triangular superior (de ``n\times n``) tal que ``A=R^TR``.*

**Observación.** Vale la pena observar lo siguiente:

- La factorización de Cholesky es una factorización ``LU`` con ``L=R^T`` y ``U=R`` (``L`` no es necesariamente unitaria).
- La matriz ``R`` se llama *factor de Cholesky* de ``A``.
- En base a uno de los comentarios anteriores al enunciado, vemos que podríamos haber enunciado el teorema de la siguiente manera: *Sea ``A\in\mathbb{R}^{n\times n}``, entonces ``A`` es sdp si y solo si existe ``R`` triangular superior (no singular) tal que ``A=R^TR``.*

**Demostración.** La demostración es constructiva. Observemos que si existiera tal matriz ``R`` debería cumplir que
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000001b
md"""
```math
A = \begin{bmatrix}\alpha&a^T\\a&A_*\end{bmatrix} = \begin{bmatrix}\rho&0^T\\r&R_*^T\end{bmatrix}\begin{bmatrix}\rho&r^T\\0&R_*\end{bmatrix} = \begin{bmatrix}\rho^2&\rho r^T\\\rho r&R_*^TR_*+rr^T\end{bmatrix},
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000001c
md"""
con ``\alpha=a_{11}``, ``\rho\in\mathbb{R}``, ``r\in\mathbb{R}^{n-1}`` y ``A_*,R_*\in\mathbb{R}^{(n-1)\times(n-1)}``, con ``R_*`` triangular superior. Para que valga la igualdad se debería cumplir:

1. ``\rho^2=\alpha``;
2. ``\rho r=a``;
3. ``R_*^TR_*+rr^T=A_*``.

Para que se cumplan las dos primeras, basta con elegir ``\rho=\sqrt\alpha>0`` (recordemos que ``\alpha=a_{11}>0``) y ``r=\frac{1}{\rho}a``. Para que se cumpla la tercera, observemos que debemos hallar ``R_*\in\mathbb{R}^{(n-1)\times(n-1)}``, triangular superior, tal que
```math
R_*^TR_* = A_*-rr^T.
```
Si podemos demostrar que ``A_*-rr^T`` es sdp, el argumento se sigue por inducción.

Veamos que ``A_*-rr^T`` es sdp. Queremos ver que dado ``y\in\mathbb{R}^{n-1}\setminus\{0\}``, se cumple ``y^T(A_*-rr^T)y>0``, y para esto veremos si es cierto que existe ``x\in\mathbb{R}^n\setminus\{0\}`` tal que ``y^T(A_*-rr^T)y=x^TAx``. Proponemos ``x=\begin{pmatrix}\eta\\y\end{pmatrix}`` con ``\eta\in\mathbb{R}`` a determinar. Veamos:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000001d
md"""
```math
x^TAx = \begin{pmatrix}\eta&y^T\end{pmatrix}\begin{bmatrix}\alpha&a^T\\a&A_*\end{bmatrix}\begin{pmatrix}\eta\\y\end{pmatrix} = \begin{pmatrix}\eta&y^T\end{pmatrix}\begin{pmatrix}\alpha\eta+a^Ty\\a\eta+A_*y\end{pmatrix} = \alpha\eta^2+\eta a^Ty+\eta y^Ta+y^TA_*y = \alpha\eta^2+2\eta a^Ty+y^TA_*y.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000001e
md"""
Debemos elegir ``\eta\in\mathbb{R}`` para que este último término sea igual a
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000004a
md"""
```math
y^T(A_*-rr^T)y = y^TA_*y-y^Trr^Ty = y^TA_*y-(r^Ty)^2.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000004b
md"""
Así, nos queda planteada la siguiente ecuación para ``\eta``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000004c
md"""
```math
0 = \alpha\eta^2+2\eta a^Ty+y^TA_*y - [y^TA_*y-(r^Ty)^2] = \alpha\eta^2+2\eta a^Ty+(r^Ty)^2 = \rho^2\eta^2+2\eta\rho r^Ty+(r^Ty)^2 = (\rho\eta+r^Ty)^2.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000004d
md"""
Por lo tanto, tomando ``\eta=-r^Ty/\rho`` se cumple esa ecuación. Finalmente, tomando ``x=\begin{pmatrix}-r^Ty/\rho\\y\end{pmatrix}`` se cumple que ``x^TAx=y^T(A_*-rr^T)y>0`` porque ``A`` es sdp.

Por lo tanto, ``A_*-rr^T`` es sdp, y el teorema queda demostrado por inducción.

**Observación.** Dada ``A`` una matriz sdp, existe una *única* matriz ``R`` triangular superior, *con elementos de la diagonal positivos*, tal que ``A=R^TR``.

**Observación.** Si ``D`` es la matriz diagonal que en la diagonal tiene los mismos elementos que ``R``, entonces
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000004e
md"""
```math
A = R^TR = \underbrace{R^TD^{-1}}_{L}\underbrace{DR}_{U} = LU,
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000004f
md"""
con ``L`` triangular inferior unitaria y ``U`` triangular superior. Así, a partir de la factorización de Cholesky se puede obtener la factorización ``LU`` y viceversa: a partir de la factorización ``LU`` se puede obtener el factor de Cholesky.

**Observación.** El algoritmo `mi_cholesky!` de más abajo podría hacerse más eficiente si, al momento de actualizar la matriz, solo se actualiza la parte triangular superior de la submatriz restante. De esta manera el número de operaciones se reduce aproximadamente a la mitad. La factorización de Cholesky ya está programada eficientemente en Julia: la función correspondiente es `cholesky` (de `LinearAlgebra`).
"""

# ╔═╡ 613fd5e5-672c-46e2-b712-69a30759a5e2
"""
	cholesky_libro(A::AbstractMatrix)

Primera implementación (en desarrollo) de la factorización de Cholesky.
```julia
function cholesky_libro(A::AbstractMatrix)
	n = size(A, 1)

	R = zeros(size(A))
	for i = 1:n-1
		R[i,i] = sqrt(A[i,i])
		R[i, i+1:n] = A[i, i+1:n] / R[i,i]
		A[i+1:n, 1+i:n] = A[i+1:n, 1+i:n] - transpose(R[i,i+1:n]) * R[i,i+1:n]
	end

	R[n,n] = sqrt(A[n,n])

	return R
end
```
"""
function cholesky_libro(A::AbstractMatrix)
	n = size(A, 1)

	R = zeros(size(A))
	for i = 1:n-1
		R[i,i] = sqrt(A[i,i])
		R[i, i+1:n] = A[i, i+1:n] / R[i,i]
		A[i+1:n, 1+i:n] = A[i+1:n, 1+i:n] - transpose(R[i,i+1:n]) * R[i,i+1:n]
	end

	R[n,n] = sqrt(A[n,n])

	return R
end

# ╔═╡ 912bf55a-1bfa-4042-97f7-cace7748212e
"""
	mi_cholesky!(A::AbstractMatrix)

Calcula ``R`` triangular superior tal que ``A = R^T R`` y devuelve el resultado sobreescribiendo ``A``.
```julia
function mi_cholesky!(A::AbstractMatrix)
	n = size(A, 1)

	for c in 1:n-1
		ρ = sqrt(A[c,c])
		ρ⁻¹ = inv(ρ) #mult es mas eficiente que div (🍎).

		A[c,c] = ρ
		lmul!(ρ⁻¹, @view A[c,c+1:n]) #(🍎) y @view para modificar in-place
		A[c+1:n, c] .= 0

		@views BLAS.syr!('U', -1.0, A[c, c+1:n], A[c+1:n, c+1:n])
		#Rutina BLAS para A + αxxᵀ con A simetrica. syr!(uplo, α, x, A)
	end

	A[n,n] = sqrt(A[n,n])


	return A
end
```
"""
function mi_cholesky!(A::AbstractMatrix)
	n = size(A, 1)
	
	for c in 1:n-1
		ρ = sqrt(A[c,c])
		ρ⁻¹ = inv(ρ) #mult es mas eficiente que div (🍎).
		
		A[c,c] = ρ
		lmul!(ρ⁻¹, @view A[c,c+1:n]) #(🍎) y @view para modificar in-place
		A[c+1:n, c] .= 0

		@views BLAS.syr!('U', -1.0, A[c, c+1:n], A[c+1:n, c+1:n])
		#Rutina BLAS para A + αxxᵀ con A simetrica. syr!(uplo, α, x, A)
	end
	
	A[n,n] = sqrt(A[n,n])

	
	return A
end

# ╔═╡ 17d39397-96f4-4b4e-8456-6541eef6cdd8
md"""
## Ejercicios

### 2.1 Demostrar las siguientes afirmaciones

**Proposición:** Si ``A`` y ``B`` son matrices triangulares inferiores (unitarias) entonces ``AB`` es triangular inferior (unitaria).

**Demostración:**
Si ``A`` y ``B`` son triangulares inferiores unitarias, entonces:

- 👽 ``a_{ii}=b_{ii}=1``.
- 🐱 para la fila ``i`` se cumple ``a_{ik} = b_{ik} = 0`` si ``k > i``
- 🐶 para la columna ``i`` se cumple ``a_{ki} = b_{ki} = 0`` si ``k < i``

Sea ``C=AB`` de componentes ``c_{ij}``. Veamos que ``C`` es triangular inferior unitaria probando que ``c_{ik} = 0`` para ``k > i`` y que ``c_{ii} = 1``.

- ``c_{ii}=1``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000050
md"""
```math
c_{ii} = \sum_{k=1}^{n}a_{ik}b_{ki} =
	\underbrace{\sum_{k=1}^{i-1}a_{ik}b_{ki}}_{🍎} +
	a_{ii}b_{ii} +
	\underbrace{\sum_{k=i+1}^{n}a_{ik}b_{ki}}_{🍐} =
	a_{ii}b_{ii} =
	1 \qquad(\text{por 👽})
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000051
md"""
La suma 🍎 se anula pues, por 🐶, las componentes ``b_{ki}=0``, mientras que la suma 🍐 se anula pues, por 🐱, las componentes ``a_{ik}=0``.

- ``c_{ik} = 0`` para ``k > i``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000052
md"""
```math
c_{ik} = \sum_{j=1}^{n}a_{ij}b_{jk} =
	\underbrace{\sum_{j=1}^{i-1}a_{ij}b_{jk}}_{🍑} +
	a_{ii}b_{ik} +
	\underbrace{\sum_{j=i}^{n}a_{ij}b_{jk}}_{🍈} =
	a_{ii}b_{ik} =
	0 \qquad(\text{por 🐱})
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000053
md"""
En la suma 🍑, dado que ``j < i`` tenemos por 🐶 que ``b_{jk}=0``. Por otro lado, en la suma 🍈, dado que ``j > i``, tenemos por 🐱 que ``a_{ij}=0``.
"""

# ╔═╡ 412cf337-4ffe-4ec8-b363-88f98bcd92f3
md"""
---------------------------------------------------
**Proposición:** Si ``A`` es triangular inferior unitaria entonces ``A`` es invertible y ``A^{−1}`` es triangular inferior unitaria.

**Demostración** Si ``A`` es triangular inferior unitaria, entonces ``\text{det}(A)=\prod_{i=1}^{n}a_{ii}=\prod_{i=1}^{n}1=1``. Luego, ``A`` tiene que ser invertible.

Veamos que ``X=A^{-1}`` es triangular inferior. 
Supogamos que esto no es así. Es decir, que tenemos al menos una componente ``x_{ij}\neq 0`` con ``i < j``. Entonces si ``I`` es la matriz identidad de componentes ``y_{ij}`` se tiene

```math
	y_{ij} = \sum_{k=1}^{n}a_{ik}x_{kj} = \underbrace{\sum_{k=1}^{i-1}a_{ik}x_{kj}}_{🍌} + a_{ii}x_{ij} + \underbrace{\sum_{k=i+1}^{n}a_{ik}x_{kj}}_{🍓} = x_{ij}\neq 0
```

En la suma 🍌, dado que ``k < i`` (💢 *no se como seguir*)

En la suma 🍓, dado que ``k > i``, las componentes ``a_{ik}=0``, por lo que la suma es cero.
"""

# ╔═╡ c7e0143a-a496-45ce-8a97-5e64706ff11f
md"""
------------------------------
### 2.2

Comparar los tiempos de ejecución de las funciones `triinf`, `triinf_filas` y `triinf_cols` para la resolución de sistemas lineales con matrices triangulares inferiores. (En el libro se usan `tic`/`toc` de MATLAB/OCTAVE; en Julia usamos `@benchmark` de `BenchmarkTools`.)
"""

# ╔═╡ dff314f4-9454-403d-9b81-c7d883ed912e
begin
	N = 1000
	NB = 50
	A = rand(N, N)
	B = rand(N, NB)
	LT = tril(A)
	UT = triu(A)
	b = rand(N)
end;

# ╔═╡ 1b667954-532e-4d9e-8c7f-113b3a8aeab0
md"""
```julia
@benchmark triinf_libro(LT, b) #versión del libro
```
"""

# ╔═╡ 5d2cd4f2-7f8a-4b4f-a68d-9b011731da7a
@benchmark triinf_libro(LT, b)

# ╔═╡ 6a1b7410-9ccd-44bf-a9c9-649669fd01db
md"""
----------------------
```julia
@benchmark triinf_filas_libro(LT, b) #versión del libro
```
"""

# ╔═╡ 24a8a755-e2a4-4d8f-8836-db3866d957b7
@benchmark triinf_filas_libro(LT, b)

# ╔═╡ cfe530ba-c6a4-4637-aead-8f047681c930
md"""
---------------------
```julia
@benchmark triinf_filas_opt(LT, b) #version mia
```
"""

# ╔═╡ 65307b35-ffb2-4495-b616-abc7b6bf57a6
@benchmark triinf_filas_opt(LT, b)

# ╔═╡ 9ab6edc0-3f44-4fbc-a616-10698326b5c5
md"""
---------------------
```julia
@benchmark triinf_cols_libro(LT, b) #version del libro
```
"""

# ╔═╡ 8ab0812f-dc18-45a3-93db-7a984c4105ad
@benchmark triinf_cols_libro(LT, b)

# ╔═╡ 049e78f8-37b1-4df2-976f-bf4f3b8d55a7
md"""
---------------------
```julia
@benchmark triinf_cols_opt(LT, b) #version mia
```
"""

# ╔═╡ 5b446296-e0b0-4cc3-8be9-248af5889b4a
@benchmark triinf_cols_opt(LT, b)

# ╔═╡ 22db9bb0-b776-4950-9053-6ad8f697f953
md"""
### 2.3

Basándose en `triinf_cols`, crear una función `triinfm` que resuelva varios sistemas simultáneamente. Más precisamente, si ``L`` es una matriz triangular inferior y ``B`` es una matriz cuyas columnas representan varios términos independientes, entonces `x = triinfm(L, B)` resuelve los distintos sistemas donde la ``j``-ésima columna de `x` contiene la solución del sistema cuyo lado derecho es la ``j``-ésima columna de ``B``.
"""

# ╔═╡ 7a5af114-ad80-4d95-a50a-5103a72ef0bb
"""
	triinfm(L, B)

Si ``L`` es una matriz triangular inferior y ``B`` es una matriz cuyas columnas representan varios términos independientes, entonces `x =  triinfm(L, B)` resuelve
los distintos sistemas donde la ``j``-ésima columna de `x` contiene la solución del sistema cuyo lado derecho es la ``j``-ésima columna de ``b``.

```julia
function triinfm(L::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	R = similar(B)
	
	for (c, b) in enumerate(eachcol(B))
		R[:, c] .= triinf_cols_libro(L, b)
	end

	return R
end
```
"""
function triinfm(L::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	R = similar(B)

	for (c, b) in enumerate(eachcol(B))
		R[:, c] .= triinf_cols_libro(L, b)
	end

	return R
end

# ╔═╡ 988078da-20a6-4e98-b090-f22703793929
md"""
```julia
@benchmark triinfm(LT, B)
```
"""

# ╔═╡ c3e94ff5-17f1-44f3-b5be-e44fb58e973b
@benchmark triinfm(LT, B)

# ╔═╡ 96c9bd13-bd3b-4ece-9652-971b957b0bb8
"""
	triinfm2(L, B)

Si ``L`` es una matriz triangular inferior y ``B`` es una matriz cuyas columnas representan varios términos independientes, entonces `x =  triinfm(L, B)` resuelve
los distintos sistemas donde la ``j``-ésima columna de `x` contiene la solución del sistema cuyo lado derecho es la ``j``-ésima columna de ``b``.

```julia
function triinfm2(L::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	R = similar(B)
	
	for (c, b) in enumerate(eachcol(B))
		R[:, c] .= triinf_cols_opt(L, b)
	end

	return R
end
```
"""
function triinfm2(L::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	R = similar(B)

	for (c, b) in enumerate(eachcol(B))
		R[:, c] .= triinf_cols_opt(L, b)
	end

	return R
end

# ╔═╡ 25802262-54df-407f-8e57-cfa2ad9df3ea
md"""
```julia
@benchmark triinfm2(LT, B)
```
"""

# ╔═╡ 0e33798e-0be8-4d4d-abec-8867b6c24285
@benchmark triinfm2(LT, B)

# ╔═╡ f763ea35-a191-4d70-a475-5166fb28814c
"""
	triinfm_opt(L, B)

Si ``L`` es una matriz triangular inferior y ``B`` es una matriz cuyas columnas representan varios términos independientes, entonces `x =  triinfm(L, B)` resuelve
los distintos sistemas donde la ``j``-ésima columna de `x` contiene la solución del sistema cuyo lado derecho es la ``j``-ésima columna de ``b``.

```julia
function triinfm_opt(L::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    n = size(B, 1)
    cb = size(B, 2)
    X = copy(B)

    @inbounds begin
        inverso = 1 / L[1,1]
        @simd for c in 1:cb
            X[1, c] *= inverso
        end

        for i in 2:n
            invpiv = 1 / L[i, i]

            for c in 1:cb
                xi = X[i-1, c]

                @simd for k in i:n
                    X[k, c] -= xi * L[k, i-1]
                end
            end

            @simd for c in 1:cb
                X[i, c] *= invpiv
            end
        end
    end
    return X
end
```
"""
function triinfm_opt(L::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    n = size(B, 1)
    cb = size(B, 2)
    X = copy(B)
    @inbounds begin
        inverso = 1 / L[1,1]
        @simd for c in 1:cb
            X[1, c] *= inverso
        end

        for i in 2:n
            invpiv = 1 / L[i, i]
            for c in 1:cb
                xi = X[i-1, c]
                @simd for k in i:n
                    X[k, c] -= xi * L[k, i-1]
                end
            end
            @simd for c in 1:cb
                X[i, c] *= invpiv
            end
        end
    end
    return X
end

# ╔═╡ 4711e388-9d5c-4e4a-9c97-3e13c3b40b0d
md"""
```julia
@benchmark triinfm_opt(LT, A)
```
"""

# ╔═╡ 6db92a10-a6e4-4671-9303-481d57006cfb
@benchmark triinfm_opt(LT, B)

# ╔═╡ e5b81c9d-318e-4c4c-b21a-ae79689ad732
md"""
### 2.4

Basándose en `triinf_cols`, crear la función `trisup` para resolver sistemas lineales con matrices triangulares superiores.
"""

# ╔═╡ 16ee2bf0-695b-4603-92e7-fcf22e1252f2
"""
	trisup(U, b)


Resuelve el sistema ``Ux = b`` suponiendo ``U`` triangular superior cuadrada con todos los elementos de la diagonal no nulos.
	
Esta implementación reutiliza `triinf_cols_opt` luego de convertir ``U`` en triangular inferior. 

No es una mejor implementación que `trisup_cols_opt`.

```julia
function trisup(U::AbstractMatrix{T}, b::AbstractVector{T}) where T
	@views x = triinf_cols_opt(U[end:-1:1, end:-1:1], b[end:-1:1])
	reverse!(x)
end
```
"""
function trisup(U::AbstractMatrix{T}, b::AbstractVector{T}) where T
	@views x = triinf_cols_opt(U[end:-1:1, end:-1:1], b[end:-1:1])
	reverse!(x)
end

# ╔═╡ 16f0305a-ac93-479b-9192-07a426a62b84
md"""
```julia
@benchmark trisup(UT, b)
```
"""

# ╔═╡ 180f9366-0b80-4cdd-8c22-f2690023e172
@benchmark trisup(UT, b)

# ╔═╡ 3a060e4c-ddcf-45fb-99dd-174ca04eb578
md"""
--------------------------
```julia
@benchmark trisup_cols_opt(UT, b)
```
"""

# ╔═╡ 9664a5b5-8f43-4432-bdc0-55bdd922ee62
@benchmark trisup_cols_opt(UT, b)

# ╔═╡ de9e6fc9-9815-4cfe-a28b-b3cb4cfed971
md"""
### 2.5

Crear una función `tridilu(d, dmenos, dmas)` que calcule la factorización ``LU`` sin pivoteo de una matriz tridiagonal ``A`` donde
```math
A = \text{diag}(dmenos(2:n),-1) + \text{diag}(d) + \text{diag}(dmas(1:n-1),1).
```
Los vectores `l` y `u` que devuelve son tales que
```math
L = \text{diag}(\text{ones}(n,1)) + \text{diag}(l(2:n),-1), \qquad U = \text{diag}(u) + \text{diag}(dmas(1:n-1),1).
```
Además:

#### 2.5.a)
Realizar el conteo de operaciones realizadas por `tridilu` y comparar con la cantidad de operaciones realizadas por la función `eleu` desarrollada en clase.

#### 2.5.b)
Crear las funciones `triinfbidi` y `trisupbidi` para resolver sistemas triangulares bidiagonales, donde los argumentos de entrada sean vectores.

#### 2.5.c)
Considerar un sistema lineal asociado a la matriz ``A`` de ``N\times N`` (``N=100,1000,10000``) donde
```math
A = \begin{bmatrix}2&-1&&&\\-1&2&-1&&\\&\ddots&\ddots&\ddots&\\&&-1&2&-1\\&&&-1&2\end{bmatrix}.
```
Comparar los tiempos necesarios para hallar la solución del sistema utilizando las funciones `eleu`, `triinf_cols` y `trisup` por un lado, y `tridilu`, `triinfbidi` y `trisupbidi` por otro lado.
"""

# ╔═╡ 39640d3d-4a86-4deb-9c47-3b614c91b83d
let 
	N = 1000
	
	d = rand(N)
	dmenos = rand(N-1)
	dmas = rand(N-1)
	
	A = Tridiagonal(dmenos, d, dmas)

	B = Bidiagonal(d, dmenos, :L)

end

# ╔═╡ 4c8dac54-fd8e-47f1-b5fc-90c3f6140388
md"""
En Julia, las matrices creadas con el constructor `Tridiagonal` poseen cuatro campos que son los vectores `dl`, `d`, `du` y `du2`. Los primeros tres campos son las diagonales. Es decir, en lugar de alojar la matriz en memoria como lo hace con las matrices densas, en este caso aloja los vectores de las diagonales. 


Dado que solo hay que hacer ceros en `l`, los unicos multiplicadores que tendrá la matriz `L` serán en la misma diagonal. 

En lo que respecta a `dmas`, como la componente de arriba es cero, cualquier combinacion lineal no cambiará su valor, por lo que `dmas` también estará presente en `U`. La diagonal principal de ``U``, por otro lado, si variará. 
"""

# ╔═╡ 033f3162-ba41-447a-90d8-31d8005dbae1
"""
	tridilu(A)

Devuelve ``l`` y ``u`` que satisfacen las condiciones del problema 2.5
```julia
function tridilu(A::Tridiagonal)
	n = size(A, 1)
	l = Vector{Float64}(undef, n-1)
	u = Vector{Float64}(undef, n)

	u[1] = A.d[1]
	@inbounds for i in 1:n-1
		l[i] = A.dl[i] / u[i]
		u[i+1] = A.d[i+1] - l[i] * A.du[i]
	end

	return l, u
end
```
"""
function tridilu(A::Tridiagonal)
	n = size(A, 1)
	l = Vector{Float64}(undef, n-1)
	u = Vector{Float64}(undef, n)

	u[1] = A.d[1]
	@inbounds for i in 1:n-1
		l[i] = A.dl[i] / u[i]
		u[i+1] = A.d[i+1] - l[i] * A.du[i]
	end

	return l, u
end

# ╔═╡ 69bb0c36-8f68-44fe-978b-b430270f27ec
"""
	eleu(A::Tridiagonal)::Tuple{Bidiagonal, Bidiagonal}

Especialización de `eleu` para matrices tridiagonales.
```julia
function eleu(A::Tridiagonal)::Tuple{Bidiagonal, Bidiagonal}
	n = length(A.d)
	l = Vector{Float64}(undef, n-1)
	u = Vector{Float64}(undef, n)

	u[1] = A.d[1]
	@inbounds for i in 1:n-1
		l[i] = A.dl[i] / u[i]
		u[i+1] = A.d[i+1] - l[i] * A.du[i] #f ⟹ dₙ - m uₙ
	end

	return Bidiagonal(ones(n), l, :L), Bidiagonal(u, A.du, :U)
end
```
"""
function eleu(A::Tridiagonal)::Tuple{Bidiagonal, Bidiagonal}
	n = length(A.d)
	l = Vector{Float64}(undef, n-1)
	u = Vector{Float64}(undef, n)

	u[1] = A.d[1]
	@inbounds for i in 1:n-1
		l[i] = A.dl[i] / u[i]
		u[i+1] = A.d[i+1] - l[i] * A.du[i] #f ⟹ dₙ - m uₙ
	end

	return Bidiagonal(ones(n), l, :L), Bidiagonal(u, A.du, :U)
end

# ╔═╡ 95183e5d-804a-4275-a669-164e4698989b
"""
	bidiinf_lu(A::Bidiagonal)::Tuple{Bidiagonal, Diagonal}

Calcula ``L`` y ``U`` para matrices triangulares inferiores bidiagonales.
```julia
function bidiinf_lu(A::Bidiagonal)::Tuple{Bidiagonal, Diagonal}
	n = length(A.dv)
	l = Vector{Float64}(undef, n-1)
	u = copy(A.dv)

	@inbounds @simd for i in 1:n-1
		l[i] = A.ev[i] / A.dv[i]
	end

	return Bidiagonal(ones(n), l, :L), Diagonal(u)
end
```
"""
function bidiinf_lu(A::Bidiagonal)::Tuple{Bidiagonal, Diagonal}
	n = length(A.dv)
	l = Vector{Float64}(undef, n-1)
	u = copy(A.dv)

	@inbounds@simd for i in 1:n-1
		l[i] = A.ev[i] / A.dv[i]
	end

	return Bidiagonal(ones(n), l, :L), Diagonal(u)
end

# ╔═╡ 0a117286-ea8d-4deb-b9f2-ba0aa0624c56
"""
	bidisup_lu(A::Bidiagonal)::Tuple{Diagonal, Bidiagonal}

Calcula ``L`` y ``U`` para matrices triangulares superiores bidiagonales.
```julia
function bidisup_lu(A::Bidiagonal)::Tuple{Diagonal, Bidiagonal}
	return I(length(A.dv)), copy(A)
end
```
"""
function bidisup_lu(A::Bidiagonal)::Tuple{Diagonal, Bidiagonal}
	return I(length(A.dv)), copy(A)
end

# ╔═╡ dd9fd48e-d288-4286-ac2a-c435ce98a232
"""
	eleu(A::Bidiagonal)::Union{Tuple{Bidiagonal,Diagonal},
								Tuple{Diagonal,Bidiagonal}}

Especialización de `eleu` para matrices bidiagonales.

```julia
function eleu(A::Bidiagonal)
	if A.uplo == 'L' 
		return bidiinf_lu(A)
	else 
		return bidisup_lu(A)
	end
end
```
"""
function eleu(A::Bidiagonal)::Union{Tuple{Bidiagonal,Diagonal}, Tuple{Diagonal,Bidiagonal}}
	if A.uplo == 'L' 
		return bidiinf_lu(A)
	else 
		return bidisup_lu(A)
	end
end

# ╔═╡ 94c2a977-ff59-4f6a-ab5d-7efe244a4cb9
md"""
#### (b)
"""

# ╔═╡ 4c3f97fb-fe9d-42da-b94a-6f183de1379f
"""
	triinfbidi(dp::Vector, ds::Vector, b::Vector)::Vector

Resuelve sistemas triangulares inferiores bidiagonales ``Lx=b``, donde `dp` es la diagonal principal de ``L`` y `ds` es la otra.

```julia
function triinfbidi(dp::Vector, ds::Vector, b::Vector)::Vector
	x = similar(b)
	
	x[1] = b[1] / dp[1]
	for i in 2:length(b)
		x[i] = (b[i] - x[i-1] * ds[i-1]) / dp[i]
	end

	return x
end
```
"""
function triinfbidi(dp::Vector, ds::Vector, b::Vector)::Vector
	x = similar(b)
	
	x[1] = b[1] / dp[1]
	for i in 2:length(b)
		x[i] = (b[i] - x[i-1] * ds[i-1]) / dp[i]
	end

	return x
end

# ╔═╡ f59c7338-d904-43df-a5ed-1a0f5e74ae81
"""
	triinfbidi(A::Bidiagonal, b::Vector)::Vector

Resuelve sistemas triangulares inferiores bidiagonales ``Ax=b``.

```julia
function triinfbidi(A::Bidiagonal, b::Vector)::Vector
	triinfbidi(A.dv, A.ev, b)
end
```
"""
function triinfbidi(A::Bidiagonal, b::Vector)::Vector
	triinfbidi(A.dv, A.ev, b)
end

# ╔═╡ 13dcf3b3-5f78-4c45-8309-f94914836d78
"""
	trisupbidi(dp::Vector, ds::Vector, b::Vector)::Vector

Resuelve sistemas triangulares superiores bidiagonales ``Ux=b``, donde `dp` es la diagonal principal de ``U`` y `ds` es la otra.

```julia
function trisupbidi(dp::Vector, ds::Vector, b::Vector)::Vector
	x = similar(b)
	n = length(b)
	
	x[n] = b[n] / dp[n]
	for i in n-1:-1:1
		x[i] = (b[i] - x[i+1] * ds[i]) / dp[i]
	end

	return x
end
```
"""
function trisupbidi(dp::Vector, ds::Vector, b::Vector)::Vector
	x = similar(b)
	n = length(b)
	
	x[n] = b[n] / dp[n]
	for i in n-1:-1:1
		x[i] = (b[i] - x[i+1] * ds[i]) / dp[i]
	end

	return x
end

# ╔═╡ 22a0d787-d4c4-463f-a5fb-a7863f25e4e6
"""
	trisupbidi(U::Bidiagonal, b::Vector)::Vector

Resuelve sistemas triangulares superiores bidiagonales ``Ux=b``.

```julia
function trisupbidi(U::Bidiagonal, b::Vector)::Vector
	trisupbidi(A.dv, A.ev, b)
end
```
"""
function trisupbidi(U::Bidiagonal, b::Vector)::Vector
	trisupbidi(U.dv, U.ev, b)
end

# ╔═╡ f9ae1b19-0d13-4270-a661-ee6d0b83a555
md"""
#### c)

En este caso, para hacer las comparaciones, dado que uso despacho multiple, tendré que crear una instancia sparse y una instancia densa para cada ``A``.
"""

# ╔═╡ 6602a167-7928-474f-867c-bb1a649a0d35
begin
	A25c100 = Tridiagonal(fill(-1, 99), fill(2, 100), fill(-1, 99))
	A25c100d = Matrix(A25c100)

	A25c1000 = Tridiagonal(fill(-1, 999), fill(2, 1000), fill(-1, 999))
	A25c1000d = Matrix(A25c1000)

	A25c10000 = Tridiagonal(fill(-1, 9999), fill(2, 10000), fill(-1, 9999))
	A25c10000d = Matrix(A25c10000)
end;

# ╔═╡ 0dd68058-8b18-49fd-aa33-32d89c451ca5
md"""
```julia
@benchmark eleu_libro(A25c100d) #como es densa, no tiene version especializada.
```
"""

# ╔═╡ 58b95828-3b31-464f-9be3-beedc994b168
@benchmark eleu_libro(A25c100d)

# ╔═╡ 2ca9d7b3-158e-47a6-9dba-57cac95975f3
md"""
```julia
#equivalente a `tridilu`
@benchmark eleu(A25c100) #como es triangular, usa version especializada.
```
"""

# ╔═╡ 76d80087-dfa7-4f1d-9424-c4934583b5fa
@benchmark eleu(A25c100)

# ╔═╡ 54b8c65d-ce69-4d64-8239-1d82bef3ac99
md"""
### 2.6

Sean ``g_j=(0,\dots,0,m_{j+1,j},\dots,m_{n,j})^T`` para ``j=1,2,\dots,n-1`` vectores de Gauss. Sean ``G_j:=I-g_je_j^T`` las transformaciones de Gauss asociadas a ``g_j``, ``j=1,2,\dots,n-1``. Demostrar que:

#### a)
**Proposición:** Si ``G^j=I-g_je_j^T`` es una matriz de Gauss, entonces ``(G^j)^{-1} = I+g_je_j^T``

**Demostración**

Veamos que ``G^j(G^j)^{-1} = I``.

```math
\begin{align*}
	G^j(G^j)^{-1} &= (I - g_je_j^T) (I + g_je_j^T)\\
	&= I + g_je_j^T - g_je_j^T - g_je_j^Tg_je_j^T\\
	&= I - g_j(e_j^Tg_j)e_j^T\\
\end{align*}
```

El vector ``e_j^T`` solo es distinto de cero en la componente ``j``, mientras que ``g_j`` en la componente ``j`` vale cero. Así, el producto escalar ``e_j^Tg_j=0`` y tenemos
```math
\begin{align*}
	G^j(G^j)^{-1} &= I - g_j(0)e_j^T\\
	&= I\\
\end{align*}
```
que es lo que queríamos demostrar.
"""

# ╔═╡ 0cbbe993-cb61-47d2-a2b6-89b6182fd891
md"""
#### b)

**Proposición:** Si ``G_i`` es una matriz de Gauss entonces
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000001f
md"""
```math
\prod_{i=1}^{n-1} (G^i)^{-1} = I + \sum_{j=1}^{n-1}g_je_j^T =
	\begin{bmatrix}
		1 & 0 & 0 & \cdots & 0\\
		m_{21} & 1 & 0 & \cdots & 0\\
		\vdots & \vdots & \ddots & \ddots & \vdots\\
		m_{n1} & m_{n2} & m_{n3} & \cdots & 1
	\end{bmatrix} = L
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000020
md"""
**Demostración**

Veamos una demostración por inducción.
Para ``n = 3`` tenemos
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000054
md"""
```math
\begin{align}
	(G^1)^{-1}(G^2)^{-1} &=  (I + g_1e_1^T) (I + g_2e_2^T)\\
	&= I + g_2e_2^T + g_1e_1^T + g_1(e_1^Tg_2)e_2^T\\
	&= I + \sum_{j=1}^{2}g_je_j^T
\end{align}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000055
md"""
Veamos ahora que, si
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000056
md"""
```math
	\prod_{i=1}^{n-2} (G^i)^{-1} = I + \sum_{j=1}^{n-2}g_je_j^T \quad\text{entonces}\quad\prod_{i=1}^{n-1} (G^i)^{-1} = I + \sum_{j=1}^{n-1}g_je_j^T.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000057
md"""
Partamos del producto de la derecha.
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000058
md"""
```math
\begin{align}
\prod_{i=1}^{n-1} (G^i)^{-1} &= \prod_{i=1}^{n-2} (G^i)^{-1} (G^{n-1})^{-1}\\
&= \bigg(I + \sum_{j=1}^{n-2}g_je_j^T\bigg)(G^{n-1})^{-1}\\
&= \bigg(I + \sum_{j=1}^{n-2}g_je_j^T\bigg)(I + g_{n-1}e_{n-1}^T)\\
&= I + g_{n-1}e_{n-1}^T + \sum_{j=1}^{n-2}g_je_j^T + \sum_{j=1}^{n-2}g_je_j^T g_{n-1}e_{n-1}^T
\end{align}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000059
md"""
El segundo y tercer sumando forman ``\sum_{j=1}^{n-1}g_je_j^T``. Además, dado que ``j < n-1``, el producto ``e_j^Tg_{n-1}=0`` del cuarto sumando se anula, anulando todo el término. Así,
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005a
md"""
```math
\prod_{i=1}^{n-1} (G^i)^{-1} = I + \sum_{j=1}^{n-1}g_je_j^T
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005b
md"""
que es lo que queríamos demostrar.
"""

# ╔═╡ 01207ff1-14a5-4d27-901f-86db6f74abb9
md"""
### 2.7

**Proposición:** La matriz
```math
A = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
```
no tiene una factorización ``LU``.

**Demostración**

Para ello, veamos si
```math
A = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix} = 
	\begin{bmatrix} 1 & 0 \\ λ & 1 \end{bmatrix} 
	\begin{bmatrix} u_{11} & u_{12} \\ 0 & u_{22} \end{bmatrix}
```

De acá, se tiene que 

```math
\left\{\begin{aligned}
	0 &= u_{11} & 👈👽\\
	1 &= u_{12}\\
	1 &= \lambda u_{11} & 👈👽\\
	0 &= \lambda u_{12} + u_{22}
\end{aligned}\right.
```

Dado que ``u_{11} = 0``, no puede existir ``\lambda`` tal que ``1 = \lambda u_{11}``. Luego, no existe ``L`` y por lo tanto, no existe factorización ``LU``.
"""

# ╔═╡ 543ea581-2c03-482c-950a-abc7ff17e601
md"""
### 2.8

Considerar ``A=\begin{bmatrix}\delta&1\\1&1\end{bmatrix}`` y ``b=(1+\delta,2)``.

#### 2.8.a)
Hallar la solución exacta de ``Ax=b`` para ``\delta\neq1``.

#### 2.8.b)
Resolver el mismo sistema utilizando las funciones `eleu`, `triinf_cols` y `trisup`, para los valores de ``\delta=10^{-2},10^{-4},10^{-6},\dots,10^{-18}``. Armar una tabla con los valores de `x[1]` y `x[2]` obtenidos, y explicar qué sucede (Ayuda: pensar qué ocurre con los errores de redondeo al hacer la eliminación de Gauss).
"""

# ╔═╡ 71ae729c-4540-4b6b-9483-50e47bd5d542
"""
	calcula_solucion_28(δ)

Muestra la solución al sistema del 2.8 para un δ en particular.

```julia
function calcula_solucion_28(δ)
	A = [δ 1; 1 1]
	b = [1 + δ, 2]

	L, U = eleu_libro(A)
	y = triinf_cols_opt(L, b)

	return trisup_cols_opt(U, y)
end
```
"""
function calcula_solucion_28(delta::Real)
	A = [delta 1; 1 1]
	b = [1+delta, 2]

	L, U = eleu_libro(A)
	y = triinf_cols_opt(L, b)
	
	return trisup_cols_opt(U, y)
end

# ╔═╡ b367a1bc-d1f8-48cc-94d4-e638a25ff29c
for i in 2:2:18
	@info "Para δ = 10^(-$(i)):" calcula_solucion_28(10.0^(-i))
end

# ╔═╡ 47bdeba6-06fc-45a6-b1af-59c1a8172550
md"""
### 2.9

Resolver el sistema del Ejercicio 2.8 utilizando la función `eleupiv` en vez de `eleu`. Realizar la misma tabla y explicar qué ocurre.
"""

# ╔═╡ 6ce62a96-7dea-4533-96f5-eb93bdb94936
"""
	calcula_solucion_29(δ)

Muestra la solución al sistema del 2.9 para un δ en particular.

```julia
function calcula_solucion_29(δ::Real)
	A = [δ 1; 1 1]
	b = [1 + δ, 2]
	
	L, U, p = eleupiv_libro(A)
	y = triinf_cols_opt(L, b[p])
	
	return trisup_cols_opt(U, y)
end
```
"""
function calcula_solucion_29(δ::Real)
	A = [δ 1; 1 1]
	b = [1 + δ, 2]
	
	L, U, p = eleupiv_libro(A)
	y = triinf_cols_opt(L, b[p])
	
	return trisup_cols_opt(U, y)
end

# ╔═╡ 122ec0db-dea0-468e-aab8-21e80616bb1e
for i in 2:2:18
	@info "Para δ = 10^(-$(i)):" calcula_solucion_29(10.0^(-i))
end

# ╔═╡ 3e3a1935-9849-4e23-9d0e-aeb67b910a6c
md"""
### 2.10

Investigar cómo funciona la función `lu` de `LinearAlgebra` en Julia (en el libro: la función `lu` de MATLAB/OCTAVE). Comparar con la función `eleupiv`.

La documentación de `LinearAlgebra.lu` en Julia.
"""

# ╔═╡ e4545a47-d487-4f7d-8f7f-d6cbdcd2fcea
"""
"""
lu

# ╔═╡ 0ac08613-e31a-45b7-8e5d-9b3deb334da2
md"""
### 2.11
**Proposición:** Si ``A`` es *estrictamente diagonalmente dominante*, entonces es no singular.


**Demostración**

Veamos, por induccion, que reducir por Gauss produce siempre pivotes no nulos. Para ello, veremos que si un pivote es no nulo y la submatriz restante es e.e.d., entonces la submatriz que se produce al eliminar debajo de ese pivote también es e.e.d.

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
\begin{align*}
	P_1 &= \begin{bmatrix}
		a_{11} & a_{12} & a_{13} & \cdots & a_{1n}\\
		0 & a_{22} - a_{21}\dfrac{a_{12}}{a_{11}} & a_{23} - a_{21}\dfrac{a_{13}}{a_{11}} & \cdots & a_{2n} - a_{21}\dfrac{a_{1n}}{a_{11}}\\
		\vdots & \vdots & \ddots & \ddots & \vdots\\
		0 & a_{n2} - a_{n1}\dfrac{a_{12}}{a_{11}} & a_{n3} - a_{n1}\dfrac{a_{13}}{a_{11}} & \cdots & a_{nn} - a_{n1}\dfrac{a_{1n}}{a_{11}}\\
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
es tambien e.e.d. Para índices ``2\leq i\leq n``, sea ``D=\{2,\dots,n\}\setminus\{i\}`` (👉👽) y
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005c
md"""
```math
	👉🐔\qquad \lambda_i = |a_{i1}|/|a_{11}|
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005d
md"""
Como ``A`` es e.e.d. en la fila 1, ``|a_{i1}| < |a_{11}|``, así que ``\lambda_i \in [0,1)``.

Por la desigualdad triangular, acotamos por separado la diagonal y la suma de ``b``:
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-00000000005e
md"""
```math
\begin{align*}
	👉🐱&& |b_{ii}| &= \bigg|a_{ii} - a_{i1}\dfrac{a_{1i}}{a_{11}}\bigg| \stackrel{🐔}{\geq} |a_{ii}| - \lambda_i|a_{1i}|\\
	👉🐶&& \sum_{j\in D}|b_{ij}| &= \sum_{j\in D}\bigg|a_{ij} - a_{i1}\dfrac{a_{1j}}{a_{11}}\bigg| \stackrel{🐔}{\leq} \sum_{j\in D}|a_{ij}| + \lambda_i\sum_{j\in D}|a_{1j}|
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
	👉🍌\qquad |b_{ii}| - \sum_{j\in D}|b_{ij}| &\geq
		|a_{ii}| - \lambda_i|a_{1i}| - \bigg[\sum_{j\in D}|a_{ij}| + \lambda_i\sum_{j\in D}|a_{1j}|\bigg]\\
	&\geq \bigg[|a_{ii}| - \sum_{j\in D}|a_{ij}|\bigg] - \lambda_i\bigg[|a_{1i}| + \sum_{j\in D}|a_{1j}|\bigg]
\end{align}
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000061
md"""
El término ``|a_{1i}| + \sum_{j\in D}|a_{1j}|`` es la suma completa de la fila 1 sin el pivote (porque ``i\notin D``), ``\sum_{j\neq 1}|a_{1j}|`` 👉🍉. Usamos las hopótesis sobre la fila 1 de ``A`` para acotar el segundo corchete, y la de de la fila ``i`` de ``A`` (sumando todos los ``j\neq i``, incluyendo ``j=1``) para acotar el primero:
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

Demostrar que si ``A`` es e.e.d., el sistema ``Ax=b`` se puede resolver por Eliminación de Gauss sin intercambios de renglones.

Demostración directa del ejercicio anterior.
"""

# ╔═╡ 525ebc35-2e94-4648-ace3-699bd22c3995
md"""
### 2.13

Una matriz cuadrada ``A \in \mathbb{R}^{n\times n}`` se dice simétrica y definida positiva (sdp) si ``A^T = A`` y satisface ``x^TAx > 0`` para todo vector ``x \in \mathbb{R}^n``, ``x \neq 0``. Decimos que es semi-definida
positiva si ``x^T Ax \geq 0`` para todo ``x \in \mathbb{R}^n``. Demostrar:

#### 2.13.(a)
**Proposición:** Si ``R\in \mathbb{R}^{m\times n}`` entonces ``R^TR`` es simétrica y semidefinida positiva.

**Demostración**
- Para ver que es simétrica, notemos que ``(R^TR)^T=((R)^T(R^T)^T)=(R^TR)``. Si ``A=R^TR``, podemos ver que ``A^T=A``.

- Para ver que es sdp tenemos que ver que ``\forall x\neq \mathbf{0}`` se cumple ``x^T(R^TR)x \geq 0``

```math
x^T(R^TR)x = (x^TR^T)(Rx) = (Rx)^T(Rx) = |Rx|^2 \geq 0.
```
"""

# ╔═╡ 4b81913b-0993-40f9-a315-7a0d153b6be1
md"""
#### 2.13.(b)
**Proposición:** Si las columnas de ``R \in \mathbf{R}^{m\times n}`` son linealmente independientes, entonces ``R^TR`` es sdp.

**Demostración**

Si las columnas de ``R`` son *li* entonces sabemos que ``Rx=\mathbf{0}`` solo si ``x=\mathbf{0}``.

Para probar que ``R^TR`` es sdp tomemos ``x\neq\mathbf{0}``. Entonces

```math
x^T R^T R x = |Rx| > 0
```
pues ``Rx\neq \mathbf{0}`` para ``x\neq\mathbf{0}``.
"""

# ╔═╡ 0a228ffd-476b-4dbf-8539-27b2dfb7c763
md"""
#### 2.13.(c).(i)

**Proposición:** Si ``A`` es sdp, entonces ``A`` es invertible.

**Demostración:**


Notemos que

```math
Ax = \mathbf{0} \implies x^TAx = x^T\mathbf{0} = 0
```

Si ``A`` es sdp entonces ``x^TAx > 0`` siempre que ``x\neq\mathbf{0}``, por lo tanto ``x`` solo puede ser ``\mathbf{0}`` para que ``Ax=\mathbf{0}``. Luego, ``A`` es no singular.
"""

# ╔═╡ d8d480fe-0d4c-44ca-8550-08f693d2a406
md"""
#### 2.13.(c).(ii)

**Proposición:** Si ``A`` es sdp, entonces todos los autovalores de ``A`` son positivos.

**Demostración:**

Como ``A`` es invertible existen ``\lambda`` y ``v\neq \mathbf{0}`` tales que

```math
Av = \lambda v \implies v^TAv = v^T\lambda v \implies v^T\lambda v > 0 \implies \lambda v^Tv > 0 \implies \lambda |v|^2 >0.
```

Para que este producto sea mayor que cero, es necesario que ``\lambda > 0``.

Luego, los autovalores de ``A`` tienen que ser positivos. 
"""

# ╔═╡ aeaf3d9c-c80e-403f-aac4-5016e38a50a5
md"""
#### 2.13.(c).(iii)

**Proposición:** Si ``A`` es sdp, entonces todos los elementos de la diagonal de ``A`` son positivos.

**Demostración:**

Tomemos ``x=e_i\neq \mathbf{0}``. Entonces, dado que ``A`` es sdf,
```math
	0 < e_i^TAe_i = e_i^Ta_{1:n,i} = a_{ii}
```

Como esto lo podemos hacer con cualquier vector canónico, tenemos lo que queríamos probar.
"""

# ╔═╡ 1bbc01c9-a77b-4771-9267-4903cae56ae0
md"""
#### 2.13.(c).(iv)

**Proposición:** Si ``A`` es sdp, entonces ``a^2_{ij} < a_{ii}a_{jj}`` para cada ``i \neq j``.

**Demostración:**

Sea ``v`` un vector con ceros, salvo en la componente ``i`` y ``j``.
"""

# ╔═╡ db3438f3-86bc-4f13-acfd-3d74c4f9c41d
md"""
```math
	0 < v^TAv = v^T(Av) = 
	v^T \begin{bmatrix}kkkk 
	\end{bmatrix}
```
"""

# ╔═╡ e2a04418-1bbb-4762-aae9-23199b5452b6
md"""
#### 2.13.(c).(v)

**Proposición:** Si ``A`` es sdp, entonces ``\max_{1\leq i,j\leq n}|a_{ij}|\leq \max_{1\leq i \leq n}|a_{ii}|`` 

**Demostración:**


"""

# ╔═╡ 20c0492a-79f9-4663-889f-262bf78ab3c8
md"""
#### 2.13.(d)

**Proposición:** Si ``A`` es sdp y escribimos ``A=[\alpha, a^T; a, A_*]``, con ``\alpha \in \mathbb{R}`` y ``A_*\in\mathbb{R}^{(n-1)\times(n-1)}``, entonces ``\alpha > 0`` y ``A_*`` es sdp.

**Demostración:**

Si ``x=e_i``, entonces ``0 < x^TAx=e_i^tAe_i=a_{ii}=\alpha``.

Veamos que ``A_*`` es sdp. Para ello, veamos que para ``y\neq \mathbf{0}`` se cumple ``y^TA_*y > 0``. Sea ``y`` y sea ``x  = (0, y)``. Como ``A`` es sdp, se cumple
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000027
md"""
```math
0 < x^TAx =
\begin{pmatrix}
	0 & y^T
\end{pmatrix}
\begin{bmatrix}
	\alpha &  a^T\\
	a & A_*
\end{bmatrix}
\begin{pmatrix}
	0 \\
	y
\end{pmatrix} =
\begin{pmatrix}
	0 & y^T
\end{pmatrix}
\begin{pmatrix}
	a^Ty \\
	A_*y
\end{pmatrix} =
y^TA_*y.
```
"""

# ╔═╡ e1a2b3c4-1111-4a1b-8c2d-000000000028
md"""
Como ``y`` era un vector arbitrario distinto del nulo, queda demostrado.
"""

# ╔═╡ a7123cea-ffbe-4fb2-8a09-146b35c761a9
md"""
### 2.14

Programar una versión recursiva de la factorización de Cholesky. Comparar el funcionamiento de ésta con el de la función `mi_cholesky!` dada en clase, para matrices de diferentes órdenes de la forma
```math
A = \begin{bmatrix}2&-1&&&\\-1&2&-1&&\\&\ddots&\ddots&\ddots&\\&&-1&2&-1\\&&&-1&2\end{bmatrix}.
```
"""

# ╔═╡ 7fc77d1e-294a-4796-8166-42836224d494
"""
	mi_cholesky(M::Tridiagonal; niv::Int)::Bidiagonal{:U}

Factorización de Cholesky para matrices tridiagonales.
```julia
function mi_cholesky!(M::Tridiagonal; niv = 1)
	if isequal(niv, size(M, 1))
		M.d[niv] = sqrt(M.d[niv])
		fill!(M.dl, 0)
		return Bidiagonal(M, :U)
	end
	
	ρ = sqrt(M.d[niv])
	M.d[niv] = ρ
	M.du[niv] = M.du[niv] / ρ
	M.d[niv+1] -= M.du[niv] ^2

	return mi_cholesky!(M, niv=niv+1)
end
```
"""
function mi_cholesky!(M::Tridiagonal; niv = 1)
	if isequal(niv, size(M, 1))
		M.d[niv] = sqrt(M.d[niv])
		fill!(M.dl, 0)
		return Bidiagonal(M, :U)
	end
	
	ρ = sqrt(M.d[niv])
	M.d[niv] = ρ
	M.du[niv] = M.du[niv] / ρ
	M.d[niv+1] -= M.du[niv] ^2

	return mi_cholesky!(M, niv=niv+1)
end

# ╔═╡ a8a4ea8d-7bb1-4734-914c-89b0f098dd21
let
	A = rand(100,100)
	Sd = A' * A + I

	@benchmark mi_cholesky!(S) setup=(S=copy($Sd)) evals=1
end

# ╔═╡ a0a3133c-2a28-4894-b68e-4c6bc5e7aa0e
"""
	mi_cholesky(A::AbstractMatrix)

Calcula ``R`` triangular superior tal que ``A = R^T R``.
```julia
function mi_cholesky(A::AbstractMatrix)
	x = copy(A)
	mi_cholesky!(x)

	return x
end
```
"""
function mi_cholesky(A::AbstractMatrix)
	x = copy(A)
	mi_cholesky!(x)

	return x
end

# ╔═╡ cc685668-aff3-4bdc-8ca3-62b94540532c
let
	A = rand(100,100)
	Sd = A' * A + I

	@benchmark mi_cholesky($Sd)
end

# ╔═╡ 6f3c1254-49cd-44a4-b7d5-0695acd52845
"""
	mi_cholesky2!(A::Tridiagonal)

Factorización de Cholesky para matrices tridiagonales, versión sin recursividad.
```julia
function mi_cholesky2!(A::Tridiagonal) #sin recursividad
	#fieldnames(Tridiagonal) ▶
		# A.d, A.du y A.dl son vectores diagonales.
		# En memoria componentes adyacentes.
	n = size(A, 1)

	@inbounds for c in 1:n-1
		ρ = sqrt(A[c,c])
		A.d[c] = ρ
		A.du[c] *= inv(ρ)
		A.d[c+1] -= A.du[c] ^ 2
	end

	fill!(A.dl, 0) #@simd
	A.d[n] = sqrt(A.d[n])

	return Bidiagonal(A, :U)
end
```
"""
function mi_cholesky2!(A::Tridiagonal) #sin recursividad
	#fieldnames(Tridiagonal) ▶
		# A.d, A.du y A.dl son vectores diagonales. 
		# En memoria componentes adyacentes.
	n = size(A, 1)
	
	@inbounds for c in 1:n-1
		ρ = sqrt(A[c,c])
		A.d[c] = ρ
		A.du[c] *= inv(ρ)
		A.d[c+1] -= A.du[c] ^ 2
	end

	fill!(A.dl, 0) #@simd
	A.d[n] = sqrt(A.d[n])

	return Bidiagonal(A, :U)
end

# ╔═╡ c6f04038-e616-431b-94a9-fdb6216b0796
begin 
	Nr_cholesky = 1000
	M = diagm(  0 => fill(2, Nr_cholesky), 
			  	 1 => fill(-1, Nr_cholesky - 1),
			 	 -1 => fill(-1, Nr_cholesky - 1)) .|> 
		Float64 |> 
		Tridiagonal

	@benchmark mi_cholesky!(Mm) setup=(Mm=copy($M)) evals=1
end

# ╔═╡ 23ae6e67-245e-48f9-8ce2-13b4a88c4ae7
begin 
	Nr2_cholesky = 1000
	M2 = diagm(  0 => fill(2, Nr_cholesky), 
			  	 1 => fill(-1, Nr_cholesky - 1),
			 	 -1 => fill(-1, Nr_cholesky - 1)) .|> 
		Float64 |> 
		Tridiagonal

	@benchmark mi_cholesky2!(Mm) setup=(Mm=copy($M2)) evals=1
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SIMD = "fdea26ae-647d-5447-a871-4b548cad5224"

[compat]
BenchmarkTools = "~1.8.0"
PlutoUI = "~0.7.83"
SIMD = "~3.7.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "8ce51ffd8ea30ef38cc7d3501c0ba6b332f0c1da"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

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

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "88352712893ec50bee3680605891eaf0e9ed6368"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.8.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

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

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "663e8b48b789916221e0765393b289ca6c88f24e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.0"

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

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
# ╟─43dcd902-9b72-11f1-b709-bf3a6b7c9b3e
# ╟─30461d29-5bf1-48d6-b7bd-6ec383304b64
# ╟─387e2688-5e2d-4b9b-991a-f99904f71c8c
# ╟─84ddd9c4-19a8-4dda-9464-764722d8e553
# ╟─8f839c88-9f4d-4e40-9349-8e7d3fb9b68e
# ╟─a0600810-db07-4e5b-8fc3-c56d0711dc13
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000001
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000002
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000000b
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000000c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000000d
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000000e
# ╟─cec62a0c-c997-4196-b067-298c19df73d0
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000003
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000000f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000010
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000011
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000012
# ╟─3eeb588f-a98b-4206-9fb0-e9aab80a58a0
# ╟─f26a9cf4-b22d-488c-b958-7ea388490af7
# ╟─b1b0d35f-ee7f-4e22-9bb8-350d1d00363f
# ╟─5d6fa8d5-d0f8-4cfb-b0ab-5b351ab29081
# ╟─885e65ff-16ff-4dd3-a759-2533a6a7a635
# ╟─83c191fc-a03f-4d35-b7c6-fba72ee03aed
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000004
# ╟─62816356-cfab-45ff-88d9-ab384e26473f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000005
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000013
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000014
# ╟─f7b1ae60-b108-4860-abfa-8ffc6438a5d4
# ╟─5d7e0747-f9c9-4794-98d0-dbbbb0d05e0a
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000006
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000029
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000002a
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000002b
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000002c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000002d
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000002e
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000002f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000030
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000031
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000032
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000033
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000034
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000035
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000036
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000037
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000038
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000039
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000003a
# ╟─a9f741f4-9e65-4fa6-b286-3a621aab6633
# ╟─2d6ac9d5-8423-4e74-9abd-3d192e1668fa
# ╟─51302fe8-2eb6-4e45-aea9-b12d7a434481
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000007
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000015
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000016
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000017
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000018
# ╟─2b09d162-6b0a-45d5-b966-119c5be4b476
# ╟─ab7c1a16-c76e-49d9-be60-feb360451a84
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000008
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000019
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000001a
# ╟─faa6b8f4-656c-4f7c-a761-e7842c89ad29
# ╟─ccebc4a5-b633-46f8-8153-e88cd0dbf717
# ╟─854bbb92-5282-4f9f-9373-a764aeae135a
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000009
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000003b
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000003c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000003d
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000003e
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000003f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000040
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000041
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000042
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000043
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000044
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000045
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000046
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000047
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000048
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000049
# ╟─208197f7-5691-480c-bad4-d4676d381a96
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000000a
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000001b
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000001c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000001d
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000001e
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000004a
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000004b
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000004c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000004d
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000004e
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000004f
# ╟─613fd5e5-672c-46e2-b712-69a30759a5e2
# ╟─912bf55a-1bfa-4042-97f7-cace7748212e
# ╟─17d39397-96f4-4b4e-8456-6541eef6cdd8
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000050
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000051
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000052
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000053
# ╟─412cf337-4ffe-4ec8-b363-88f98bcd92f3
# ╟─c7e0143a-a496-45ce-8a97-5e64706ff11f
# ╟─dff314f4-9454-403d-9b81-c7d883ed912e
# ╟─1b667954-532e-4d9e-8c7f-113b3a8aeab0
# ╟─5d2cd4f2-7f8a-4b4f-a68d-9b011731da7a
# ╟─6a1b7410-9ccd-44bf-a9c9-649669fd01db
# ╟─24a8a755-e2a4-4d8f-8836-db3866d957b7
# ╟─cfe530ba-c6a4-4637-aead-8f047681c930
# ╟─65307b35-ffb2-4495-b616-abc7b6bf57a6
# ╟─9ab6edc0-3f44-4fbc-a616-10698326b5c5
# ╟─8ab0812f-dc18-45a3-93db-7a984c4105ad
# ╟─049e78f8-37b1-4df2-976f-bf4f3b8d55a7
# ╟─5b446296-e0b0-4cc3-8be9-248af5889b4a
# ╟─22db9bb0-b776-4950-9053-6ad8f697f953
# ╟─7a5af114-ad80-4d95-a50a-5103a72ef0bb
# ╟─988078da-20a6-4e98-b090-f22703793929
# ╟─c3e94ff5-17f1-44f3-b5be-e44fb58e973b
# ╟─96c9bd13-bd3b-4ece-9652-971b957b0bb8
# ╟─25802262-54df-407f-8e57-cfa2ad9df3ea
# ╟─0e33798e-0be8-4d4d-abec-8867b6c24285
# ╟─f763ea35-a191-4d70-a475-5166fb28814c
# ╟─4711e388-9d5c-4e4a-9c97-3e13c3b40b0d
# ╟─6db92a10-a6e4-4671-9303-481d57006cfb
# ╟─e5b81c9d-318e-4c4c-b21a-ae79689ad732
# ╟─16ee2bf0-695b-4603-92e7-fcf22e1252f2
# ╟─16f0305a-ac93-479b-9192-07a426a62b84
# ╟─180f9366-0b80-4cdd-8c22-f2690023e172
# ╟─3a060e4c-ddcf-45fb-99dd-174ca04eb578
# ╟─9664a5b5-8f43-4432-bdc0-55bdd922ee62
# ╟─de9e6fc9-9815-4cfe-a28b-b3cb4cfed971
# ╟─39640d3d-4a86-4deb-9c47-3b614c91b83d
# ╟─4c8dac54-fd8e-47f1-b5fc-90c3f6140388
# ╟─033f3162-ba41-447a-90d8-31d8005dbae1
# ╟─69bb0c36-8f68-44fe-978b-b430270f27ec
# ╟─95183e5d-804a-4275-a669-164e4698989b
# ╟─0a117286-ea8d-4deb-b9f2-ba0aa0624c56
# ╟─dd9fd48e-d288-4286-ac2a-c435ce98a232
# ╟─94c2a977-ff59-4f6a-ab5d-7efe244a4cb9
# ╟─4c3f97fb-fe9d-42da-b94a-6f183de1379f
# ╟─f59c7338-d904-43df-a5ed-1a0f5e74ae81
# ╟─13dcf3b3-5f78-4c45-8309-f94914836d78
# ╟─22a0d787-d4c4-463f-a5fb-a7863f25e4e6
# ╟─f9ae1b19-0d13-4270-a661-ee6d0b83a555
# ╠═6602a167-7928-474f-867c-bb1a649a0d35
# ╟─0dd68058-8b18-49fd-aa33-32d89c451ca5
# ╟─58b95828-3b31-464f-9be3-beedc994b168
# ╟─2ca9d7b3-158e-47a6-9dba-57cac95975f3
# ╟─76d80087-dfa7-4f1d-9424-c4934583b5fa
# ╟─54b8c65d-ce69-4d64-8239-1d82bef3ac99
# ╟─0cbbe993-cb61-47d2-a2b6-89b6182fd891
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000001f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000020
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000054
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000055
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000056
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000057
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000058
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000059
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005a
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005b
# ╟─01207ff1-14a5-4d27-901f-86db6f74abb9
# ╟─543ea581-2c03-482c-950a-abc7ff17e601
# ╟─71ae729c-4540-4b6b-9483-50e47bd5d542
# ╟─b367a1bc-d1f8-48cc-94d4-e638a25ff29c
# ╟─47bdeba6-06fc-45a6-b1af-59c1a8172550
# ╟─6ce62a96-7dea-4533-96f5-eb93bdb94936
# ╟─122ec0db-dea0-468e-aab8-21e80616bb1e
# ╟─3e3a1935-9849-4e23-9d0e-aeb67b910a6c
# ╟─e4545a47-d487-4f7d-8f7f-d6cbdcd2fcea
# ╟─0ac08613-e31a-45b7-8e5d-9b3deb334da2
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000021
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000022
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000023
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000024
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000025
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000026
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005c
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005d
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005e
# ╟─e1a2b3c4-1111-4a1b-8c2d-00000000005f
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000060
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000061
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000062
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000063
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000064
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000065
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000066
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000067
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000068
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000069
# ╟─c6f5fd8d-505b-4b55-bd73-e5c028b771e1
# ╟─525ebc35-2e94-4648-ace3-699bd22c3995
# ╟─4b81913b-0993-40f9-a315-7a0d153b6be1
# ╟─0a228ffd-476b-4dbf-8539-27b2dfb7c763
# ╟─d8d480fe-0d4c-44ca-8550-08f693d2a406
# ╟─aeaf3d9c-c80e-403f-aac4-5016e38a50a5
# ╟─1bbc01c9-a77b-4771-9267-4903cae56ae0
# ╟─db3438f3-86bc-4f13-acfd-3d74c4f9c41d
# ╟─e2a04418-1bbb-4762-aae9-23199b5452b6
# ╟─20c0492a-79f9-4663-889f-262bf78ab3c8
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000027
# ╟─e1a2b3c4-1111-4a1b-8c2d-000000000028
# ╟─a7123cea-ffbe-4fb2-8a09-146b35c761a9
# ╟─7fc77d1e-294a-4796-8166-42836224d494
# ╟─a8a4ea8d-7bb1-4734-914c-89b0f098dd21
# ╟─a0a3133c-2a28-4894-b68e-4c6bc5e7aa0e
# ╟─cc685668-aff3-4bdc-8ca3-62b94540532c
# ╟─6f3c1254-49cd-44a4-b7d5-0695acd52845
# ╟─c6f04038-e616-431b-94a9-fdb6216b0796
# ╟─23ae6e67-245e-48f9-8ce2-13b4a88c4ae7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
