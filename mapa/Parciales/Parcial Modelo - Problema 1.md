---
tags: [parcial, capitulo4]
---

# Parcial Modelo - Problema 1 (15 pts)

## Enunciado

> Consideremos la matriz
$$
 M = \begin{bmatrix}0.8&1&0\\0&-0.8&1\\0&0&0.8\end{bmatrix}.
$$
> **(a)** Calcule las siguientes normas de $M$: $\|M\|_1$, $\|M\|_2$, $\|M\|_\infty$. (Puede usar la computadora para calcular alguna de ellas.)
>
> **(b)** ¿A qué converge la iteración $x_{k+1}=Mx_k$, comenzando de $x_0=(1,1,1)^T$?
>
> **(c)** ¿Se anima a decir a qué converge la iteración $x_{k+1}=Mx_k$ a partir de cualquier $x_0$? Fundamente.

## Resolución

### (a) Normas de $M$

Para $\|M\|_1$ y $\|M\|_\infty$ no hace falta la computadora: son las [[Normas de vectores y matrices|fórmulas explícitas]] de las normas matriciales inducidas por $\|\cdot\|_1$ y $\|\cdot\|_\infty$ en $\mathbb{R}^n$ -  $\|A\|_1$ es la máxima suma de **columnas** en valor absoluto, $\|A\|_\infty$ la máxima suma de **filas** (demostradas en [[4.5 - Formulas normas matriciales|Problema 4.5]], no hay que redemostrarlas, solo aplicarlas):
$$
\|M\|_1 = \max\{\underbrace{0.8}_{\text{col 1}},\ \underbrace{1+0.8}_{\text{col 2}},\ \underbrace{1+0.8}_{\text{col 3}}\} = 1.8, \qquad \|M\|_\infty = \max\{\underbrace{0.8+1}_{\text{fila 1}},\ \underbrace{0.8+1}_{\text{fila 2}},\ \underbrace{0.8}_{\text{fila 3}}\} = 1.8.
$$
$\|M\|_2=\sqrt{\lambda_{\max}(M^TM)}$ sí tiene forma cerrada en general, pero para esta $M$ (ni simétrica ni normal) no simplifica a nada más corto que calcular los autovalores de $M^TM$ a mano - de ahí que el enunciado **habilite explícitamente** usar la computadora para esta parte: es la única de las tres que conviene delegar.
```julia
using LinearAlgebra
M = [0.8 1.0 0.0; 0.0 -0.8 1.0; 0.0 0.0 0.8]
opnorm(M, 1), opnorm(M, 2), opnorm(M, Inf)
```
Da $\|M\|_1=1.8$ y $\|M\|_\infty=1.8$ (coincide con lo calculado a mano, confirma que no hay error de cuentas) y $\|M\|_2\approx1.61802$.

### (b) Convergencia de $x_{k+1}=Mx_k$ desde $x_0=(1,1,1)^T$

Primero, reescribir la iteración sin recursión: aplicando la fórmula $k+1$ veces, $x_{k+1}=Mx_k=M(Mx_{k-1})=\dots=M^{k+1}x_0$. Entonces la pregunta "¿a qué converge $x_k$?" es la pregunta "¿a qué converge $M^kx_0$?", y para eso hace falta saber a qué converge $M^k$ como matriz.

$M$ es **triangular superior**, y los autovalores de una matriz triangular son exactamente los elementos de su diagonal (demostrado en [[4.6 - Autovalores de matriz triangular|Problema 4.6]] - no es una coincidencia que $M$ sea triangular en el enunciado, es lo que hace posible leer $\sigma(M)$ sin calcular nada): $\sigma(M)=\{0.8,-0.8,0.8\}$, así que $\rho(M)=\max|\lambda|=0.8$.

Por qué $\rho(M)<1$ implica $M^k\to0$ (la matriz, no solo un vector): en [[Espectro y radio espectral]] se prueba la identidad $\rho(A)=\lim_{k\to\infty}\|A^k\|^{1/k}$ para cualquier norma inducida. Como $\rho(M)=0.8<1$, para $k$ suficientemente grande $\|M^k\|^{1/k}<1$ también (se acerca a $0.8$), es decir $\|M^k\|<1$ y de hecho $\|M^k\|\to0$ geométricamente. Entonces, para **cualquier** vector $x_0$ (en particular $(1,1,1)^T$):
$$
\|x_{k+1}\| = \|M^{k+1}x_0\| \leq \|M^{k+1}\|\,\|x_0\| \xrightarrow[k\to\infty]{} 0.
$$
Luego $\lim_k x_k=0$.

### (c) Convergencia desde cualquier $x_0$

El argumento de (b), leído de nuevo, **no usó en ningún paso que $x_0=(1,1,1)^T$** - la desigualdad $\|x_{k+1}\|\leq\|M^{k+1}\|\|x_0\|\to0$ vale para cualquier $x_0$ fijo, porque $\|M^{k+1}\|\to0$ no depende de $x_0$ en absoluto. Esa es toda la "fundamentación" que pide el enunciado: el límite $0$ no es una propiedad del punto de partida $(1,1,1)^T$, es una propiedad de $M$ sola ($\rho(M)<1$).

Dicho con el marco general del Capítulo 4 (para conectarlo con la teoría ya vista, no porque haga falta un argumento distinto): $x_{k+1}=Mx_k$ es la iteración estacionaria $x_{k+1}=Mx_k+c$ con $c=0$ ([[Metodos iterativos estacionarios - marco general]]), y el [[Teorema de convergencia de metodos iterativos estacionarios|teorema de convergencia]] de esa familia dice exactamente: converge **para todo $x_0$** si y solo si $\rho(M)<1$. Es el mismo resultado que (b) ya probó a mano, aplicado ahora explícitamente "para todo $x_0$" en vez de a uno fijo.

## Ver también
- [[Parcial Modelo - Indice|Parcial Modelo]]
- [[Normas de vectores y matrices]]
- [[Espectro y radio espectral]]
- [[4.6 - Autovalores de matriz triangular|Problema 4.6]]
