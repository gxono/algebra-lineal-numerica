---
tags: [capitulo5, teorema]
---

# Lema 5.43 (breakdown = buena noticia)

Sea $A$ no singular y sean $v_1,v_2,\dots$ los vectores generados por el [[GMRes - implementacion y rotaciones de Givens|algoritmo de Arnoldi]]. Sea $i$ el menor entero para el cual
$$
z_{i+1} = Av_i-\sum_{j=1}^i\big((Av_i)^Tv_j\big)v_j = 0.
$$
Entonces $x:=A^{-1}b\in x_0+\mathcal K_i$ - es decir, el espacio de Krylov actual ya contiene la solución exacta.

> [!demostracion]-
> **Paso 1: $\mathcal K_i$ es $A$-invariante.** Por construcción de Arnoldi, $Av_j\in\mathcal K_{j+1}\subset\mathcal K_i$ para $j=1,\dots,i-1$. Como $z_{i+1}=0$, también $Av_i=\sum_{j=1}^i((Av_i)^Tv_j)v_j\in\mathcal K_i$. Así, $A\mathcal K_i\subset\mathcal K_i$.
>
> **Paso 2: $AV_i=V_iH$ con $H$ no singular.** Sea $V_i=\begin{pmatrix}v_1&\cdots&v_i\end{pmatrix}$ (columnas ortonormales, base de $\mathcal K_i$). Por el Paso 1, cada $Av_j$ es combinación lineal de $v_1,\dots,v_i$, así que existe una matriz $H\in\mathbb R^{i\times i}$ (con $H_{jk}=(Av_k)^Tv_j$, igual que en el proceso de Arnoldi) tal que
> $$AV_i=V_iH.$$
> $H$ es no singular: si $Hy=0$ para algún $y\neq0$, entonces $AV_iy=V_iHy=0$; pero $V_iy\neq0$ (porque $V_i$ tiene columnas ortonormales e $y\neq0$) y $A$ es invertible, así que $AV_iy\neq0$ - contradicción. Luego $Hy=0\implies y=0$.
>
> **Paso 3: el residuo mínimo sobre $x_0+\mathcal K_i$ es cero.** Sea $\beta=\|r_0\|_2$ y $e_1=(1,0,\dots,0)^T\in\mathbb R^i$; como $v_1=r_0/\beta$, se tiene $r_0=\beta V_ie_1$. Para cualquier $x\in x_0+\mathcal K_i$, escribiendo $x-x_0=V_iy$ con $y\in\mathbb R^i$ (posible porque las columnas de $V_i$ son base de $\mathcal K_i$):
> $$\|b-Ax\|_2 = \|r_0-A(x-x_0)\|_2 = \|\beta V_ie_1-AV_iy\|_2 = \|\beta V_ie_1-V_iHy\|_2 = \|V_i(\beta e_1-Hy)\|_2 = \|\beta e_1-Hy\|_2,$$
> donde la última igualdad usa que $V_i$ tiene columnas ortonormales (preserva norma 2). Como $H$ es no singular (Paso 2), la elección $y=\beta H^{-1}e_1$ anula exactamente esta expresión: $\|\beta e_1-Hy\|_2=0$. Y $0$ es, trivialmente, el mínimo posible de una norma - así que $x_i:=x_0+V_iy$ es el minimizador de $\|b-Ax\|_2$ sobre $x_0+\mathcal K_i$ (la definición de GMRes) con residuo exactamente cero, es decir $x_i=A^{-1}b$.

> [!idea]
> El nombre "*lucky breakdown*" es literal: en GC o en Arnoldi, que una recurrencia se anule antes de tiempo suena a falla del algoritmo - pero acá es la mejor noticia posible, porque significa que la base de Krylov ya alcanzó para representar exactamente la solución. Es el mismo fenómeno que $r_k=0$ en GC ([[GC como metodo directo]]): el espacio de Krylov deja de crecer porque ya "encontró" todo lo que necesitaba.

## Ver también
- [[GMRes - implementacion y rotaciones de Givens]] (nota de tema completa)
- [[GMRes - definicion y propiedad de minimizacion]]
- [[GC como metodo directo]] (el análogo en GC)
- [[Demostraciones - Capitulo 5|Demostraciones]]
