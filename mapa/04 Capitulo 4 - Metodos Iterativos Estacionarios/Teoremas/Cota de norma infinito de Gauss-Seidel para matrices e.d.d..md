---
tags: [capitulo4, teorema]
---

# Cota de $\|M_{GS}\|_\infty$ para matrices e.d.d. (misma razón que Jacobi)

Sea $A\in\mathbb R^{n\times n}$ e.d.d., con razón de dominancia uniforme
$$
\alpha := \max_{1\leq i\leq n}\frac{\sum_{j\neq i}|a_{ij}|}{|a_{ii}|} < 1
$$
(la misma $\alpha$ que acota $\|M_J\|_\infty$ para [[Metodo de Jacobi|Jacobi]], por definición de la norma $\infty$ inducida - ver [[4.5(b) - Formula norma infinito|Problema 4.5(b)]]). Entonces
$$
\|M_{GS}\|_\infty \leq \alpha.
$$

> [!idea]
> Para Jacobi esto es **inmediato**: $M_J=-D^{-1}(L+U)$, así que la fila $i$ de $M_J$ tiene suma de valores absolutos $\sum_{j\neq i}|a_{ij}|/|a_{ii}|\leq\alpha$ directamente, por la propia definición de e.d.d. - no hace falta demostrar nada más. Para Gauss-Seidel, $M_{GS}=-(L+D)^{-1}U$ involucra una matriz **inversa**, así que sus entradas no se leen directamente de $A$ - por eso hace falta un argumento (la inducción de abajo) en vez de una sola línea. Lo notable es que, pese a esa complicación, **el resultado final es exactamente la misma cota** que para Jacobi.

> [!demostracion]-
> Sea $e\in\mathbb R^n$ con $\|e\|_\infty\leq1$, y sea $\hat e:=M_{GS}e$. De $(L+D)\hat e=-Ue$, la fila $i$ da
> $$a_{ii}\hat e_i = -\sum_{j<i}a_{ij}\hat e_j - \sum_{j>i}a_{ij}e_j.$$
> Demostramos $|\hat e_i|\leq\alpha$ para todo $i$, por inducción en $i$.
>
> **Caso base $i=1$:** no hay términos $j<1$, así que $a_{11}\hat e_1=-\sum_{j>1}a_{1j}e_j$, y usando $|e_j|\leq1$:
> $$|a_{11}||\hat e_1| \leq \sum_{j>1}|a_{1j}| \leq \sum_{j\neq1}|a_{1j}| \leq \alpha|a_{11}| \implies |\hat e_1|\leq\alpha.$$
>
> **Paso inductivo:** supongamos $|\hat e_j|\leq\alpha$ para todo $j<i$. Usando $|\hat e_j|\leq\alpha$ (para $j<i$) y $|e_j|\leq1$ (para $j>i$, componentes de $e$ todavía no actualizadas):
> $$|a_{ii}||\hat e_i| \leq \sum_{j<i}|a_{ij}||\hat e_j| + \sum_{j>i}|a_{ij}||e_j| \leq \alpha\sum_{j<i}|a_{ij}| + \sum_{j>i}|a_{ij}|.$$
> Por e.d.d. con razón uniforme $\alpha$, $\sum_{j\neq i}|a_{ij}|\leq\alpha|a_{ii}|$. Como además $\alpha\sum_{j<i}|a_{ij}|\leq\alpha\sum_{j\neq i}|a_{ij}|\leq\alpha\cdot\alpha|a_{ii}|$ y $\sum_{j>i}|a_{ij}|\leq\sum_{j\neq i}|a_{ij}|\leq\alpha|a_{ii}|$, alcanza con la cota más directa:
> $$\alpha\sum_{j<i}|a_{ij}| + \sum_{j>i}|a_{ij}| \leq \sum_{j<i}|a_{ij}| + \sum_{j>i}|a_{ij}| = \sum_{j\neq i}|a_{ij}| \leq \alpha|a_{ii}|$$
> (en el primer paso se usó $\alpha<1$ para acotar $\alpha\sum_{j<i}|a_{ij}|\leq\sum_{j<i}|a_{ij}|$ - una cota más floja que la que da la hipótesis inductiva, pero que ya alcanza para cerrar el argumento). Por lo tanto $|\hat e_i|\leq\alpha$, completando la inducción.
>
> Así, $\|\hat e\|_\infty=\max_i|\hat e_i|\leq\alpha$ para todo $e$ con $\|e\|_\infty\leq1$, es decir $\|M_{GS}\|_\infty\leq\alpha$.

> [!pendiente]
> Verificado numéricamente con matrices e.d.d. aleatorias (no tridiagonales, filas con razones de dominancia distintas entre sí) de tamaños $n=2,\dots,20$: la cota $\|M_{GS}\|_\infty\leq\alpha$ se cumple siempre, con margen estrictamente positivo salvo casos degenerados - confirma que el argumento de arriba no depende de ninguna estructura particular de $A$ (tridiagonal, simétrica, etc.), solo de e.d.d. con razón uniforme $\alpha$.

## Relación con el teorema de convergencia general

Esto es una **cota cuantitativa** adicional, más fuerte que el [[Teorema de convergencia de Gauss-Seidel con e.d.d.|teorema cualitativo]] ($\rho(M_{GS})<1$, demostrado ahí con un argumento distinto, vía autovectores): acá se obtiene directamente una tasa explícita $\|M_{GS}\|_\infty\leq\alpha$, útil para estimar **cuántas** iteraciones hacen falta (no solo que el método converge). Es exactamente lo que permite calcular el $\mu_n$ de [[Parcial Modelo - Problema 4|Parcial Modelo - Problema 4(a)]].

## Ver también
- [[Metodo de Gauss-Seidel]] (nota de tema completa)
- [[Teorema de convergencia de Gauss-Seidel con e.d.d.]] (resultado cualitativo, demostración distinta)
- [[Metodo de Jacobi]] (cota análoga, trivial por definición)
- [[Parcial Modelo - Problema 4]] (aplicación con $n=1000$)
- [[Demostraciones - Capitulo 4|Demostraciones]]
