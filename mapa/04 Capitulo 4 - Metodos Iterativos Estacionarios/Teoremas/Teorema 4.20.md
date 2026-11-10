---
tags: [capitulo4, teorema]
---

# Teorema 4.20 (caracterización del radio espectral)

Sea $A\in\mathbb R^{n\times n}$. Entonces
$$
\rho(A) = \inf\{\|A\| : \|\cdot\|\text{ es una norma inducida}\}.
$$

> [!demostracion]-
> **La parte fácil ($\rho(A)\leq\|A\|$ para toda norma inducida):** si $\lambda\in\sigma(A)$ con autovector $x\neq0$, $Ax=\lambda x$, así que $|\lambda|\|x\|=\|\lambda x\|=\|Ax\|\leq\|A\|\|x\|$, y como $\|x\|>0$, $|\lambda|\leq\|A\|$. Como $\lambda$ era arbitrario, $\rho(A)\leq\|A\|$. Esto ya da $\rho(A)\leq\inf\{\|A\|\}$.
>
> **La parte difícil (dado $\varepsilon>0$, existe una norma inducida con $\|A\|<\rho(A)+\varepsilon$):** por el [[Teorema 4.21]] (descomposición de Schur) existe $U$ unitaria tal que $UA\bar U^T=T$, con $T$ triangular superior y $\text{diag}(T)=(\lambda_1,\dots,\lambda_n)$, los autovalores de $A$:
> $$T=\begin{pmatrix}\lambda_1&t_{12}&\cdots&t_{1n}\\0&\lambda_2&\cdots&t_{2n}\\\vdots&&\ddots&\vdots\\0&\cdots&0&\lambda_n\end{pmatrix}.$$
> Para $z>0$, sea $D_z=\text{diag}(z,z^2,\dots,z^n)$. Entonces $D_zT D_z^{-1}$ tiene las mismas entradas diagonales que $T$ (los $\lambda_i$) y, fuera de la diagonal, en la posición $(i,j)$ con $i<j$, la entrada $t_{ij}$ queda multiplicada por $z^{i-j}=z^{-(j-i)}$:
> $$D_zT D_z^{-1} = \begin{pmatrix}\lambda_1&z^{-1}t_{12}&z^{-2}t_{13}&\cdots\\0&\lambda_2&z^{-1}t_{23}&\cdots\\ \vdots&&\ddots&\\0&\cdots&0&\lambda_n\end{pmatrix}.$$
> Como $z^{-(j-i)}\to0$ cuando $z\to\infty$ (para $j>i$), eligiendo $z$ suficientemente grande se logra
> $$\max_{1\leq j\leq n}\sum_{i=1}^{j-1}z^{i-j}|t_{ij}| < \varepsilon,$$
> y entonces, usando que $\|B\|_1=\max_j\sum_i|b_{ij}|$ (máxima suma de valores absolutos por columna, [[4.5 - Formulas normas matriciales|Problema 4.5]]):
> $$\big\|D_zTD_z^{-1}\big\|_1 = \max_{1\leq j\leq n}\Big(|\lambda_j|+\sum_{i=1}^{j-1}z^{i-j}|t_{ij}|\Big) < \rho(A)+\varepsilon.$$
> Falta ver que $\|D_zTD_z^{-1}\|_1$ **es** la norma inducida de $A$ por alguna norma vectorial concreta. Definimos $\|x\|:=\|D_zUx\|_1$ para $x\in\mathbb R^n$ (es norma: composición de la norma $\|\cdot\|_1$ con la aplicación lineal invertible $D_zU$). Para esta norma, la norma inducida de $A$ es
> $$\|A\| = \max_{x\neq0}\frac{\|Ax\|}{\|x\|} = \max_{x\neq0}\frac{\|D_zUAx\|_1}{\|D_zUx\|_1} \overset{y:=D_zUx}{=} \max_{y\neq0}\frac{\|D_zUA(D_zU)^{-1}y\|_1}{\|y\|_1} = \big\|D_zUA\bar U^TD_z^{-1}\big\|_1 = \big\|D_zTD_z^{-1}\big\|_1,$$
> donde el cambio de variable $y=D_zUx$ es válido porque $D_zU$ es invertible (biyección $\mathbb R^n\to\mathbb R^n$ cuando $x$ recorre $\mathbb R^n\setminus\{0\}$). Así, $\|A\|=\|D_zTD_z^{-1}\|_1<\rho(A)+\varepsilon$, que es exactamente lo que había que mostrar.

> [!idea]
> La parte difícil explica por qué, aunque $\rho(A)\leq\|A\|$ siempre, no hay una brecha "estructural": para **cualquier** matriz, por más que $\rho(A)$ sea mucho menor que alguna norma particular, siempre hay *otra* norma inducida que se acerca arbitrariamente a $\rho(A)$. Esto es justo lo que permite, en el [[Teorema 4.23]], pasar de "$\rho(A)<1$" a "existe una norma inducida con $\|A\|<1$" - el paso (1)$\Rightarrow$(2) de esa cadena de equivalencias.

## Ver también
- [[Teorema 4.21]] (descomposición de Schur, usada acá)
- [[Teorema 4.23]] (usa este teorema)
- [[Espectro y radio espectral]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
