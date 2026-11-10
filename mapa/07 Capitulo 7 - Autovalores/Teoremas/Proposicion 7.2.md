---
tags: [capitulo7, teorema]
---

# Proposición 7.2 (convergencia del método de las potencias)

Sea $A\in\mathbb{C}^{n\times n}$ con $|\lambda_1|>|\lambda_2|\geq|\lambda_3|\geq\cdots\geq|\lambda_n|$, $x_1$ autovector a derecha de $\lambda_1$ con $\|x_1\|_2=1$, e $y_1^H$ autovector a izquierda de $\lambda_1$. Si $y_1^Hu_0\neq0$, la sucesión $u_k=\sigma_{k-1}Au_{k-1}$, $\sigma_{k-1}=1/\|Au_{k-1}\|$, converge en dirección a $x_1$: existen $\omega_k\in\mathbb{C}$, $|\omega_k|=1$, con $\omega_ku_k\to x_1$. La tasa es asintóticamente $(|\lambda_2/\lambda_1|+\varepsilon)^k$ para cualquier $\varepsilon>0$.

> [!demostracion]-
> **Descomposición.** $\lambda_1$ es dominante, luego simple, y su bloque de Jordan es $[\lambda_1]$. Existe $X=\begin{bmatrix}x_1&X_2\end{bmatrix}$ no singular, $\|x_1\|_2=1$, con
> $$X^{-1}AX=\begin{pmatrix}\lambda_1&0^T\\0&B\end{pmatrix},\qquad X^{-1}=\begin{pmatrix}y_1^H\\Y_2^H\end{pmatrix},$$
> donde $x_1$ es autovector a derecha y $y_1^H$ (primera fila de $X^{-1}$) autovector a izquierda de $\lambda_1$. Además $\sigma(B)=\{\lambda_2,\dots,\lambda_n\}$.
>
> **Potencias.** Escribimos $X^{-1}u_0=\begin{pmatrix}\nu_0\\w_0\end{pmatrix}$, $\nu_0=y_1^Hu_0$. Como $A^k=X\begin{pmatrix}\lambda_1^k&0^T\\0&B^k\end{pmatrix}X^{-1}$:
> $$A^ku_0=X\begin{pmatrix}\lambda_1^k\nu_0\\B^kw_0\end{pmatrix}=\lambda_1^k\big(\nu_0x_1+z_k\big),\qquad z_k:=X_2(B/\lambda_1)^kw_0. \tag{7.2}$$
> Como $\rho(B/\lambda_1)=|\lambda_2|/|\lambda_1|<1$, $(B/\lambda_1)^k\to0$ ([[Teorema 4.23]]) y $z_k\to0$.
>
> **Buena definición.** Si $\nu_0\neq0$, la primera componente de $X^{-1}A^ku_0$ es $\lambda_1^k\nu_0\neq0$, así que $A^ku_0\neq0$ y vale el [[Lema 7.1]].
>
> **Convergencia en dirección.** Normalizando (7.2) con $\lambda_1^k=|\lambda_1|^k(\lambda_1/|\lambda_1|)^k$:
> $$u_k=\frac{A^ku_0}{\|A^ku_0\|}=\Big(\frac{\lambda_1}{|\lambda_1|}\Big)^k\frac{\nu_0x_1+z_k}{\|\nu_0x_1+z_k\|}.$$
> Como $z_k\to0$ y $\|x_1\|_2=1$: $\dfrac{\nu_0x_1+z_k}{\|\nu_0x_1+z_k\|}\to\dfrac{\nu_0}{|\nu_0|}x_1$. Con $\omega_k:=\Big(\dfrac{|\lambda_1|}{\lambda_1}\Big)^k\dfrac{|\nu_0|}{\nu_0}$ (módulo 1), $\omega_ku_k\to x_1$.
>
> **Tasa.** Dado $\varepsilon>0$, existe $C$ (que depende de $B$ y $\varepsilon$) con $\|(B/\lambda_1)^k\|\leq C(\rho(B/\lambda_1)+\varepsilon)^k=C(|\lambda_2/\lambda_1|+\varepsilon)^k$ ([[Teorema 4.24|fórmula de Gelfand]]), y $\|z_k\|$ hereda esa cota.

## Ver también
- [[Metodo de las potencias]]
- [[Proposicion 7.7]] (la misma prueba aplicada a $(A-\mu I)^{-1}$)
- [[Demostraciones - Capitulo 7|Demostraciones]]
