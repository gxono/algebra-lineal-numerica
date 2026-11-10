---
tags: [capitulo4, teorema]
---

# Teorema 4.21 (descomposición de Schur)

Sea $A\in\mathbb C^{n\times n}$ y sean $\lambda_1,\dots,\lambda_n$ sus autovalores contados con multiplicidad algebraica (los ceros de $p(\lambda):=\det(\lambda I-A)$). Entonces existe $U$ unitaria ($\bar U^TU=I$) tal que
$$
UA\bar U^T=T,
$$
con $T$ triangular superior y $\text{diag}(T)=(\lambda_1,\dots,\lambda_n)$.

> [!idea]
> Si los autovalores de $A$ son todos reales, $U$ puede tomarse real y ortogonal ($U^TU=I$), con $T$ triangular superior real. Vale **aunque $A$ no sea diagonalizable** - es más general que la diagonalización.

> [!demostracion]- Por inducción en $n$
> **Caso base $n=1$:** $A=[a_{11}]$, $\lambda_1=a_{11}$; el teorema vale con $T=[a_{11}]$, $U=[1]$.
>
> **Paso inductivo:** supongamos el teorema válido para matrices de orden $n-1$, y veámoslo para $A$ de orden $n>1$.
>
> Sean $\lambda_1,\dots,\lambda_n$ los autovalores de $A$ y $x$ un autovector asociado a $\lambda_1$, normalizado ($\|x\|_2^2=\bar x^Tx=1$) con $x_1\geq0$ (si $\lambda_1$ es complejo y $x_1$ no es real no-negativo, se divide $x$ por $x_1/|x_1|$ y se renormaliza - esto no cambia que sea autovector).
>
> Sea $H=I-u\bar u^T$ con $u=\dfrac{1}{\sqrt{1+x_1}}(x+e_1)$ la reflexión de Householder que satisface $\bar H^T=H=H^{-1}$ (es unitaria) y $Hx=-e_1$, $He_1=-x$ (misma construcción que en [[Matrices de Householder]], Capítulo 3, extendida a $\mathbb C^n$).
>
> Entonces $HA\bar H^T=HAH$ (por ser $H$ real y unitaria, $\bar H=H$), y su primera columna es
> $$\text{col}_1(HAH) = HAHe_1 = -HAx = -\lambda_1Hx = \lambda_1e_1.$$
> Es decir,
> $$HAH = \begin{pmatrix}\lambda_1 & v^T\\ 0 & \tilde A\end{pmatrix}, \qquad v\in\mathbb C^{n-1},\ \tilde A\in\mathbb C^{(n-1)\times(n-1)}.$$
> Además, como $\det(\lambda I-HAH)=(\lambda-\lambda_1)\det(\lambda I-\tilde A)$ (determinante de una matriz triangular por bloques), $\lambda\in\sigma(A)\iff\lambda\in\sigma(HAH)\iff\lambda=\lambda_1$ o $\lambda\in\sigma(\tilde A)$. Por lo tanto $\sigma(\tilde A)=\{\lambda_2,\dots,\lambda_n\}$.
>
> Por hipótesis inductiva, existe $\tilde U$ unitaria de orden $n-1$ tal que $\tilde U\tilde A\bar{\tilde U}^T=\tilde T$, con $\tilde T$ triangular superior y $\text{diag}(\tilde T)=(\lambda_2,\dots,\lambda_n)$, es decir $\tilde A=\bar{\tilde U}^T\tilde T\tilde U$. Sustituyendo:
> $$HAH = \begin{pmatrix}\lambda_1 & v^T\\ 0 & \bar{\tilde U}^T\tilde T\tilde U\end{pmatrix} = \begin{pmatrix}1&0^T\\0&\bar{\tilde U}^T\end{pmatrix}\begin{pmatrix}\lambda_1&v^T\bar{\tilde U}^T\\0&\tilde T\end{pmatrix}\begin{pmatrix}1&0^T\\0&\tilde U\end{pmatrix}.$$
> Multiplicando a izquierda por $U:=\begin{pmatrix}1&0^T\\0&\tilde U\end{pmatrix}$ y a derecha por $\bar U^T=\begin{pmatrix}1&0^T\\0&\bar{\tilde U}^T\end{pmatrix}$:
> $$U(HAH)\bar U^T = \begin{pmatrix}\lambda_1&v^T\bar{\tilde U}^T\\0&\tilde T\end{pmatrix} =: T,$$
> triangular superior con $\text{diag}(T)=(\lambda_1,\lambda_2,\dots,\lambda_n)$. Como $U$ y $H$ son ambas unitarias, $UH$ también lo es, y $T=(UH)A\overline{(UH)}^T$ - exactamente lo que había que probar, con la matriz unitaria $UH$ en lugar de $U$.

## Ver también
- [[Teorema 4.20]] (usa este teorema para caracterizar $\rho(A)$)
- [[Matrices de Householder]] (Capítulo 3)
- [[Demostraciones - Capitulo 4|Demostraciones]]
