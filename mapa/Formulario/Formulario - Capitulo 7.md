---
tags: [capitulo7, moc]
---

# Formulario - Capítulo 7 (Autovalores)

## Básicos

- $p(\lambda)=\det(\lambda I-A)$; $\sigma(A^T)=\sigma(A)$; $\sigma(A^H)=\overline{\sigma(A)}$; autovector a izquierda: $y^HA=\lambda y^H$.
- Mult. geométrica $=n-\text{rank}(A-\lambda I)\leq$ algebraica; si $<$: **defectuoso** (no diagonalizable).
- Hermítica: $\lambda_i\in\mathbb{R}$, $U^HAU=\text{diag}(\lambda_i)$ con $U$ unitaria.

## [[Metodo de las potencias|Potencias]]

$$
u_k=\frac{Au_{k-1}}{\|Au_{k-1}\|},\quad\mu_k=u_k^HAu_k,\quad r_k=Au_k-\mu_ku_k,\quad\text{parar si }\|r_k\|_2\leq\text{tol}\,|\mu_k|.
$$
Hipótesis: $|\lambda_1|>|\lambda_2|$, $y_1^Hu_0\neq0$. Tasa del autovector $|\lambda_2/\lambda_1|$; del autovalor $|\lambda_2/\lambda_1|$ (general) o $|\lambda_2/\lambda_1|^2$ (hermítica).

## [[Cociente de Rayleigh y criterio de parada|Cociente de Rayleigh]]

$$
\mu=\frac{x^HAx}{x^Hx}=\arg\min_\mu\|Ax-\mu x\|_2;\qquad|\lambda-\mu|\leq2\|A\|_2\|x-u\|_2;\qquad\text{hermítica: }|\lambda_1-\mu|\leq\max_j|\lambda_1-\lambda_j|\sin^2\theta.
$$
Hermítica: $\min_j|\lambda_j-\mu|\leq\|Au-\mu u\|_2$ ($\|u\|=1$).

## [[Metodo de las potencias inversas|Potencias inversas]]

$$
(A-\mu I)v_{k+1}=u_k,\quad u_{k+1}=\frac{v_{k+1}}{\|v_{k+1}\|},\quad\lambda\approx\frac{1}{u_k^Hv_{k+1}}+\mu;\qquad\rho_\mu=\max_{\lambda_i\neq\lambda_*}\left|\frac{\lambda_*-\mu}{\lambda_i-\mu}\right|.
$$
Shift actualizado $\mu_k=u_k^HAu_k$: orden 3 (hermítica), 2 (general).

## [[Algoritmo QR - definicion y convergencia|Algoritmo QR]]

$$
A_k=Q_kR_k,\ A_{k+1}=R_kQ_k=Q_k^HA_kQ_k;\qquad A_{k+1}=P_k^HAP_k,\quad A^{k+1}=P_kU_k.
$$
Converge a Schur (diagonal ordenada por módulo) con tasa $\max_i|\lambda_{i+1}/\lambda_i|$; $(A_k)_{11}\to\lambda_1$ con $|\lambda_2/\lambda_1|$, $(A_k)_{nn}\to\lambda_n$ con $|\lambda_n/\lambda_{n-1}|$. Hermítica: $A_k\to\text{diag}(\lambda_i)$, $P_k\to$ autovectores.

## [[Algoritmo QR - Hessenberg, shifts y deflacion|QR práctico]]

- Hessenberg ($h_{ij}=0$ si $i>j+1$): reducción $O(n^3)$; paso QR $O(n^2)$; se conserva.
- Shift: $A_k-\mu_kI=Q_kR_k$, $A_{k+1}=R_kQ_k+\mu_kI$; $\prod(A-\mu_jI)=P_kU_k$.
- Rayleigh $\mu_k=a^{(k)}_{nn}$; Wilkinson: autovalor de $\begin{pmatrix}a_{n-1,n-1}&a_{n-1,n}\\a_{n,n-1}&a_{nn}\end{pmatrix}$ más cercano a $a_{nn}$.
- Deflación: $|a^{(k)}_{n,n-1}|<\text{tol}$ ⟹ $\lambda_n\approx a_{nn}$, seguir con $(n-1)\times(n-1)$.

## Matriz $A_N=\text{tridiag}(-1,2,-1)$

$$
\lambda_i=4\sin^2\frac{i\pi}{2(N+1)},\quad(x_i)_j=\sin\frac{ij\pi}{N+1};\qquad\frac{\lambda_{N-1}}{\lambda_N}=\frac{\cos^2\frac{\pi}{N+1}}{\cos^2\frac{\pi}{2(N+1)}}\approx1-\frac{3\pi^2}{4(N+1)^2},\qquad\frac{\lambda_1}{\lambda_2}=\frac{1}{4\cos^2\frac{\pi}{2(N+1)}}\to\frac14.
$$
