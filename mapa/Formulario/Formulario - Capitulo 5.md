---
tags: [capitulo5, moc]
---

# Formulario - Capítulo 5 (Métodos de descenso)

## El funcional de energía ($A$ sdp)

$$
\phi(x) = \tfrac12x^TAx-b^Tx, \qquad \phi(y) = \phi(x_*)+\tfrac12\|y-x_*\|_A^2, \qquad \|z\|_A:=\sqrt{z^TAz}.
$$
Minimizar $\phi$ sobre $S$ $\iff$ minimizar $\|\cdot-x_*\|_A$ sobre $S$.

## [[Descenso mas pronunciado|Descenso más pronunciado]] (SD)

$$
x_{k+1}=x_k+\alpha r_k, \qquad \alpha=\frac{r_k^Tr_k}{r_k^TAr_k}, \qquad \|x_{k+1}-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A.
$$

## Gradiente conjugado (GC)

$$
d_{k+1}=r_k+\beta_{k+1}d_k, \qquad \beta_{k+1}=\frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2}, \qquad \alpha_{k+1}=\frac{\|r_k\|_2^2}{d_{k+1}^TAd_{k+1}}, \qquad x_{k+1}=x_k+\alpha_{k+1}d_{k+1}.
$$
Costo/iteración: 1 mat-vec, 2 prod. escalares, 3 $axpy$. Termina en $\leq n$ pasos ($\leq k$ si $A$ tiene $k$ autovalores distintos).

$$
\|x_k-x_*\|_A \leq 2\Big(\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\Big)^k\|x_0-x_*\|_A \qquad \text{(cota de Chebyshev)}.
$$
Cota espectral general: $\|x_k-x_*\|_A\leq\big(\max_{\lambda\in\sigma(A)}|p(\lambda)|\big)\|x_0-x_*\|_A$ para cualquier $p\in\overline{\mathcal P}_k$ ($p(0)=1$).

## GC precondicionado (GCP)

Con $M\approx A^{-1}$, $z_k=Mr_k$:
$$
d_{k+1}=z_k+\beta_{k+1}d_k, \qquad \beta_{k+1}=\frac{r_k^Tz_k}{r_{k-1}^Tz_{k-1}}, \qquad \alpha_{k+1}=\frac{r_k^Tz_k}{d_{k+1}^TAd_{k+1}}.
$$
$c_1\leq\dfrac{x^TMx}{x^TA^{-1}x}\leq c_2$, tasa $\Big(\dfrac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\Big)^k$. Precondicionadores: $M_J=D^{-1}$; $M_{SGS}=(L^T+D)^{-1}D(L+D)^{-1}$; [[Precondicionadores genericos|Cholesky incompleta]] $M=(R^TR)^{-1}$.

## Multigrid en 1D

$$
A_h=\frac1{h^2}\text{tridiag}(-1,2,-1), \qquad \omega_{\text{óptimo}}=\tfrac23, \qquad \mu_{\text{suav}}(2/3)=\tfrac13.
$$
Ciclo en V: $O(n)$ por ciclo, $\rho<1$ independiente de $h$. Como precondicionador: $\kappa_2(SA_hS)\leq\dfrac{1+\rho}{1-\rho}$ (independiente de $h$).

## No simétricas: CGNR / CGNE

$$
\text{CGNR: } A^TAx=A^Tb \ (\text{minimiza } \|Ax-b\|_2), \qquad \text{CGNE: } AA^Ty=b,\ x=A^Ty \ (\text{minimiza } \|x-x_*\|_2).
$$
$\kappa_2(A^TA)=\kappa_2(A)^2$.

## [[GMRes - definicion y propiedad de minimizacion|GMRes]]

$$
x_k=\operatorname*{argmin}_{x\in x_0+\mathcal K_k}\|b-Ax\|_2, \qquad \|r_k\|_2=\min_{p\in\overline{\mathcal P}_k}\|p(A)r_0\|_2.
$$
Arnoldi: $AV_k=V_{k+1}H_k$ ($H_k$ Hessenberg). [[Rotaciones de Givens|Givens]] incremental: $\nu=\text{hypot}(h_{kk},h_{k+1,k})$, $c_k=h_{kk}/\nu$, $s_k=-h_{k+1,k}/\nu$. Termina en $\leq n$ pasos.

Precondicionado: izquierda $MAx=Mb$ ($\kappa_2(MA)\rho_I^k$); derecha $(AM)y=b$, $x=My$ ($\kappa_2(A)\rho_D^k$, pero minimiza el residuo verdadero).

## Costo por iteración: GC vs. GMRes

| | GC | GMRes (iter. $k$) |
|---|---|---|
| Mat-vec | 1 | 1 |
| Memoria | $O(n)$, fija | $O(nk)$, crece |
| Costo acumulado hasta $k$ | $O(nk)$ | $O(nk^2)$ |

## Ver también
- [[Formulario - Indice]]
- [[Capitulo 5 - MOC]]
