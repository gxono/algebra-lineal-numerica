---
tags: [capitulo4, moc]
---

# Formulario - Capítulo 4 (Iterativos estacionarios)

## Marco general

$A=L+D+U$, $A=N-P$, iteración $x_{k+1}=Mx_k+c$ con $M=N^{-1}P$, $c=N^{-1}b$. Converge $\forall x_0 \iff \rho(M)<1$.

## Los 4 métodos

| Método | $M$ | Condición suficiente |
|---|---|---|
| [[Metodo de Jacobi\|Jacobi]] | $-D^{-1}(L+U)$ | $A$ e.d.d. |
| [[Metodo de Gauss-Seidel\|Gauss-Seidel]] | $-(L+D)^{-1}U$ | $A$ e.d.d., o $A$ sdp |
| [[Metodo SOR\|SOR]] | $\big(\frac{D}{\omega}+L\big)^{-1}\big(\frac{1-\omega}{\omega}D-U\big)$ | $A$ sdp, $\omega\in(0,2)$ (Ostrowski-Reich) |
| [[Metodo de Richardson\|Richardson]] | $I-\alpha A$ | $A$ sdp, $0<\alpha<2/\lambda_{\max}$ |

$$
\alpha_{\text{opt}} = \frac2{\lambda_{\min}+\lambda_{\max}}, \qquad \rho(M_{\alpha_{\text{opt}}}) = \frac{\kappa_2(A)-1}{\kappa_2(A)+1}.
$$
Condición **necesaria** para SOR (cualquier $A$): $\omega\in(0,2)$.

## Normas y espectro

$$
\rho(A) = \max_i|\lambda_i| \leq \|A\| \ \text{(cualquier norma inducida)}, \qquad \rho(A)=\lim_k\|A^k\|^{1/k}.
$$
$e_{\mathrm{d.d.}}$ ([[2.11 - Matrices e.d.d. son no singulares|estrictamente diagonal dominante]]): $|a_{ii}|>\sum_{j\neq i}|a_{ij}| \ \forall i$.

## Número de condición

$$
\kappa(A)=\|A\|\|A^{-1}\|\geq1, \qquad \frac1{\kappa(A)}\frac{\|\Delta b\|}{\|b\|} \leq \frac{\|\Delta x\|}{\|x\|} \leq \kappa(A)\frac{\|\Delta b\|}{\|b\|}.
$$
$A$ simétrica: $\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$. $A$ rectangular: $\kappa_2(A)=\sigma_1/\sigma_m$ (Capítulo 3).

## Ver también
- [[Formulario - Indice]]
- [[Capitulo 4 - MOC]]
