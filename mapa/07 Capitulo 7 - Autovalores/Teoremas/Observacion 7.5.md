---
tags: [capitulo7, observacion]
---

# Observación 7.5 (el residuo como criterio de parada)

Sea $u_k$ con $\|u_k\|_2=1$, $\mu_k=u_k^HAu_k$ y $r_k=Au_k-\mu_ku_k$. Entonces:
1. $(\mu_k,u_k)$ es un autopar **exacto** de $A-r_ku_k^H$, y $\|r_ku_k^H\|_2=\|r_k\|_2$.
2. Si $A$ es hermítica, existe $\lambda_j\in\sigma(A)$ con $|\lambda_j-\mu_k|\leq\|r_k\|_2$.

> [!demostracion]-
> **1.** $(A-r_ku_k^H)u_k=Au_k-r_k(u_k^Hu_k)=Au_k-r_k=\mu_ku_k$. Para una matriz de rango 1, $\|ab^H\|_2=\|a\|_2\|b\|_2$, así que $\|r_ku_k^H\|_2=\|r_k\|_2$.
>
> **2.** Con una base ortonormal de autovectores $\{x_j\}$, $u_k=\sum_jc_jx_j$, $\sum_j|c_j|^2=1$:
> $$r_k=\sum_j(\lambda_j-\mu_k)c_jx_j,\qquad\|r_k\|_2^2=\sum_j|\lambda_j-\mu_k|^2|c_j|^2\geq\min_j|\lambda_j-\mu_k|^2\sum_j|c_j|^2=\min_j|\lambda_j-\mu_k|^2.$$

Esto justifica el test relativo $\|r_k\|_2\leq\text{tol}\,|\mu_k|$: por la [[Observacion 7.4]], $\mu_k\to\lambda_1$, y para hermíticas $|\lambda_1|=\rho(A)=\|A\|_2$, así que el test equivale a $\|r_k\|\lesssim\text{tol}\,\|A\|_2$ sin calcular $\|A\|_2$.

## Ver también
- [[Cociente de Rayleigh y criterio de parada]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
