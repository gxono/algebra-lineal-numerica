---
tags: [capitulo7, teorema]
---

# Autovalores y autovectores de $A_N=\text{tridiag}(-1,2,-1)$

$(x_i)_j=\sin\frac{\pi ij}{N+1}$ es autovector de $A_N$ con autovalor $\lambda_i=4\sin^2\frac{\pi i}{2(N+1)}$, $i=1,\dots,N$.

![[7.8 - QR en la matriz A_N#(a)]]

**Fórmulas útiles derivadas** (todas usadas en los problemas):
- $\kappa_2(A_N)=\dfrac{\lambda_N}{\lambda_1}=\cot^2\dfrac{\pi}{2(N+1)}\sim\dfrac{4(N+1)^2}{\pi^2}$.
- Potencias: $\dfrac{\lambda_{N-1}}{\lambda_N}=\dfrac{\cos^2\frac{\pi}{N+1}}{\cos^2\frac{\pi}{2(N+1)}}\approx1-\dfrac{3\pi^2}{4(N+1)^2}$.
- Potencias inversas ($\mu=0$): $\dfrac{\lambda_1}{\lambda_2}=\dfrac{1}{4\cos^2\frac{\pi}{2(N+1)}}\to\dfrac14$.
- Simetría del espectro: $\lambda_i+\lambda_{N+1-i}=4$.

## Ver también
- [[7.8 - QR en la matriz A_N]], [[7.2 - Potencias en la matriz A_N]], [[7.4 - Potencias inversas en la matriz A_N]]
- [[Multigrid en 1D - problema modelo y suavizado]] (la misma matriz, Cap. 5)
- [[Demostraciones - Capitulo 7|Demostraciones]]
