---
tags: [capitulo7, lema]
---

# Lema 7.1

Sea $u_0\in\mathbb{C}^n$ con $\|u_0\|=1$ y supongamos $A^ku_0\neq0$ para todo $k$. Entonces la recursión $u_k=\sigma_{k-1}Au_{k-1}$, $\sigma_{k-1}=1/\|Au_{k-1}\|$, está bien definida y
$$
u_k=\frac{A^ku_0}{\|A^ku_0\|},\qquad k=0,1,2,\dots
$$

> [!demostracion]- (inducción)
> $k=0$: es $\|u_0\|=1$. Si vale para $k$,
> $$Au_k=\frac{A^{k+1}u_0}{\|A^ku_0\|}\neq0,$$
> así que $\sigma_k$ está bien definido y, al normalizar, $\|A^ku_0\|$ se cancela:
> $$u_{k+1}=\frac{Au_k}{\|Au_k\|}=\frac{A^{k+1}u_0/\|A^ku_0\|}{\|A^{k+1}u_0\|/\|A^ku_0\|}=\frac{A^{k+1}u_0}{\|A^{k+1}u_0\|}.$$

Consecuencia: el análisis de convergencia se reduce a estudiar $A^ku_0$ y normalizar recién al final.

## Ver también
- [[Metodo de las potencias]]
- [[Proposicion 7.2]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
