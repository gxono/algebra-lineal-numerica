---
tags: [capitulo7, observacion]
---

# Observación 7.4 (error del cociente de Rayleigh)

Sea $(\lambda,x)$ autopar de $A$ con $\|x\|_2=1$, y $u$ con $\|u\|_2=1$, $\tilde\lambda=u^HAu$. Entonces
$$
|\lambda-\tilde\lambda|\leq2\|A\|_2\|x-u\|_2.
$$
**Además (no está en el apunte):** si $A$ es hermítica y $\lambda=\lambda_1$,
$$
|\lambda_1-\tilde\lambda|\leq\max_j|\lambda_1-\lambda_j|\,\sin^2\theta,\qquad\theta=\angle(u,x_1),
$$
es decir, el error es **cuadrático** en el error del autovector.

> [!demostracion]-
> **Cota general.** Sumando y restando $x^HAu$:
> $$\lambda-\tilde\lambda=x^HAx-u^HAu=x^HA(x-u)+(x-u)^HAu,$$
> y por Cauchy-Schwarz, $|\lambda-\tilde\lambda|\leq(\|A^Hx\|_2+\|Au\|_2)\|x-u\|_2\leq2\|A\|_2\|x-u\|_2$. *(El apunte escribe $\|Ax\|_2$; con $\|x\|=\|u\|=1$ ambos términos están acotados por $\|A\|_2$.)*
>
> **En el método de las potencias.** $(\omega u_k)^HA(\omega u_k)=\bar\omega\omega u_k^HAu_k=u_k^HAu_k$ si $|\omega|=1$. Aplicando la cota a $\omega_ku_k$:
> $$|\lambda_1-\mu_k|=|x_1^HAx_1-(\omega_ku_k)^HA(\omega_ku_k)|\leq2\|A\|_2\|x_1-\omega_ku_k\|_2\to0.$$
>
> **Caso hermítico (cuadrático).** Sea $\{x_j\}$ base ortonormal de autovectores ([[Autovalores y diagonalizacion de matrices hermiticas|teorema espectral]]) y $u=c_1x_1+w$ con $w\perp x_1$, $|c_1|^2+\|w\|^2=1$, $\|w\|=\sin\theta$. Como $Aw\perp x_1$ (porque $A$ deja invariante $x_1^\perp$):
> $$\tilde\lambda=u^HAu=\lambda_1|c_1|^2+w^HAw\implies\lambda_1-\tilde\lambda=\lambda_1\|w\|^2-w^HAw=w^H(\lambda_1I-A)w,$$
> y $|w^H(\lambda_1I-A)w|\leq\max_{j\geq2}|\lambda_1-\lambda_j|\,\|w\|^2$.

**Verificación numérica** ($\lambda_1=10$, $\lambda_2=8$, $u\propto x_1+t\,x_2$, o sea $\sin^2\theta\approx t^2$):

| $t$ | simétrica: $\vert \lambda_1-\tilde\lambda\vert $ | no simétrica |
|---|---|---|
| $10^{-1}$ | $1.98\times10^{-2}$ | $5.66\times10^{-2}$ |
| $10^{-2}$ | $2.00\times10^{-4}$ | $4.12\times10^{-3}$ |
| $10^{-3}$ | $2.00\times10^{-6}$ | $3.96\times10^{-4}$ |

Cuadrático ($\approx2t^2=(\lambda_1-\lambda_2)t^2$) contra lineal. Por eso en simétricas el autovalor converge a tasa $|\lambda_2/\lambda_1|^2$ en el método de las potencias.

## Ver también
- [[Cociente de Rayleigh y criterio de parada]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
