---
tags: [capitulo7, teorema]
---

# Proposición 7.7 (convergencia de las potencias inversas)

Sea $A\in\mathbb{C}^{n\times n}$, $\lambda_1$ autovalor simple con autovectores $x_1$ ($\|x_1\|_2=1$) a derecha e $y_1^H$ a izquierda, y $\mu\notin\sigma(A)$ con $|\lambda_1-\mu|<|\lambda_i-\mu|$ para todo $\lambda_i\neq\lambda_1$. Si $y_1^Hu_0\neq0$, las potencias inversas convergen en dirección a $x_1$ con tasa $(\rho_\mu+\varepsilon)^k$,
$$
\rho_\mu:=\max_{\lambda_i\neq\lambda_1}\left|\frac{\lambda_1-\mu}{\lambda_i-\mu}\right|,
$$
y $\mu_k\to\lambda_1$.

> [!demostracion]-
> **Autopares de $B:=(A-\mu I)^{-1}$.** Si $Ax=\lambda x$ con $\lambda\neq\mu$, $(A-\mu I)x=(\lambda-\mu)x$; aplicando $(A-\mu I)^{-1}$ y dividiendo por $\lambda-\mu$:
> $$Bx=\frac{1}{\lambda-\mu}x.$$
> $B$ tiene autovalores $\theta_i=(\lambda_i-\mu)^{-1}$ con **los mismos autovectores** que $A$ (lo mismo a izquierda).
>
> **$\theta_1$ es dominante.** La hipótesis sobre $\mu$ es exactamente $\left|\frac{1}{\lambda_1-\mu}\right|>\left|\frac{1}{\lambda_i-\mu}\right|$ para $\lambda_i\neq\lambda_1$: $\theta_1$ es el autovalor dominante de $B$, con autovectores $x_1$ e $y_1^H$.
>
> **Aplicar la Proposición 7.2.** Las potencias inversas son el método de las potencias aplicado a $B$. Por la [[Proposicion 7.2]], convergen en dirección a $x_1$ con tasa $(|\theta_2/\theta_1|+\varepsilon)^k$, donde $\theta_2$ es el segundo de mayor módulo:
> $$\left|\frac{\theta_2}{\theta_1}\right|=\max_{\lambda_i\neq\lambda_1}\frac{|\lambda_1-\mu|}{|\lambda_i-\mu|}=\rho_\mu.$$
> **El autovalor.** $u_k^Hv_{k+1}=u_k^HBu_k$ es el cociente de Rayleigh de $B$ en $u_k$, que converge a $\theta_1$ ([[Observacion 7.4]]). Luego
> $$\mu_k=\big(u_k^Hv_{k+1}\big)^{-1}+\mu\to\theta_1^{-1}+\mu=\lambda_1.$$

## Ver también
- [[Metodo de las potencias inversas]]
- [[7.3 - Programar potencias inversas]] (tasa verificada: $\rho_\mu=0.111\to$ medida $\approx0.11$)
- [[Demostraciones - Capitulo 7|Demostraciones]]
