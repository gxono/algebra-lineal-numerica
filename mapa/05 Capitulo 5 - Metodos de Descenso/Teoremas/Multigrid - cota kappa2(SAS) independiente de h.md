---
tags: [capitulo5, teorema]
---

# Multigrid: $\kappa_2(SA_hS)\leq(1+\rho)/(1-\rho)$ (Paso 2)

**Paso 2 ($M$ es un buen precondicionador: de $\sigma(MA_h)$ a $\kappa_2(SA_hS)$).** Sea $S:=M^{1/2}$.

> [!demostracion]- $\kappa_2(SA_hS)\leq\dfrac{1+\rho}{1-\rho}$, independiente de $h$
> Por el [[5.9 - Espectro de S2A y SAS|Problema 5.9]], $SA_hS$ y $MA_h=S^2A_h$ son semejantes (mismo espectro): $SA_hS=S^{-1}(S^2A_h)S$.
>
> Por otro lado, las constantes $c_1,c_2$ de la definición de precondicionador ($c_1\leq\frac{x^TMx}{x^TA_h^{-1}x}\leq c_2$) son también autovalores de $MA_h$: con $x=A_h^{1/2}z$,
> $$\frac{x^TMx}{x^TA_h^{-1}x} = \frac{z^TA_h^{1/2}MA_h^{1/2}z}{z^Tz},$$
> y este cociente de Rayleigh recorre exactamente $[\lambda_{\min}(B),\lambda_{\max}(B)]$ con $B:=A_h^{1/2}MA_h^{1/2}$. Como $B=A_h^{1/2}(MA_h)A_h^{-1/2}$, $B$ es semejante a $MA_h$, así que $\sigma(B)=\sigma(MA_h)$ y por lo tanto
> $$c_1=\lambda_{\min}(MA_h), \qquad c_2=\lambda_{\max}(MA_h), \qquad \kappa_2(SA_hS)=\frac{c_2}{c_1}.$$
> Por el Paso 1, $\sigma(MA_h)\subset[1-\rho,1+\rho]$ con $\rho<1$, así que $c_1\geq1-\rho>0$ y $c_2\leq1+\rho$, de donde
> $$\kappa_2(SA_hS) = \frac{c_2}{c_1} \leq \frac{1+\rho}{1-\rho}. \tag{5.27}$$
> El lado derecho depende solo de $\rho$ (independiente de $h$, por el Teorema 5.34), nunca de $n$.

## Ver también
- [[Multigrid en 1D - transferencia e implementacion]] (nota de tema completa)
- [[Multigrid - autovalores de E en -rho,rho|Multigrid: autovalores de E (Paso 1)]]
- [[5.9 - Espectro de S2A y SAS|Problema 5.9]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
