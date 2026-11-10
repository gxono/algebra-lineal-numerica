---
tags: [capitulo7, definicion, observacion]
---

# Cociente de Rayleigh y criterio de parada (Obs. 7.4–7.6)

> [!idea]
> Dado un autovector aproximado $x$, el mejor autovalor aproximado (en mínimos cuadrados) es el **cociente de Rayleigh** $\mu=\frac{x^HAx}{x^Hx}$. Su error está controlado por el del autovector: $|\lambda-\mu|\leq2\|A\|_2\|x-u\|_2$, y si $A$ es hermítica es **cuadrático** en ese error. El residuo $r=Au-\mu u$ sale gratis y da el criterio de parada: para $A$ hermítica, **hay un autovalor exacto a distancia $\leq\|r\|_2$ de $\mu$**.

## De dónde sale

Dado $x\approx$ autovector, buscar $\mu\in\mathbb{C}$ que minimice $\|Ax-\mu x\|_2$ es un problema de [[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]] con "matriz" $x$ ($n\times1$) y dato $Ax$. Ecuaciones normales: $(x^Hx)\mu=x^H(Ax)$, o sea
$$
\mu=\frac{x^HAx}{x^Hx}\qquad\text{(cociente de Rayleigh)}.
$$

## Observación 7.4 (error del autovalor vs. error del autovector)

$(\lambda,x)$ autopar con $\|x\|_2=1$; $u$ con $\|u\|_2=1$ y $\tilde\lambda=u^HAu$. Entonces
$$
|\lambda-\tilde\lambda|\leq\big(\|Au\|_2+\|Ax\|_2\big)\|x-u\|_2\leq2\|A\|_2\|x-u\|_2.
$$
Como $u_k^HAu_k$ no cambia si se multiplica $u_k$ por $\omega_k$ con $|\omega_k|=1$, en el método de las potencias
$$
|\lambda_1-\mu_k|\leq2\|A\|_2\,\|x_1-\omega_ku_k\|_2\to0.
$$
Ver [[Observacion 7.4]]. La misma nota prueba lo que **no** está en el apunte pero se ve en todos los experimentos: si $A$ es hermítica,
$$
|\lambda_1-\mu|\leq\max_j|\lambda_1-\lambda_j|\,\sin^2\theta\qquad(\theta=\angle(u,x_1)),
$$
**error cuadrático**. Verificado: perturbando $x_1$ con $t\,x_2$, el error del cociente es $2.0\times10^{-2},\,2.0\times10^{-4},\,2.0\times10^{-6}$ para $t=10^{-1},10^{-2},10^{-3}$ (simétrica: $\approx(\lambda_1-\lambda_2)t^2$), contra $5.7\times10^{-2},\,4.1\times10^{-3},\,4.0\times10^{-4}$ para una no simétrica (lineal).

## Observación 7.5 (criterio de parada)

El residuo $r_k=Au_k-\mu_ku_k=v_{k+1}-\mu_ku_k$ **no cuesta nada** (el producto $Au_k$ ya se hizo).
- Siempre: $(A-r_ku_k^H)u_k=\mu_ku_k$, o sea $(\mu_k,u_k)$ es autopar **exacto** de una matriz que difiere de $A$ en $\|r_ku_k^H\|_2=\|r_k\|_2$ (**error hacia atrás**).
- Si $A$ es **hermítica** (base ortonormal de autovectores, $u_k=\sum c_jx_j$):
$$
\|r_k\|_2^2=\sum_j|\lambda_j-\mu_k|^2|c_j|^2\geq\min_j|\lambda_j-\mu_k|^2\implies\min_j|\lambda_j-\mu_k|\leq\|r_k\|_2.
$$
Justifica el test **relativo** $\|r_k\|_2\leq\text{tol}\,|\mu_k|$, que no requiere conocer $\|A\|_2$ (para hermítica $|\lambda_1|=\rho(A)=\|A\|_2$ y $\mu_k\to\lambda_1$). Ver [[Observacion 7.5]].

## Observación 7.6 (potencias "al revés")

Si en el algoritmo se reemplaza $v_{k+1}=Au_k$ por **resolver** $Av=u_k$ ($v_{k+1}=A^{-1}u_k$), se aproxima el autovalor dominante de $A^{-1}$, que es $1/\lambda_{\min}$, el **inverso del autovalor de módulo mínimo** de $A$. Es la puerta de entrada al [[Metodo de las potencias inversas|método de las potencias inversas]].

## Ver también
- [[Metodo de las potencias]]
- [[Metodo de las potencias inversas]] (iteración del cociente de Rayleigh)
- [[Glosario - Indice|Glosario]]
- [[Capitulo 7 - MOC]]
