---
tags: [capitulo5, teorema]
---

# Teorema 5.19 (cota general vía Chebyshev)

$$
\|x_k-x_*\|_A \leq 2\Big(\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\Big)^k\|x_0-x_*\|_A.
$$

> [!demostracion]-
> Sean $0<\lambda_{\min}=\lambda_1\leq\cdots\leq\lambda_n=\lambda_{\max}$ los autovalores de $A$, $c:=\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$. Recordemos los polinomios de Chebyshev $T_0=1$, $T_1=x$, $T_k=2xT_{k-1}-T_{k-2}$, que satisfacen $T_k(x)=\cos(k\arccos x)$ para $x\in[-1,1]$ y $\max_{x\in[-1,1]}|T_k(x)|\leq1$.
>
> Con la biyección $x(\lambda)=\dfrac{\lambda-\frac{\lambda_{\max}+\lambda_{\min}}2}{\frac{\lambda_{\max}-\lambda_{\min}}2}:[\lambda_{\min},\lambda_{\max}]\to[-1,1]$, definimos
> $$p_k(\lambda) := \frac{T_k(x(\lambda))}{T_k(x(0))} \in\overline{\mathcal P}_k$$
> (grado $k$ y $p_k(0)=1$ por construcción, siempre que $T_k(x(0))\neq0$, lo que se verifica abajo).
>
> **Cota del numerador.** Como $\lambda\in\sigma(A)\subset[\lambda_{\min},\lambda_{\max}]$ y $\max_{[-1,1]}|T_k|\leq1$:
> $$\max_{\lambda\in\sigma(A)}|p_k(\lambda)| \leq \frac{1}{|T_k(x(0))|}. \tag{5.7}$$
>
> **Cota del denominador.** $x(0)=-\dfrac{\lambda_{\max}+\lambda_{\min}}{\lambda_{\max}-\lambda_{\min}}=-\dfrac{c+1}{c-1}=\dfrac{z+z^{-1}}2$, con $z:=-\dfrac{\sqrt c-1}{\sqrt c+1}$ (verificar algebraicamente). Usando la identidad $T_k\big(\frac{z+z^{-1}}2\big)=\frac{z^k+z^{-k}}2$ (válida para cualquier $z\neq0$, real o complejo - en particular fuera de $[-1,1]$; demostración por inducción abajo):
> $$T_k(x(0)) = \frac12\Big[\Big(-\frac{\sqrt c-1}{\sqrt c+1}\Big)^k+\Big(-\frac{\sqrt c+1}{\sqrt c-1}\Big)^k\Big] = (-1)^k\frac12\Big[\Big(\frac{\sqrt c-1}{\sqrt c+1}\Big)^k+\Big(\frac{\sqrt c+1}{\sqrt c-1}\Big)^k\Big],$$
> y como ambos términos entre corchetes son positivos,
> $$|T_k(x(0))| \geq \frac12\Big(\frac{\sqrt c+1}{\sqrt c-1}\Big)^k. \tag{5.8}$$
>
> **Conclusión.** De (5.7) y (5.8), $\max_{\lambda\in\sigma(A)}|p_k(\lambda)|\leq2\big(\frac{\sqrt c-1}{\sqrt c+1}\big)^k$, y aplicando el [[Corolario 5.14]] con este $p_k$ se obtiene el resultado.

> [!demostracion]- Identidad auxiliar: $T_k\big(\frac{z+z^{-1}}2\big)=\frac{z^k+z^{-k}}2$ para todo $z\neq0$
> Por inducción en $k$, usando la misma recurrencia que define a $T_k$. Sea $x:=\frac{z+z^{-1}}2$.
>
> **Casos base.** $k=0$: $T_0(x)=1=\frac{z^0+z^0}2=\frac{1+1}2=1$ ✓. $k=1$: $T_1(x)=x=\frac{z+z^{-1}}2$ ✓ (por definición de $x$).
>
> **Paso inductivo.** Suponiendo la identidad para $k-1$ y $k-2$:
> $$T_k(x) = 2x\,T_{k-1}(x)-T_{k-2}(x) = 2\cdot\frac{z+z^{-1}}2\cdot\frac{z^{k-1}+z^{-(k-1)}}2 - \frac{z^{k-2}+z^{-(k-2)}}2.$$
> Expandiendo el primer producto, $(z+z^{-1})(z^{k-1}+z^{1-k}) = z^k+z^{2-k}+z^{k-2}+z^{-k}$, así que
> $$T_k(x) = \frac{z^k+z^{2-k}+z^{k-2}+z^{-k}}2 - \frac{z^{k-2}+z^{2-k}}2 = \frac{z^k+z^{-k}}2,$$
> donde los términos $z^{2-k}$ y $z^{k-2}$ se cancelan entre el primer y el segundo sumando. Por inducción, la identidad vale para todo $k\geq0$.

Comparar con steepest descent: $\frac{\kappa_2(A)-1}{\kappa_2(A)+1}$ - la diferencia entre $\kappa$ y $\sqrt\kappa$ en el exponente es la mejora clave (ver [[Tablas Comparativas - GC vs SD]]).

## Ver también
- [[Cotas de error del gradiente conjugado]] (nota de tema completa)
- [[Teorema 5.13]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
