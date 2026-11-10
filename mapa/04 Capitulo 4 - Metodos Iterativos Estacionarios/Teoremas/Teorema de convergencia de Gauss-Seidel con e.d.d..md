---
tags: [capitulo4, teorema]
---

# Teorema: e.d.d. $\implies$ Gauss-Seidel converge

> [!demostracion]- Si $A$ es e.d.d., entonces $ρ(M_{GS})<1$
> Sea $\lambda$ autovalor de $M_{GS}$ con autovector $v\neq0$, y escalemos $v$ para que $\|v\|_\infty=1$, con $|v_k|=1$ para algún $k$ (la componente donde se alcanza el máximo). Si $\lambda=0$, trivialmente $|\lambda|<1$; supongamos $\lambda\neq0$.
>
> De $-(L+D)^{-1}Uv=\lambda v$ se sigue $\lambda(L+D)v+Uv=0$, y la fila $k$ da
> $$\lambda a_{kk}v_k = -\lambda\sum_{j<k}a_{kj}v_j - \sum_{j>k}a_{kj}v_j.$$
> Tomando módulo, con $|v_k|=1$ y $|v_j|\leq1$ para todo $j$:
> $$|\lambda||a_{kk}| \leq |\lambda|\sum_{j<k}|a_{kj}|+\sum_{j>k}|a_{kj}|.$$
> Supongamos, para llegar a un absurdo, que $|\lambda|\geq1$. Entonces $\sum_{j>k}|a_{kj}|\leq|\lambda|\sum_{j>k}|a_{kj}|$, así que
> $$|\lambda||a_{kk}| \leq |\lambda|\Big(\sum_{j<k}|a_{kj}|+\sum_{j>k}|a_{kj}|\Big) = |\lambda|\sum_{j\neq k}|a_{kj}|.$$
> Dividiendo por $|\lambda|>0$: $|a_{kk}|\leq\sum_{j\neq k}|a_{kj}|$, que **contradice** que $A$ sea e.d.d. en la fila $k$. Luego $|\lambda|<1$. Como $\lambda$ era un autovalor arbitrario, $\rho(M_{GS})<1$.

> [!pendiente]
> Esta demostración es más sutil que la de Jacobi: la prueba de **esta** página no alcanza con acotar una norma de $M_{GS}$ directamente (no es tan inmediato como $‖M_J‖_∞$, que se lee de un tirón de la definición de e.d.d.), hace falta trabajar con el autovector y la fila $k$ donde se alcanza el máximo. Vale la pena repasarla paso a paso antes del parcial.
>
> **Aclaración importante:** esto no quiere decir que sea *imposible* acotar $\|M_{GS}\|_\infty$ directamente - sí se puede, con un argumento distinto (inducción fila por fila), y da **la misma cota** $\alpha$ que Jacobi. Ver [[Cota de norma infinito de Gauss-Seidel para matrices e.d.d.]]. La demostración de esta página prueba *solo* $\rho(M_{GS})<1$ (convergencia cualitativa); si además hace falta una *tasa* explícita (cuántas iteraciones), usar esa otra página.

## Ver también
- [[Metodo de Gauss-Seidel]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
