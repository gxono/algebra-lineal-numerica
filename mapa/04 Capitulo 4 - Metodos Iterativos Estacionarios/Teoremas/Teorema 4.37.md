---
tags: [capitulo4, teorema]
---

# Teorema 4.37 (cota en norma $A$ de Richardson óptimo)

Si $A$ es sdp y se aplica Richardson con $\alpha_{\text{opt}}=\dfrac{2}{\lambda_{\min}+\lambda_{\max}}$ ([[Lema 4.35]]), a partir de cualquier $x_0$:
$$
\|x-x_{k+1}\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\,\|x-x_k\|_A, \qquad k=0,1,2,\dots,
$$
donde $\|y\|_A=(y^TAy)^{1/2}$ ([[Producto escalar A y norma A]]).

> [!demostracion]-
> Como $A$ es sdp, es diagonalizable ortogonalmente: existe una base ortonormal de autovectores $\{v_1,\dots,v_n\}$ con $Av_i=\lambda_iv_i$. Para $\alpha$ cualquiera, $M_\alpha:=I-\alpha A$ tiene la **misma** base de autovectores, con $M_\alpha v_i=(1-\alpha\lambda_i)v_i$ - es decir, $A$ y $M_\alpha$ se diagonalizan simultáneamente (conmutan, porque $M_\alpha$ es un polinomio en $A$).
>
> Sea $e:=x-x_k$ un vector cualquiera, escrito en esa base como $e=\sum_ic_iv_i$. Entonces:
> $$\|e\|_A^2 = e^TAe = \Big(\sum_ic_iv_i\Big)^TA\Big(\sum_jc_jv_j\Big) = \sum_i\lambda_ic_i^2$$
> (usando $v_i^TAv_j=\lambda_jv_i^Tv_j=\lambda_j\delta_{ij}$ por ortonormalidad). Y $M_\alpha e=\sum_i(1-\alpha\lambda_i)c_iv_i$, así que
> $$\|M_\alpha e\|_A^2 = \sum_i\lambda_i(1-\alpha\lambda_i)^2c_i^2 \leq \Big(\max_i(1-\alpha\lambda_i)^2\Big)\sum_i\lambda_ic_i^2 = \rho(M_\alpha)^2\,\|e\|_A^2.$$
> Es decir, **para cualquier** $\alpha$ y cualquier vector $e$: $\|M_\alpha e\|_A\leq\rho(M_\alpha)\|e\|_A$.
>
> Aplicando esto con $\alpha=\alpha_{\text{opt}}$ y $e=x-x_k$ (recordando que $x-x_{k+1}=M_{\alpha_{\text{opt}}}(x-x_k)$, igual que en el [[Teorema 4.36]]):
> $$\|x-x_{k+1}\|_A = \|M_{\alpha_{\text{opt}}}(x-x_k)\|_A \leq \rho(M_{\alpha_{\text{opt}}})\,\|x-x_k\|_A = \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\,\|x-x_k\|_A,$$
> usando el valor de $\rho(M_{\alpha_{\text{opt}}})$ del [[Lema 4.35]].

> [!idea]
> El punto clave es que $A$ y $M_\alpha$ **comparten autovectores** (porque $M_\alpha$ es un polinomio en $A$) - eso permite trabajar las dos normas ($\|\cdot\|_2$ para el Teorema 4.36 y $\|\cdot\|_A$ para este teorema) con el mismo argumento espectral, sin necesitar construir $A^{1/2}$ explícitamente (que es el camino que toma el apunte). Es exactamente el mismo fenómeno que permite, en el Capítulo 5, acotar el error de GC vía polinomios evaluados en los autovalores de $A$ ([[Cotas de error del gradiente conjugado]]): cualquier polinomio en $A$ actúa diagonalmente en la base de autovectores de $A$.

## Ver también
- [[Lema 4.35]]
- [[Teorema 4.36]] (misma cota, en norma 2)
- [[Metodo de Richardson]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
