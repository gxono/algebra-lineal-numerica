---
tags: [capitulo5, teorema]
---

# Teorema 5.36 (minimización vía polinomios residuales)

$$
\|r_k\|_2 = \min_{p\in\overline{\mathcal P}_k}\|p(A)r_0\|_2 \implies \|r_k\|_2\leq\|p(A)\|_2\|r_0\|_2 \quad\forall p\in\overline{\mathcal P}_k.
$$

> [!demostracion]-
> Para $x\in x_0+\mathcal K_k$: $x=x_0+\sum_{j=0}^{k-1}\gamma_jA^jr_0$, así que
> $$r=b-Ax=b-Ax_0-\sum_{j=0}^{k-1}\gamma_jA^{j+1}r_0=r_0-\sum_{j=1}^k\gamma_{j-1}A^jr_0=p(A)r_0,$$
> con $p(z):=1-\sum_{j=1}^k\gamma_{j-1}z^j\in\overline{\mathcal P}_k$. Recíprocamente todo $p\in\overline{\mathcal P}_k$ se obtiene así. Como $x_k$ minimiza $\|b-Ax\|_2=\|r\|_2$ sobre $x\in x_0+\mathcal K_k$, minimizar sobre $x$ es lo mismo que minimizar $\|p(A)r_0\|_2$ sobre $p\in\overline{\mathcal P}_k$, de donde $\|r_k\|_2=\min_p\|p(A)r_0\|_2$. La desigualdad es inmediata: para cualquier $p$ particular, $\|p(A)r_0\|_2\leq\|p(A)\|_2\|r_0\|_2$, y el mínimo es $\leq$ cualquier valor particular.

(Mismo tipo de argumento que el [[Teorema 5.13]] de GC, pero con la norma 2 del residuo en vez de la norma $A$ del error - acá no hace falta $A$ sdp porque no se usa el Lema 5.3.)

## Ver también
- [[GMRes - definicion y propiedad de minimizacion]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
