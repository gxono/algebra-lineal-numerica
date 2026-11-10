---
tags: [capitulo5, teorema]
---

# Teorema 5.38 (terminación finita)

GMRes encuentra la solución exacta en a lo sumo $n$ pasos.

> [!demostracion]-
> Sea $p(z):=\det(A-zI)/\det(A)$ (el polinomio característico de $A$, normalizado: $p(0)=\det(A)/\det(A)=1$, usando que $A$ no singular $\Rightarrow\det(A)\neq0$). Por Cayley-Hamilton, $\det(A-zI)$ evaluado en $A$ da la matriz nula, así que $p(A)=0$. Aplicando el [[Teorema 5.36]] con este $p\in\overline{\mathcal P}_n$: $\|r_n\|_2\leq\|p(A)\|_2\|r_0\|_2=0$, es decir $r_n=0$ y $x_n=x_*$.

## Ver también
- [[GMRes - definicion y propiedad de minimizacion]] (nota de tema completa)
- [[Teorema 5.36]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
