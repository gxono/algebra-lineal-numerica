---
tags: [capitulo4, teorema]
---

# Lema 4.10 (número de condición y error relativo)

Si $\tilde x$ resuelve $A\tilde x=\tilde b$ (aproximadamente, $\tilde b\approx b$) y $x$ resuelve $Ax=b$:
$$
\frac{1}{\kappa(A)}\frac{\|\tilde b-b\|}{\|b\|} \leq \frac{\|\tilde x-x\|}{\|x\|} \leq \kappa(A)\frac{\|\tilde b-b\|}{\|b\|}.
$$

> [!demostracion]-
> De $Ax=b$ y $A\tilde x=\tilde b$: $A(\tilde x-x)=\tilde b-b$, así que $\tilde x-x=A^{-1}(\tilde b-b)$, de donde $\|\tilde x-x\|\leq\|A^{-1}\|\|\tilde b-b\|$. Por otro lado $\|b\|=\|Ax\|\leq\|A\|\|x\|$. Combinando:
> $$\frac{\|\tilde x-x\|}{\|x\|} \leq \frac{\|A^{-1}\|\|\tilde b-b\|}{\|x\|} \leq \|A^{-1}\|\|A\|\frac{\|\tilde b-b\|}{\|b\|} = \kappa(A)\frac{\|\tilde b-b\|}{\|b\|}.$$
> La cota inferior es análoga, intercambiando los roles de $x$ y $\tilde x$.

**Interpretación:** si el residuo relativo $\|r_k\|/\|b\|$ es chico, el error relativo $\|x_k-x_*\|/\|x_*\|$ puede seguir siendo grande si $κ(A)$ es grande - *residuo chico no implica error chico* si la matriz está mal condicionada. Este es el motivo por el que todas las cotas de convergencia del Capítulo 5 (GC, [[GMRes - definicion y propiedad de minimizacion|GMRes]]) terminan con un factor $κ(A)$ o $κ_2(A)$ multiplicando.

## Ver también
- [[Numero de condicion y error relativo]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
