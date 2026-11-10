---
tags: [capitulo5, teorema]
---

# Teorema 5.9 (cota de convergencia de descenso más pronunciado)

$$
\|x_{k+1}-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A, \qquad k\geq0.
$$

> [!demostracion]-
> Sea $k\geq0$ fijo. La iteración de Richardson con $\alpha_{\text{opt}}$ (Capítulo 4) da $x_{k+1}^R=x_k-\alpha_{\text{opt}}(Ax_k-b)=x_k+\alpha_{\text{opt}}r_k$. Para SD, $x_{k+1}=x_k+\alpha r_k$ con $\alpha$ elegido para minimizar $h(\alpha):=\phi(x_k+\alpha r_k)$ - es decir, el $\alpha$ de SD es el **mejor posible** sobre esa recta, mientras que $\alpha_{\text{opt}}$ es solo una elección particular. Por el [[Funcional de energia y caracterizacion|Lema 5.3]], minimizar $\phi$ sobre la recta equivale a minimizar $\|\cdot-x_*\|_A$ sobre ella, así que
> $$\|\underbrace{x_k+\alpha r_k}_{x_{k+1}}-x_*\|_A \leq \|\underbrace{x_k+\alpha_{\text{opt}}r_k}_{x_{k+1}^R}-x_*\|_A.$$
> Por el teorema de convergencia de Richardson (Capítulo 4), $\|x_{k+1}^R-x_*\|_A\leq\frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A$, y combinando ambas desigualdades se obtiene el resultado.

Misma cota que Richardson con $α$ óptimo global - la "ventaja" de recalcular $α$ en cada paso no mejora la cota en el peor caso, pero sí el comportamiento típico. (La demostración en sí lo deja claro: SD es, en cada paso, **al menos tan bueno** como Richardson con el mejor $\alpha$ fijo posible, porque SD elige el óptimo para la dirección actual.)

## Ver también
- [[Descenso mas pronunciado]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
