---
tags: [capitulo4, definicion, teorema]
---

# Normas de vectores y matrices

> [!idea]
> Una norma mide "tamaño"; axiomas: $‖x‖≥0$ (y $=0$ solo si $x=0$), $‖αx‖=|α|‖x‖$, $‖x+y‖≤‖x‖+‖y‖$. La norma **inducida** en matrices ($‖A‖=\max_{x\neq0}‖Ax‖/‖x‖$) hereda submultiplicatividad: $‖AB‖≤‖A‖‖B‖$, $‖Ax‖≤‖A‖‖x‖$ - esto es lo que permite acotar cómo un método iterativo amplifica (o no) el error en cada paso.

## Normas de vectores en $ℝⁿ$

$$
\|x\|_1=\sum_i|x_i|, \qquad \|x\|_2=\Big(\sum_ix_i^2\Big)^{1/2}, \qquad \|x\|_\infty=\max_i|x_i|.
$$

**Relaciones (sharp, ver [[4.2 - Normas en Rn - equivalencia|Problema 4.2]]):**
$$
\|x\|_\infty\leq\|x\|_2\leq\|x\|_1, \qquad \|x\|_1\leq\sqrt n\|x\|_2\leq n\|x\|_\infty.
$$

## Norma inducida de una matriz

$$
\|A\| := \max_{x\neq0}\frac{\|Ax\|}{\|x\|} = \max_{\|x\|=1}\|Ax\|.
$$

Es norma en $ℝⁿˣⁿ$ (ver [[4.3 - Norma matricial inducida es norma|Problema 4.3]]), y además:
$$
\|Ax\|\leq\|A\|\|x\|, \qquad \|AB\|\leq\|A\|\|B\| \quad(\text{submultiplicativa, ver [[4.4 - Submultiplicatividad norma matricial|Problema 4.4]]}).
$$

## Fórmulas explícitas (las tres normas inducidas más usadas)

> [!comparacion]
> | Norma | Fórmula | Nombre |
> |---|---|---|
> | $‖A‖₁$ | $\max_j\sum_i\|a_{ij}\|$ | máxima suma por columnas |
> | $‖A‖_∞$ | $\max_i\sum_j\|a_{ij}\|$ | máxima suma por filas |
> | $‖A‖₂$ | $\max\{\sqrt\lambda:\lambda\in\sigma(A^TA)\}$ | mayor valor singular |
> | $‖A‖₂$ (si $A$ simétrica) | $\max\{|\lambda|:\lambda\in\sigma(A)\}$ | $=\rho(A)$ |

Demostraciones completas en [[4.5 - Formulas normas matriciales|Problema 4.5]].

## Ver también
- [[Espectro y radio espectral]]
- [[Numero de condicion y error relativo]]
