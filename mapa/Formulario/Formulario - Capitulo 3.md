---
tags: [capitulo3, moc]
---

# Formulario - Capítulo 3 (Factorización QR)

## Mínimos cuadrados y condicionamiento

$$
\min_x\|Ax-b\|_2, \qquad \kappa_2(A)=\frac{\sigma_1}{\sigma_m}, \qquad \kappa_2(A^TA)=\kappa_2(A)^2.
$$
Error hacia adelante $\cong \kappa_2(A)\times$ error hacia atrás.

## Ecuaciones normales

$$
A^TAx = A^Tb \qquad \text{(resolver con Cholesky, costo} \cong nm^2\text{, error} \cong\kappa_2(A)^2\varepsilon_{maq}).
$$

## [[Ortogonalizacion de Gram-Schmidt|Gram-Schmidt]]

$$
u_1=\frac{v_1}{\|v_1\|}, \qquad w_k = v_k-\sum_{i<k}(v_k\cdot u_i)u_i, \qquad u_k=\frac{w_k}{\|w_k\|}.
$$
Usar la versión **modificada** (restar proyecciones de a una, recalculando).

## [[Matrices de Householder|Householder]]

$$
H = I-uu^T, \quad \|u\|_2=\sqrt2, \qquad u = \frac{1}{\sqrt{1+\sigma x_1}}(x+\sigma e_1), \quad \sigma=\text{sgn}(x_1), \qquad Hx=-\sigma\|x\|_2e_1.
$$
Costo $\cong2nm^2$, error $\cong\kappa_2(A)\varepsilon_{maq}$ (óptimo).

## [[Rotaciones de Givens|Givens]] (para Hessenberg)

$$
G=\begin{bmatrix}c&-s\\s&c\end{bmatrix}, \qquad c=\frac{h_{11}}{\sqrt{h_{11}^2+h_{21}^2}}, \qquad s=\frac{-h_{21}}{\sqrt{h_{11}^2+h_{21}^2}}.
$$
Anula una entrada por rotación, $O(1)$ por columna que toca; QR de Hessenberg en $O(n^2)$ (vs. $O(n^3)$ denso).

## Costo y estabilidad, comparados

| Método | Costo ($n\gg m$) | Error hacia adelante |
|---|---|---|
| Ecuaciones normales | $\cong nm^2$ | $\cong\kappa_2(A)^2\varepsilon_{maq}$ |
| Gram-Schmidt | $\cong2nm^2$ | intermedio |
| Householder | $\cong2nm^2$ | $\cong\kappa_2(A)\varepsilon_{maq}$ |

## Ver también
- [[Formulario - Indice]]
- [[Capitulo 3 - MOC]]
