---
tags: [capitulo5, teorema]
---

# Lema 5.29 (equivalencia de normas)

$$
\sqrt{c_1}\|v\|_{M^{-1}} \leq \|v\|_A \leq \sqrt{c_2}\|v\|_{M^{-1}}.
$$

> [!demostracion]-
> Para $B$ sdp, $v^TB^{-1}v=\max_{x\neq0}\dfrac{(x^Tv)^2}{x^TBx}$ (Cauchy-Schwarz con $x^TB^{-1/2}B^{1/2}$). Aplicando esto con $B=M$ y con $B=A^{-1}$, y usando la hipótesis $c_1x^TA^{-1}x\leq x^TMx\leq c_2x^TA^{-1}x$, se obtiene $\frac1{c_2}\|v\|_A^2\leq\|v\|_{M^{-1}}^2\leq\frac1{c_1}\|v\|_A^2$, que reordenado da la desigualdad de arriba.

Esto permite trasladar la cota de convergencia (naturalmente en norma $M^{-1}$, porque el criterio de parada usa $\|r_k\|_M$) de vuelta a la norma $A$ (en la que está escrita toda la teoría): $\|x_k-x_*\|_A\leq2\big(\frac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\big)^k\|x_0-x_*\|_A$.

## Ver también
- [[Precondicionamiento y GCP]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
