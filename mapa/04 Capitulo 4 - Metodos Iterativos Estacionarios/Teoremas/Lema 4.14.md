---
tags: [capitulo4, teorema]
---

# Lema 4.14

Si $A\in\mathbb R^{n\times n}$ es simétrica, entonces
$$
\kappa_2(A) = \frac{\lambda_{\max}}{\lambda_{\min}},
$$
donde $\lambda_{\max}$ y $\lambda_{\min}$ denotan el autovalor de módulo máximo y mínimo de $A$, respectivamente (en valor absoluto).

> [!demostracion]-
> Como $A$ es simétrica, $\|A\|_2=\rho(A)=\max\{|\lambda|:\lambda\in\sigma(A)\}=:\lambda_{\max}$ ([[4.5(d) - Formula norma 2 matriz simetrica|Problema 4.5(d)]]).
>
> $A^{-1}$ también es simétrica (la inversa de una matriz simétrica invertible es simétrica: $(A^{-1})^T=(A^T)^{-1}=A^{-1}$), así que, por la misma fórmula aplicada a $A^{-1}$:
> $$\|A^{-1}\|_2 = \rho(A^{-1}) = \max\{|\mu| : \mu\in\sigma(A^{-1})\}.$$
> Los autovalores de $A^{-1}$ son los recíprocos de los de $A$: si $Av=\lambda v$ con $v\neq0$, $\lambda\neq0$ (porque $A$ es invertible), entonces $A^{-1}v=\lambda^{-1}v$, así que $\sigma(A^{-1})=\{\lambda^{-1}:\lambda\in\sigma(A)\}$. Por lo tanto,
> $$\|A^{-1}\|_2 = \max\{|\lambda|^{-1}:\lambda\in\sigma(A)\} = \frac{1}{\min\{|\lambda|:\lambda\in\sigma(A)\}} = \frac{1}{\lambda_{\min}}.$$
> Juntando ambas partes,
> $$\kappa_2(A) = \|A\|_2\,\|A^{-1}\|_2 = \lambda_{\max}\cdot\frac{1}{\lambda_{\min}} = \frac{\lambda_{\max}}{\lambda_{\min}}.$$

> [!idea]
> Para $A$ sdp (el caso que importa en GC, Capítulo 5) los autovalores ya son positivos, así que $\lambda_{\max}$ y $\lambda_{\min}$ son simplemente el mayor y el menor autovalor - de ahí la notación $\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$ que aparece, sin más comentario, en [[Corolario 5.14]] y en todas las cotas de convergencia del Capítulo 5.

## Ver también
- [[Numero de condicion y error relativo]] (nota de tema completa)
- [[Lema 4.17]] ($\kappa_2(A^TA)=\kappa_2(A)^2$ para $A$ rectangular)
- [[Demostraciones - Capitulo 4|Demostraciones]]
