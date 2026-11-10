---
tags: [capitulo4, teorema]
---

# Teorema 4.36 (cota en norma 2 de Richardson óptimo)

Si $A$ es sdp y se aplica Richardson con $\alpha_{\text{opt}}=\dfrac{2}{\lambda_{\min}+\lambda_{\max}}$ ([[Lema 4.35]]), a partir de cualquier $x_0$:
$$
\|x-x_{k+1}\|_2 \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\,\|x-x_k\|_2, \qquad k=0,1,2,\dots
$$

> [!demostracion]-
> El error satisface $x-x_{k+1}=M_{\alpha_{\text{opt}}}(x-x_k)$ con $M_{\alpha_{\text{opt}}}=I-\alpha_{\text{opt}}A$ (igual que en [[Teorema de convergencia de metodos iterativos estacionarios]]: $x$ es el punto fijo, $x-x_{k+1}=M(x-x_k)$). Así,
> $$\|x-x_{k+1}\|_2 = \|M_{\alpha_{\text{opt}}}(x-x_k)\|_2 \leq \|M_{\alpha_{\text{opt}}}\|_2\,\|x-x_k\|_2.$$
> Como $A$ es simétrica, $M_{\alpha_{\text{opt}}}=I-\alpha_{\text{opt}}A$ también es simétrica, y por [[4.5(d) - Formula norma 2 matriz simetrica|Problema 4.5(d)]], $\|M_{\alpha_{\text{opt}}}\|_2=\rho(M_{\alpha_{\text{opt}}})$. Por el [[Lema 4.35]], $\rho(M_{\alpha_{\text{opt}}})=\dfrac{\kappa_2(A)-1}{\kappa_2(A)+1}$. Sustituyendo,
> $$\|x-x_{k+1}\|_2 \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\,\|x-x_k\|_2.$$

> [!idea]
> Para matrices simétricas, la norma 2 inducida **coincide** con el radio espectral ($\|M\|_2=\rho(M)$) - por eso esta cota en norma 2 sale directo del [[Lema 4.35]] sin trabajo adicional. Eso deja de ser cierto en norma $A$ (donde $M_{\alpha_{\text{opt}}}$ no es simétrica respecto al producto $A$-ortogonal de la misma manera trivial) - de ahí que el [[Teorema 4.37]] necesite un paso extra.

## Ver también
- [[Lema 4.35]]
- [[Teorema 4.37]] (misma cota, en norma $A$)
- [[Metodo de Richardson]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
