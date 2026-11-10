---
tags: [capitulo5, teorema]
---

# Corolario 5.14 (cota espectral)

$$
\|x_k-x_*\|_A \leq \Big(\max_{\lambda\in\sigma(A)}|p(\lambda)|\Big)\|x_0-x_*\|_A, \qquad \forall p\in\overline{\mathcal P}_k.
$$

> [!demostracion]-
> Del [[Teorema 5.13]], $\|x_k-x_*\|_A\leq\|p(A)(x_0-x_*)\|_A\leq\|p(A)\|_A\|x_0-x_*\|_A$ para todo $p\in\overline{\mathcal P}_k$ (norma matricial inducida por $\|\cdot\|_A$). Basta entonces probar $\|p(A)\|_A\leq\max_{\lambda\in\sigma(A)}|p(\lambda)|$.
>
> **Paso 1: $\|p(A)\|_A\leq\|p(A)\|_2$.** Para $x\in\mathbb{R}^n$, usando que $p(A)$ (polinomio en $A$) conmuta con $A^{1/2}$:
> $$\|p(A)x\|_A^2 = (p(A)x)^TA(p(A)x) = \|A^{1/2}p(A)x\|_2^2 = \|p(A)A^{1/2}x\|_2^2 \leq \|p(A)\|_2^2\|A^{1/2}x\|_2^2 = \|p(A)\|_2^2\|x\|_A^2,$$
> así que $\|p(A)\|_A\leq\|p(A)\|_2$.
>
> **Paso 2: $\|p(A)\|_2\leq\max_{\lambda\in\sigma(A)}|p(\lambda)|$.** Como $A$ es simétrica, $A=V\Lambda V^T$ con $V$ ortogonal y $\Lambda=\text{diag}(\lambda_1,\dots,\lambda_n)$. Como $p$ es polinomio, $p(A)=Vp(\Lambda)V^T$ con $p(\Lambda)=\text{diag}(p(\lambda_1),\dots,p(\lambda_n))$ (verificar expandiendo $A^{i}=V\Lambda^iV^T$ en la serie de $p$). Entonces $p(A)$ es simétrica con $\sigma(p(A))=\{p(\lambda):\lambda\in\sigma(A)\}$, y para una matriz simétrica $\|p(A)\|_2=\rho(p(A))=\max_{\lambda\in\sigma(A)}|p(\lambda)|$.
>
> Combinando ambos pasos se obtiene el resultado.

> [!idea]
> **Estrategia general** (usada en todos los ejemplos siguientes): elegir **cualquier** $p\in\overline{\mathcal P}_k$ que sea chico sobre el espectro de $A$, sin preocuparse de qué hace fuera de él. El GC, sin "saberlo", siempre hace al menos tan bien como la mejor elección de $p$ - por eso todas estas cotas son cotas *superiores*, nunca el comportamiento exacto.

## Ver también
- [[Cotas de error del gradiente conjugado]] (nota de tema completa)
- [[Teorema 5.13]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
