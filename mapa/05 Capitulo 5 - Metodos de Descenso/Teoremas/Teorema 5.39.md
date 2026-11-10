---
tags: [capitulo5, teorema]
---

# Teorema 5.39 (cotas para matrices diagonalizables)

Si $A=V\Lambda V^{-1}$ (posiblemente compleja):
$$
\frac{\|r_k\|_2}{\|r_0\|_2} \leq \kappa_2(V)\max_{z\in\sigma(A)}|p(z)|, \qquad \forall p\in\overline{\mathcal P}_k.
$$

> [!demostracion]-
> Sea $p\in\overline{\mathcal P}_k$. Como $A=V\Lambda V^{-1}$, $p(A)=Vp(\Lambda)V^{-1}$ (expandiendo la serie de $p$ término a término: $A^j=V\Lambda^jV^{-1}$). Entonces, usando submultiplicatividad de $\|\cdot\|_2$:
> $$\|p(A)\|_2 = \|Vp(\Lambda)V^{-1}\|_2 \leq \|V\|_2\|V^{-1}\|_2\|p(\Lambda)\|_2 = \kappa_2(V)\max_{z\in\sigma(A)}|p(z)|,$$
> donde $\|p(\Lambda)\|_2=\max_z|p(z)|$ porque $p(\Lambda)$ es diagonal con esas entradas. Por el Corolario 5.37 (consecuencia inmediata del Teorema 5.36), $\|r_k\|_2/\|r_0\|_2\leq\|p(A)\|_2$, y combinando con lo anterior se obtiene el resultado.

El **precio de la no simetría**: aparece $\kappa_2(V)$ (cuán lejos está $A$ de ser normal). Si $A$ es normal ($AA^T=A^TA$, en particular simétrica u ortogonal), $V$ es ortogonal y $\kappa_2(V)=1$ - se recupera una cota tan limpia como la de GC (de hecho, para $A$ sdp esto reconstruye exactamente el [[Corolario 5.14]], vía $\|p(A)\|_2=\max_\lambda|p(\lambda)|$). Si $A$ está cerca de no ser diagonalizable, $\kappa_2(V)$ puede ser grande y la cota, poco informativa (aunque GMRes puede seguir funcionando bien en la práctica).

## Ver también
- [[GMRes - definicion y propiedad de minimizacion]] (nota de tema completa)
- [[Teorema 5.36]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
