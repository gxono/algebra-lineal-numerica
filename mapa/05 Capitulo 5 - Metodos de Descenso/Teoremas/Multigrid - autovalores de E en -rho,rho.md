---
tags: [capitulo5, teorema]
---

# Multigrid: autovalores de $E$ en $[-\rho,\rho]$ (Paso 1)

**Paso 1 (multigrid es un método estacionario).** Un ciclo en V con datos $(u,f)$ se puede escribir $u\mapsto u+M(f-A_hu)$ para cierta matriz $M$ (tomando $u=0$, $f=\xi$ identifica $M$ como "un ciclo en V con dato inicial nulo"). Si $u_*$ es la solución exacta, restando $u_*=u_*+M(f-A_hu_*)$ de la iteración se obtiene que el error nuevo es $E=I-MA_h$ por el error viejo - exactamente un método estacionario del Capítulo 4, con matriz de iteración $E$.

> [!demostracion]- $E$ tiene autovalores reales en $[-\rho,\rho]$, y por lo tanto $\sigma(MA_h)\subset[1-\rho,1+\rho]$
> $MA_h$ es autoadjunta respecto del producto interno $\langle u,v\rangle_{A_h}:=u^TA_hv$: si $M$ y $A_h$ son simétricas,
> $$\langle MA_hu,v\rangle_{A_h} = (MA_hu)^TA_hv = u^TA_hMA_hv = \langle u,MA_hv\rangle_{A_h}$$
> (usando $A_h^T=A_h$ y $M^T=M$). Lo mismo vale para $E=I-MA_h$. Una matriz autoadjunta respecto de un producto interno es diagonalizable con autovalores reales, y su radio espectral coincide con su norma inducida por ese producto interno. Por (5.25), $\|E\|_{A_h}\leq\rho$, así que $|\lambda|\leq\rho$ para todo $\lambda\in\sigma(E)$. Como $\sigma(E)=\sigma(I-MA_h)=\{1-\nu:\nu\in\sigma(MA_h)\}$, esto equivale a $\sigma(MA_h)\subset[1-\rho,1+\rho]$.

## Ver también
- [[Multigrid en 1D - transferencia e implementacion]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
