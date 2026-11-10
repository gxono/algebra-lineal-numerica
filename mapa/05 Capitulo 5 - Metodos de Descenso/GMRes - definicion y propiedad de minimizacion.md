---
tags: [capitulo5, teorema, definicion]
---

# GMRes: definición y propiedad de minimización

> [!idea]
> La generalización natural del GC a matrices **no simétricas**: en vez de minimizar el error en norma $A$ (que exige $A$ sdp para que $‖·‖_A$ sea norma), GMRes minimiza el **residuo en norma 2**, que tiene sentido para cualquier $A$ invertible. El precio: sin simetría, no hay recurrencia de dos términos - hay que guardar toda la base de Krylov.

## Definición

$$
x_k = \operatorname*{argmin}_{x\in x_0+\mathcal K_k}\|b-Ax\|_2, \qquad \mathcal K_k=\text{span}\{r_0,Ar_0,\dots,A^{k-1}r_0\}.
$$

## Teorema 5.36 (minimización vía polinomios residuales)

$$
\|r_k\|_2 = \min_{p\in\overline{\mathcal P}_k}\|p(A)r_0\|_2 \implies \|r_k\|_2\leq\|p(A)\|_2\|r_0\|_2 \quad\forall p\in\overline{\mathcal P}_k.
$$

Mismo tipo de argumento que el [[Teorema 5.13]] de GC, pero con la norma 2 del residuo en vez de la norma $A$ del error - acá no hace falta $A$ sdp porque no se usa el Lema 5.3. Demostración completa: [[Teorema 5.36]].

## Terminación finita (Teorema 5.38)

GMRes encuentra la solución exacta en a lo sumo $n$ pasos. Demostración completa (Cayley-Hamilton + Teorema 5.36): [[Teorema 5.38]].

## Cotas para matrices diagonalizables (Teorema 5.39)

Si $A=V\Lambda V^{-1}$ (posiblemente compleja):
$$
\frac{\|r_k\|_2}{\|r_0\|_2} \leq \kappa_2(V)\max_{z\in\sigma(A)}|p(z)|, \qquad \forall p\in\overline{\mathcal P}_k.
$$

Demostración completa: [[Teorema 5.39]].

El **precio de la no simetría**: aparece $\kappa_2(V)$ (cuán lejos está $A$ de ser normal). Si $A$ es normal ($AA^T=A^TA$, en particular simétrica u ortogonal), $V$ es ortogonal y $\kappa_2(V)=1$ - se recupera una cota tan limpia como la de GC (de hecho, para $A$ sdp esto reconstruye exactamente el Corolario 5.14, vía $\|p(A)\|_2=\max_\lambda|p(\lambda)|$). Si $A$ está cerca de no ser diagonalizable, $\kappa_2(V)$ puede ser grande y la cota, poco informativa (aunque GMRes puede seguir funcionando bien en la práctica).

## Terminación con pocos autovalores distintos (Corolarios 5.40, 5.41)

Análogos exactos de los [[GC como metodo directo|Teoremas 5.21/5.22]] de GC: $k$ autovalores distintos ⟹ termina en $≤k$ pasos; $b$ combinación de $k$ autovectores, $x_0=0$ ⟹ termina en $≤k$ pasos.

## Criterio de parada

$$
\|r_k\|_2<\text{tol}\|b\|_2 \implies \frac{\|e_k\|_2}{\|x_*\|_2}\leq\kappa_2(A)\,\text{tol}.
$$

## Ver también
- [[Gradiente conjugado - definicion y propiedades]]
- [[GMRes - precondicionamiento]]
- [[Tablas Comparativas - GC vs GMRes]]
- [[5.19 - Demostrar Corolario 5.40|Problema 5.19]], [[5.20 - Demostrar Corolario 5.41|Problema 5.20]]
