---
tags: [capitulo3, definicion, teorema]
---

# Matrices ortogonales

> [!idea]
> $Q$ ortogonal ($QᵀQ=I$) preserva ángulos y longitudes: multiplicar por $Q$ no amplifica nada. Por eso los algoritmos basados en matrices ortogonales ([[Matrices de Householder|Householder]], [[Rotaciones de Givens|Givens]]) son estables, y los basados en $AᵀA$ (ecuaciones normales) no lo son tanto.

## Definición 3.1

$Q\in\mathbb{R}^{n\times m}$ es ortogonal si $QᵀQ=I$.

## Propiedades (Observación 3.2)

- Las columnas de $Q$ son ortonormales (de $(QᵀQ)_{ij}=\text{col}_i(Q)\cdot\text{col}_j(Q)=\delta_{ij}$); en particular son linealmente independientes, así que $n\geq m$.
- Producto de ortogonales es ortogonal: $(PQ)ᵀ(PQ)=QᵀPᵀPQ=QᵀQ=I$.
- **$Q$ preserva producto escalar y norma**: $(Qx)\cdot(Qy)=x^TQ^TQy=x^Ty$; en particular $\|Qx\|_2=\|x\|_2$.

## Por qué esto resuelve [[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]] sin perder estabilidad

Si $Q\in\mathbb{R}^{n\times n}$ es ortogonal, minimizar $\|Ax-b\|_2$ equivale a minimizar $\|Q^T(Ax-b)\|_2$ (misma norma, por la propiedad anterior). Si se logra $Q^TA=R=\begin{bmatrix}R_1\\0\end{bmatrix}$ con $R_1$ triangular superior $m\times m$, y se escribe $c=Q^Tb=\binom{c_1}{c_2}$:
$$
\|Q^T(Ax-b)\|_2^2 = \|Rx-c\|_2^2 = \|R_1x-c_1\|_2^2+\|c_2\|_2^2.
$$
El mínimo se logra resolviendo el sistema triangular $R_1x=c_1$, y vale exactamente $\|c_2\|_2^2$. **Alcanza con las primeras $m$ columnas de $Q$ y las primeras $m$ filas de $R$**: el problema se reduce a hallar $Q\in\mathbb{R}^{n\times m}$ ortogonal y $R\in\mathbb{R}^{m\times m}$ triangular superior con $A=QR$.

## Ver también
- [[Ortogonalizacion de Gram-Schmidt]]
- [[Matrices de Householder]]
- [[Minimos cuadrados y ecuaciones normales]]
