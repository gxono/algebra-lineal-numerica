---
tags: [capitulo3, teorema]
---

# Teorema 3.13: existencia de QR vía Householder

Si $A\in\mathbb{R}^{n\times m}$ tiene columnas li, existen $Q$ ortogonal y $R=\begin{bmatrix}R_1\\0\end{bmatrix}$ (con $R_1$ triangular superior no singular) tales que $A=QR$.

> [!demostracion]- Idea de la inducción en $m$
> Caso $m=1$: $A=x$, se aplica el Teorema 3.9 directamente, $H$ hace de $Q$.
>
> Paso inductivo: con $x=\text{col}_1(A)$, el $H$ del Teorema 3.9 da $HA=\begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat A\end{bmatrix}$. Las columnas de $\hat A$ son li (si no lo fueran, $A$ tampoco tendría $m$ columnas li, ver detalle abajo). Por hipótesis inductiva $\hat A=\hat Q\hat R$, y entonces
> $$A = H\begin{bmatrix}1&0^T\\0&\hat Q\end{bmatrix}\begin{bmatrix}-\sigma\|x\|_2&r^T\\0&\hat R\end{bmatrix}.$$
>
> *(Detalle: si $\hat A\alpha=0$ para algún $\alpha\neq0$, tomando $\beta_1=r^T\alpha/(\sigma\|x\|_2)$ y $\beta=(\beta_1,\alpha)\neq0$, se verifica $HA\beta=0$, y como $H$ es invertible, $A\beta=0$ con $\beta\neq0$: contradice que las columnas de $A$ sean li.)*

## Ver también
- [[Matrices de Householder]] (nota de tema completa)
- [[Teorema 3.9]]
- [[Demostraciones - Capitulo 3|Demostraciones]]
