---
tags: [capitulo2, moc]
---

# Formulario - Capítulo 2 (Métodos directos)

## Sustitución

$$
x_i = \frac{1}{\ell_{ii}}\Big(b_i-\sum_{j<i}\ell_{ij}x_j\Big) \quad\text{(adelante)}, \qquad x_i = \frac{1}{u_{ii}}\Big(b_i-\sum_{j>i}u_{ij}x_j\Big) \quad\text{(atrás)}.
$$

## Gauss / [[Factorizacion LU|LU]]

$$
m_{kj} = \frac{a_{kj}^{(j-1)}}{a_{jj}^{(j-1)}}, \qquad A = LU \ \ (\ell_{kj}=m_{kj}), \qquad PA = LU \ \text{(con pivoteo)}.
$$
Existe $LU$ sii todos los menores principales $A_{1:k,1:k}$ son no singulares.

## [[Factorizacion de Cholesky|Cholesky]] ($A$ sdp)

$$
A = R^TR, \qquad \rho = \sqrt{a_{11}}, \qquad r = a/\rho, \qquad A_* \leftarrow A_* - rr^T.
$$
$A$ sdp $\iff A^T=A$ y $x^TAx>0 \ \forall x\neq0$.

## Costo (orden dominante)

| | Costo |
|---|---|
| Sustitución | $n^2$ |
| Gauss / LU | $2n^3/3$ |
| Cholesky | $n^3/3$ |
| Tridiagonal (`tridilu`) | $O(n)$ |
| $k$ sistemas con LU ya factorizado | $2n^3/3+k\cdot n^2$ |

## Ver también
- [[Formulario - Indice]]
- [[Capitulo 2 - MOC]]
