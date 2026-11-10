---
tags: [capitulo2, teorema]
---

# Teorema de existencia de Cholesky (factorización $A=R^TR$)

> [!demostracion]- Si $A$ es sdp, existe $R$ triangular superior con $A=RᵀR$
> **Demostración constructiva, por inducción en $n$.** Escribiendo $A=\begin{bmatrix}α&aᵀ\\a&A_*\end{bmatrix}$ y buscando $R=\begin{bmatrix}ρ&rᵀ\\0&R_*\end{bmatrix}$ tal que $RᵀR=A$:
> $$R^TR = \begin{bmatrix}\rho^2&\rho r^T\\\rho r&R_*^TR_*+rr^T\end{bmatrix} \stackrel{!}{=} \begin{bmatrix}\alpha&a^T\\a&A_*\end{bmatrix}.$$
> Las dos primeras condiciones se satisfacen con $ρ=√α>0$ (recordar $α=a₁₁>0$ por ser $A$ sdp) y $r=a/ρ$. Falta hallar $R_*$ triangular superior con $R_*^TR_* = A_*-rr^T$; por inducción, alcanza con ver que $A_*-rr^T$ es sdp.
>
> Dado $y≠0$, se busca $x=(η,y)$ tal que $xᵀAx = yᵀ(A_*-rrᵀ)y$. Expandiendo $xᵀAx = αη²+2ηaᵀy+yᵀA_*y$ e igualando a $yᵀA_*y-(rᵀy)²$, la ecuación para $η$ resulta $(ρη+rᵀy)²=0$, que se resuelve con $η=-rᵀy/ρ$. Con ese $x$, $xᵀAx=yᵀ(A_*-rrᵀ)y$, y como $A$ es sdp y $x≠0$ (porque $y≠0$), $xᵀAx>0$. Luego $A_*-rrᵀ$ es sdp, y por hipótesis inductiva admite factor de Cholesky $R_*$.

**Observación:** la factorización con diagonal positiva es **única**. Además, $A=RᵀR = (RᵀD⁻¹)(DR) = LU$ con $D=\text{diag}(R)$: Cholesky es un caso particular de LU (con $L=Rᵀ$, salvo el reescalado por $D$).

## Ver también
- [[Factorizacion de Cholesky]] (nota de tema completa)
- [[Demostraciones - Capitulo 2|Demostraciones]]
