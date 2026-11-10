---
tags: [capitulo2, teorema]
---

# Teorema de existencia y unicidad de la factorización LU

> [!demostracion]- Teorema: si $A_{1:k,1:k}$ es no singular para todo $k=1,…,n-1$, existe factorización $A=LU$ (única si $A$ es invertible)
> **Enunciado.** Sea $A$ de $n×n$. Si $\det(A_{1:k,1:k})≠0$ para todo $k=1,…,n-1$, entonces existe $L$ triangular inferior unitaria y $U$ triangular superior con $A=LU$. Si además $A$ es invertible, la factorización es única.
>
> **Demostración (existencia).** Por inducción se ve que $a_{jj}^{(j-1)}≠0$ para todo $j$: para $j=1$, $a_{11}=\det(A_{1:1,1:1})≠0$. Si vale hasta $j=k$, como las $Gʲ$ tienen determinante 1,
> $$\det(A^{(k)}_{1:k+1,1:k+1}) = \det(A_{1:k+1,1:k+1}) = a_{11}^{(0)}a_{22}^{(1)}\cdots a_{(k+1)(k+1)}^{(k)},$$
> y como $A^{(k)}_{1:k+1,1:k+1}$ es triangular superior salvo la última fila, esa igualdad fuerza $a_{(k+1)(k+1)}^{(k)}≠0$.
>
> **Demostración (unicidad).** Si $L₁U₁=L₂U₂=A$ con $Lᵢ$ unitarias y $Uᵢ$ triangulares superiores, entonces $L₂⁻¹L₁=U₂U₁⁻¹$. El lado izquierdo es triangular inferior unitaria y el derecho es triangular superior, así que ambos deben ser $I$, de donde $L₁=L₂$, $U₁=U₂$.

> [!pendiente]
> **No siempre existe factorización LU.** Contraejemplo clásico: $A=\begin{bmatrix}0&1\\1&0\end{bmatrix}$ no tiene factorización LU (ver [[2.7 - Matriz sin factorizacion LU|Problema 2.7]]). La hipótesis de los menores principales no nulos es necesaria, y por eso hace falta el [[Pivoteo y matrices de permutacion|pivoteo]].

## Ver también
- [[Factorizacion LU]] (nota de tema completa)
- [[Demostraciones - Capitulo 2|Demostraciones]]
