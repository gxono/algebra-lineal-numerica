---
tags: [capitulo2, teorema]
---

# Teorema: con pivoteo, $PA=LU$

> [!demostracion]- Con pivoteo, $G^{n-1}P^{n-1}\cdots G^1P^1A=U$ triangular superior, y $PA=LU$
> Escribiendo cada paso de Gauss con pivoteo como $A^{(j)}=G^jP^jA^{(j-1)}$, y usando que toda permutación cumple $PP=I$, se puede reagrupar
> $$U = G^{n-1}P^{n-1}\cdots G^2P^2G^1P^1A = \tilde G^{n-1}\cdots\tilde G^1\,\underbrace{P^{n-1}\cdots P^1}_{P}\,A,$$
> con $\tilde Gʲ$ matrices de Gauss "permutadas" (conjugadas por las $P$ posteriores). De ahí $PA=LU$, con $L$ triangular inferior unitaria y $|ℓᵢⱼ|≤1$ (porque el pivote elegido es siempre el de mayor módulo). La columna $k$ de $L$ debajo de la diagonal es una permutación de los multiplicadores del paso $k$.

## Ver también
- [[Pivoteo y matrices de permutacion]] (nota de tema completa)
- [[Demostraciones - Capitulo 2|Demostraciones]]
