---
tags: [capitulo5, teorema]
---

# Lema 5.23

Para todo $k\in\mathbb N$,
$$
x_{k+1}-x_k\in\mathcal K_{k+1}, \qquad (x_{k+1}-x_k)^TA\xi=0 \quad \forall\,\xi\in\mathcal K_k.
$$

> [!demostracion]-
> Sea $k\in\mathbb N$. Sabemos que $x_k=x_0+w_k$ con $w_k\in\mathcal K_k$, y $x_{k+1}=x_0+w_{k+1}$ con $w_{k+1}\in\mathcal K_{k+1}$ (definición del gradiente conjugado: $x_j\in x_0+\mathcal K_j$). Como $\mathcal K_k\subset\mathcal K_{k+1}$,
> $$x_{k+1}-x_k = \underbrace{w_{k+1}}_{\in\mathcal K_{k+1}} - \underbrace{w_k}_{\in\mathcal K_k\subset\mathcal K_{k+1}} \in \mathcal K_{k+1}.$$
> Además, si $\xi\in\mathcal K_k$,
> $$(x_{k+1}-x_k)^TA\xi = \xi^TA(x_{k+1}-x_k) = \xi^T\big((Ax_{k+1}-b)+(b-Ax_k)\big) = \xi^T(-r_{k+1}+r_k) = \xi^Tr_k-\xi^Tr_{k+1}=0,$$
> ya que tanto $r_k$ como $r_{k+1}$ son ortogonales a $\mathcal K_k=\text{span}\{r_0,\dots,r_{k-1}\}$ (por [[Lema 5.10]]: $r_j\perp\{r_0,\dots,r_{j-1}\}$ para todo $j$, y $k-1<k$).

> [!idea]
> Dicho de otro modo: $x_{k+1}-x_k$ es la **única** dirección en $\mathcal K_{k+1}$ que es $A$-ortogonal a todo $\mathcal K_k$ (existe y es única porque $\dim\mathcal K_{k+1}-\dim\mathcal K_k=1$: hay un único grado de libertad "nuevo" al pasar de $\mathcal K_k$ a $\mathcal K_{k+1}$, módulo $A$-ortogonalidad al resto). Esa dirección es, por definición, $d_{k+1}$ - de ahí sale directamente el [[Corolario 5.24]].

## Ver también
- [[Implementacion del gradiente conjugado]] (nota de tema completa)
- [[Corolario 5.24]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
