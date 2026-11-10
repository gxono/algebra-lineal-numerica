---
tags: [capitulo5, teorema]
---

# Lema 5.10 (residuos ortogonales)

$r_k\perp\{r_0,\dots,r_{k-1}\}$ para todo $k\geq1$.

> [!demostracion]-
> $x_k$ minimiza $\phi$ sobre $x_0+\mathcal{K}_k$, así que para cualquier $\xi\in\mathcal{K}_k$, $t\mapsto\phi(x_k+t\xi)$ tiene mínimo en $t=0$, de donde $\xi^T\nabla\phi(x_k)=-\xi^Tr_k=0$. Como $\xi\in\mathcal{K}_k=\text{span}\{r_0,\dots,r_{k-1}\}$ es arbitrario, $r_k\perp\{r_0,\dots,r_{k-1}\}$.

**Observación 5.11:** si $r_k\neq0$, $\{r_0,\dots,r_k\}$ es linealmente independiente (vectores no nulos, ortogonales dos a dos).

## Ver también
- [[Gradiente conjugado - definicion y propiedades]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
