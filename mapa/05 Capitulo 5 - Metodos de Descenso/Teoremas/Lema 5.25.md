---
tags: [capitulo5, teorema]
---

# Lema 5.25

Sea $k\in\mathbb N$. Si $r_{k-1}\neq0$, entonces $\{d_1,\dots,d_k\}$ es $A$-ortogonal y
$$
\mathcal K_k=\text{span}\{d_1,\dots,d_k\}.
$$

> [!demostracion]- Por inducción en $k$
> **Caso base $k=1$:** $d_1=r_0$, así que $\mathcal K_1=\text{span}\{r_0\}=\text{span}\{d_1\}$, trivial.
>
> **Paso inductivo:** supongamos $\mathcal K_k=\text{span}\{d_1,\dots,d_k\}$ con $\{d_1,\dots,d_k\}$ $A$-ortogonal (H.I.), y veamos que si $r_k\neq0$, $\mathcal K_{k+1}=\text{span}\{d_1,\dots,d_{k+1}\}$ con $\{d_1,\dots,d_{k+1}\}$ $A$-ortogonal.
>
> Por H.I., $\{d_1,\dots,d_k\}\subset\mathcal K_k\subset\mathcal K_{k+1}$. Por definición de $d_{k+1}$ ([[Corolario 5.24]]), $d_{k+1}\in\mathcal K_{k+1}$ y $d_{k+1}$ es $A$-ortogonal a $\mathcal K_k=\text{span}\{d_1,\dots,d_k\}$ (H.I.) - en particular $A$-ortogonal a cada $d_i$, $i\leq k$. Junto con la $A$-ortogonalidad de $\{d_1,\dots,d_k\}$ (H.I.), esto da que $\{d_1,\dots,d_{k+1}\}$ es $A$-ortogonal, y
> $$\text{span}\{d_1,\dots,d_{k+1}\}\subset\mathcal K_{k+1}.$$
> Como $r_k\neq0$, los vectores $d_1,\dots,d_{k+1}$ son todos no nulos (cada $d_{i+1}$ es $A$-ortogonal a $d_i\neq0$ con $A$ sdp, así que $d_{i+1}=0$ forzaría... en rigor: $d_{i+1}\neq0$ porque es la dirección que lleva de $x_i$ a $x_{i+1}\neq x_i$, ya que $r_i\neq0$ implica que $\phi$ todavía no está minimizada en $x_0+\mathcal K_i$). Vectores $A$-ortogonales no nulos son linealmente independientes (si $A$ es sdp, $A$-ortogonalidad es ortogonalidad respecto al producto interno $(\cdot,\cdot)_A$ - ver [[Producto escalar A y norma A]]), así que $\text{span}\{d_1,\dots,d_{k+1}\}$ tiene dimensión $k+1$. Pero $\mathcal K_{k+1}$ tiene dimensión a lo sumo $k+1$ y lo contiene - luego son iguales.

## Ver también
- [[Implementacion del gradiente conjugado]] (nota de tema completa)
- [[Corolario 5.24]]
- [[Teorema 5.26]]
- [[Producto escalar A y norma A]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
