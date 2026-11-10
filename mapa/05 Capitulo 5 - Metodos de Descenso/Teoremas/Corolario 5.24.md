---
tags: [capitulo5, teorema]
---

# Corolario 5.24

Sea $d_1=r_0$ y, para $k\geq1$, sea $d_{k+1}$ la única dirección en $\mathcal K_{k+1}$ que es $A$-ortogonal a $\mathcal K_k$ (existe y es única por [[Lema 5.23]], ya que $\dim\mathcal K_{k+1}-\dim\mathcal K_k=1$). Entonces, si $k\in\mathbb N_0$ y $r_k\neq0$,
$$
x_{k+1}=x_k+\alpha_{k+1}d_{k+1}, \qquad \alpha_{k+1}=\frac{d_{k+1}^Tr_k}{d_{k+1}^TAd_{k+1}}.
$$

> [!demostracion]-
> Sea $k\in\mathbb N$. Por el [[Lema 5.23]], $x_{k+1}-x_k\in\mathcal K_{k+1}$ y $x_{k+1}-x_k\perp_A\mathcal K_k$ - es decir, $x_{k+1}-x_k$ es **la** dirección en $\mathcal K_{k+1}$ que es $A$-ortogonal a $\mathcal K_k$ (salvo escala). Por unicidad de esa dirección (a escala), existe $\alpha_{k+1}\in\mathbb R$ tal que
> $$x_{k+1}-x_k=\alpha_{k+1}d_{k+1}, \qquad\text{es decir,}\qquad x_{k+1}=x_k+\alpha_{k+1}d_{k+1}.$$
> Como $x_{k+1}=\operatorname*{argmin}_{x\in x_0+\mathcal K_{k+1}}\phi(x)$ y $x_k+\alpha d_{k+1}$ recorre $x_0+\mathcal K_{k+1}$ al variar $\alpha$ (porque $x_k-x_0\in\mathcal K_k\subset\mathcal K_{k+1}$ y $d_{k+1}\in\mathcal K_{k+1}$), $\alpha_{k+1}$ debe minimizar la función de una variable $\alpha\mapsto\phi(x_k+\alpha d_{k+1})$.
>
> Esa minimización es puro cálculo de una variable: con $\phi(x)=\tfrac12x^TAx-b^Tx$,
> $$\frac{d}{d\alpha}\phi(x_k+\alpha d_{k+1}) = d_{k+1}^T\big(A(x_k+\alpha d_{k+1})-b\big) = \alpha\, d_{k+1}^TAd_{k+1} - d_{k+1}^Tr_k,$$
> que es una recta en $\alpha$ con pendiente $d_{k+1}^TAd_{k+1}>0$ (porque $A$ es sdp y $d_{k+1}\neq0$) - así que se anula en un único punto, que es el mínimo:
> $$\alpha_{k+1}=\frac{d_{k+1}^Tr_k}{d_{k+1}^TAd_{k+1}}.$$

## Ver también
- [[Implementacion del gradiente conjugado]] (nota de tema completa)
- [[Lema 5.23]]
- [[Lema 5.27]] (fórmula simplificada, solo en términos de residuos)
- [[Demostraciones - Capitulo 5|Demostraciones]]
