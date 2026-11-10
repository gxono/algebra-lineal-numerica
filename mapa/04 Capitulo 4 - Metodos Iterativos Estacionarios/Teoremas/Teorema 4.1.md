---
tags: [capitulo4, teorema]
---

# Teorema 4.1 (desigualdad de Cauchy-Schwarz)

Sea $(\cdot,\cdot)$ un producto escalar en un espacio vectorial $X$ y sea $\|\cdot\|$ la norma que genera ($\|x\|:=(x,x)^{1/2}$). Entonces
$$
|(x,y)|\leq\|x\|\,\|y\|, \qquad \forall\,x,y\in X.
$$

> [!demostracion]-
> Si $y=0$, ambos lados son $0$ y la desigualdad es trivial.
>
> Si $y\neq0$, para todo $t\in\mathbb R$, por no-negatividad del producto escalar aplicado a $x-ty$:
> $$0\leq\|x-ty\|^2=(x-ty,x-ty)=\|x\|^2-2t(x,y)+t^2\|y\|^2.$$
> Esto es un polinomio cuadrático en $t$ (con coeficiente principal $\|y\|^2>0$) que nunca es negativo. Evaluándolo en el punto donde alcanza su mínimo, $t=\dfrac{(x,y)}{\|y\|^2}$:
> $$0\leq\|x\|^2-2\cdot\frac{(x,y)}{\|y\|^2}\cdot(x,y)+\Big(\frac{(x,y)}{\|y\|^2}\Big)^2\|y\|^2 = \|x\|^2-\frac{2(x,y)^2}{\|y\|^2}+\frac{(x,y)^2}{\|y\|^2} = \|x\|^2-\frac{(x,y)^2}{\|y\|^2}.$$
> Despejando, $(x,y)^2\leq\|x\|^2\|y\|^2$, es decir $|(x,y)|\leq\|x\|\|y\|$.

> [!idea]
> Este es el ingrediente que falta en la prueba de que $\|x\|:=(x,x)^{1/2}$ cumple la desigualdad triangular ([[4.1 - Norma inducida por producto escalar|Problema 4.1]]): sin Cauchy-Schwarz, $\|x+y\|^2=\|x\|^2+2(x,y)+\|y\|^2$ no se puede acotar por $(\|x\|+\|y\|)^2$. Para el producto escalar estándar $(x,y)=x^Ty$ en $\mathbb R^n$, esto da la versión más conocida: $|x^Ty|\leq\|x\|_2\|y\|_2$.

## Ver también
- [[4.1 - Norma inducida por producto escalar|Problema 4.1]] (usa este teorema)
- [[Normas de vectores y matrices]]
- [[Demostraciones - Capitulo 4|Demostraciones]]
