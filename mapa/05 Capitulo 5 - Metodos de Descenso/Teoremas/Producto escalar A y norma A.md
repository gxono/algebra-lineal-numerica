---
tags: [capitulo5, teorema]
---

# El producto escalar $A$ y la norma $A$

## Definición

Para $A$ sdp, $(x,y)_A := x^TAy$.

## Teorema

Si $A$ es sdp, $(x,y)_A:=x^TAy$ es un producto escalar en $\mathbb{R}^n$, y por lo tanto $\|z\|_A:=(z,z)_A^{1/2}=\sqrt{z^TAz}$ es una norma (por [[4.1 - Norma inducida por producto escalar|Problema 4.1]], todo producto escalar induce una norma).

> [!demostracion]-
> Hay que verificar los tres axiomas de producto escalar.
>
> **Bilinealidad.** Para $\alpha,\beta\in\mathbb{R}$ y $x,y,z\in\mathbb{R}^n$: $(\alpha x+\beta y,z)_A=(\alpha x+\beta y)^TAz=\alpha x^TAz+\beta y^TAz=\alpha(x,z)_A+\beta(y,z)_A$ - lineal en la primera componente porque la transposición y el producto matricial son lineales. Por simetría (ver el punto siguiente) también es lineal en la segunda.
>
> **Simetría.** $(x,y)_A=x^TAy$ es un escalar, así que coincide con su propia transpuesta: $(x,y)_A=(x^TAy)^T=y^TA^Tx=y^TAx=(y,x)_A$, usando $A^T=A$ ($A$ simétrica).
>
> **Definida positiva.** $(x,x)_A=x^TAx>0$ para todo $x\neq0$ - exactamente la definición de que $A$ sea sdp. Y $(0,0)_A=0$.
>
> Los tres axiomas se cumplen, así que $(\cdot,\cdot)_A$ es un producto escalar.

**Por qué importa:** esta es la base silenciosa de **todo** el Capítulo 5 - cada vez que aparece $\|z\|_A$ (en el [[Funcional de energia y caracterizacion|funcional de energía]], en las cotas de [[Cotas de error del gradiente conjugado|GC]], en todas las demostraciones de convergencia), se está usando implícitamente que esto es una norma de verdad (no negativa, homogénea, con desigualdad triangular) - sin este resultado, expresiones como "$x_k$ minimiza $\|\cdot-x_*\|_A$" no tendrían sentido como problema de distancia.

## Ver también
- [[Funcional de energia y caracterizacion]]
- [[4.1 - Norma inducida por producto escalar|Problema 4.1]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
