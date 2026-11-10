---
tags: [capitulo6, definicion, teorema]
---

# Contracciones y teorema de punto fijo (Contraction Mapping Theorem)

> [!idea]
> Si $G$ manda un cerrado $D$ en sí mismo y **acerca** puntos ($\|G(x)-G(y)\|\leq\gamma\|x-y\|$, $\gamma<1$), entonces tiene **un único** punto fijo en $D$ y la iteración converge desde **cualquier** $x_0\in D$, linealmente con factor $\gamma$. A diferencia de [[Iteracion de punto fijo y teorema de Ostrowski|Ostrowski]], no hace falta saber de antemano que existe $x_*$ ni que $G$ sea diferenciable: el teorema **demuestra** la existencia.

## Definiciones 6.10 y 6.11

- **Lipschitz** en $D$ con constante $\gamma>0$: $\|G(x)-G(y)\|\leq\gamma\|x-y\|$ para todo $x,y\in D$.
- **Contracción** en $D$ (con la norma $\|\cdot\|$): Lipschitz con $0<\gamma<1$.

Ser contracción **depende de la norma**: $G$ puede ser contracción en $\|\cdot\|_\infty$ y no en $\|\cdot\|_2$. Alcanza con encontrar **una** norma en la que lo sea.

## Teorema 6.12

> Sea $D\subset\mathbb{R}^n$ **cerrado** y $G:D\to D$ una contracción con constante $\gamma\in(0,1)$. Entonces $G$ tiene un único punto fijo $x_*\in D$, que es punto de atracción, y para todo $x_0\in D$:
> - convergencia lineal: $\|x_{k+1}-x_*\|\leq\gamma\|x_k-x_*\|$;
> - **estimación del error a priori**: $\|x_k-x_*\|\leq\dfrac{\gamma^k}{1-\gamma}\|x_1-x_0\|$.

Demostración completa: [[Teorema 6.12]] (la sucesión es de Cauchy; $D$ cerrado da $x_*\in D$; continuidad de $G$ da $x_*=G(x_*)$; contracción da unicidad).

> [!idea] Las tres hipótesis se usan
> - **$D$ cerrado:** para que el límite de la sucesión de Cauchy quede en $D$.
> - **$G(D)\subset D$:** para que la iteración esté bien definida (es lo que más se olvida verificar en los ejercicios).
> - **$\gamma<1$ uniforme:** "$\|G(x)-G(y)\|<\|x-y\|$" sin una constante uniforme **no alcanza**: $G(x)=\sqrt{1+x^2}$ en $\mathbb{R}$ cumple eso y no tiene punto fijo ([[6.3 - Contraccion en acotados|Problema 6.3]]).

**Uso práctico de la cota a priori:** para garantizar $\|x_k-x_*\|\leq\text{tol}$ alcanza con
$$
k\geq\frac{\log\big(\text{tol}\,(1-\gamma)/\|x_1-x_0\|\big)}{\log\gamma}.
$$
En el [[6.4 - Region para el teorema de contraccion|Problema 6.4]] esto da $k=44$ para $\text{tol}=10^{-8}$; en la práctica alcanza con muchas menos, porque la tasa real ($\approx0.3$) es mejor que la cota $\gamma\approx0.65$.

## Observación 6.13 (el caso lineal, otra vez)

Para $G(x)=Mx+c$ en $D=\mathbb{R}^n$: $G$ es contracción en **alguna** norma $\iff\rho(M)<1$. Así el Teorema 6.12 generaliza el [[Teorema de convergencia de metodos iterativos estacionarios|teorema de convergencia del Capítulo 4]]. Demostración: [[Observacion 6.13]].

## Cómo probar que algo es contracción

Casi nunca con la definición directa: se acota el jacobiano y se usa el [[Criterio de contraccion via la derivada|Corolario 6.19]] (si $D$ es convexo y $\|G'(x)\|\leq\alpha<1$ en $D$, $G$ es contracción con constante $\alpha$). La receta completa está en [[Criterio de contraccion via la derivada#Receta para los ejercicios|esa nota]].

## Ver también
- [[Iteracion de punto fijo y teorema de Ostrowski]]
- [[Criterio de contraccion via la derivada]]
- [[6.4 - Region para el teorema de contraccion]], [[6.5 - Contraccion para h chico]] (aplicaciones)
- [[Capitulo 6 - MOC]]
