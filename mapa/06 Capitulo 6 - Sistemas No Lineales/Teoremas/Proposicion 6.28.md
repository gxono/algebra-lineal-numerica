---
tags: [capitulo6, teorema]
---

# Proposición 6.28 (caracterización de funciones convexas)

Sea $f:\mathbb{R}^n\to\mathbb{R}$ diferenciable. Entonces $f$ es convexa si y solo si
$$
f(y)-f(x)\geq\nabla f(x)^T(y-x),\qquad\forall x,y\in\mathbb{R}^n.
$$
Para $n=1$: $f(y)\geq f(x)+f'(x)(y-x)$, el gráfico está por encima de todas sus rectas tangentes.

> [!demostracion]-
> **($\Rightarrow$)** Sean $x,y$ y $0<h<1$. Por convexidad, $f((1-h)x+hy)\leq(1-h)f(x)+hf(y)$. Como $(1-h)x+hy=x+h(y-x)$,
> $$\frac{f(x+h(y-x))-f(x)}{h}\leq f(y)-f(x).$$
> Tomando $h\to0^+$, el lado izquierdo tiende a la derivada direccional $\nabla f(x)^T(y-x)$, y queda $\nabla f(x)^T(y-x)\leq f(y)-f(x)$.
>
> **($\Leftarrow$)** Sean $x,y$, $0<t<1$ y $z:=tx+(1-t)y$. Por hipótesis aplicada en $z$:
> $$f(x)-f(z)\geq\nabla f(z)^T(x-z),\qquad f(y)-f(z)\geq\nabla f(z)^T(y-z).$$
> Multiplicando la primera por $t$, la segunda por $1-t$ y sumando:
> $$tf(x)+(1-t)f(y)-f(z)\geq\nabla f(z)^T\big(tx+(1-t)y-z\big)=\nabla f(z)^T0=0,$$
> es decir $f(tx+(1-t)y)\leq tf(x)+(1-t)f(y)$.

## Ver también
- [[Teorema 6.29]] (la usa dos veces)
- [[Convergencia global de Newton - Newton-Baluev]]
- [[Demostraciones - Capitulo 6|Demostraciones]]
