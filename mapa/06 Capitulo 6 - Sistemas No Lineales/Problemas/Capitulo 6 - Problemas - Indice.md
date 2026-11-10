---
tags: [capitulo6, moc]
---

# Capítulo 6 - Índice de problemas

| # | Tema |
|---|---|
| [[6.1 - Rho igual a 1 no decide\|6.1]] | $G'(0)=1$: $x-x^3$ atrae, $x+x^3$ no |
| [[6.2 - Contraccion local si rho menor que 1\|6.2]] | $\rho(G'(x_*))<1$ ⟹ contracción en una bola cerrada |
| [[6.3 - Contraccion en acotados\|6.3]] | Contracción en cada acotado: unicidad y convergencia sí, existencia no ($\sqrt{1+x^2}$) |
| [[6.4 - Region para el teorema de contraccion\|6.4]] | $D=[0,\frac12]^2$, $\gamma=\frac{9}{8\sqrt3}\approx0.65$ |
| [[6.5 - Contraccion para h chico\|6.5]] | $\Vert G'\Vert _\infty\leq2.28h$: contracción si $h<0.438$ |
| [[6.6 - No hay teorema del valor medio vectorial\|6.6]] | $G=(x_1^3,x_1^2)$: no existe $z$ intermedio |
| [[6.7 - Demostrar Teorema 6.17\|6.7]] | $\Vert \int H\Vert \leq\int\Vert H\Vert $ por sumas de Riemann |
| [[6.8 - Derivada de G en Newton es nula\|6.8]] | Newton como punto fijo: $G'(x_*)=0$ |
| [[6.9 - Programar Newton\|6.9]] | `newton(F, DF, x0, miter, tol)` |
| [[6.10 - Orden experimental de Newton\|6.10]] | Orden experimental 2 (con `BigFloat`) |
| [[6.11 - Cuatro ceros y fminsearch\|6.11]] | Curvas de nivel, Newton, Nelder–Mead (`Optim.jl`) |
| [[6.12 - Todas las raices de dos sistemas\|6.12]] | (a) 2 raíces, (b) 4 raíces, con valores exactos |
| [[6.13 - Minimo con Newton\|6.13]] | Newton sobre $\nabla f=0$; ciclo atractor de período 5 |
| [[6.14 - Newton para x al cuadrado\|6.14]] | $F'(x_*)=0$ ⟹ lineal (tasa $\frac12$) |
| [[6.15 - Newton para x mas x a la 1 mas alfa\|6.15]] | $F'$ no Lipschitz ⟹ orden $1+\alpha$ |
| [[6.16 - Newton simplificado\|6.16]] | $F'(x_0)$ fija: converge localmente, lineal |
| [[6.17 - Newton en problema discretizado\|6.17]] | $Ax+h^2x^3$: convergencia local, orden **3** |
| [[6.18 - Newton-Baluev sin convexidad\|6.18]] | $\arctan x$ diverge si $\vert x_0\vert >1.39$; $2x+\sin x$: 2-ciclo $\pm\pi$ superatractor |

Volver a [[Capitulo 6 - MOC]].
