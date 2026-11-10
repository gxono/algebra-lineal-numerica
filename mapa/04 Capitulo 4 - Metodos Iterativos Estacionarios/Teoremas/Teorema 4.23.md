---
tags: [capitulo4, teorema]
---

# Teorema 4.23

Sea $A\in\mathbb R^{n\times n}$. Las siguientes afirmaciones son equivalentes:
1. $\rho(A)<1$.
2. $\|A\|<1$ para alguna norma inducida.
3. $\lim_{k\to\infty}\|A^k\|=0$, para toda norma inducida $\|\cdot\|$.
4. $\lim_{k\to\infty}A^kx=0$, para todo vector $x\in\mathbb R^n$.

> [!demostracion]- $(1)\Rightarrow(2)\Rightarrow(3)\Rightarrow(4)\Rightarrow(1)$
> **(1)$\Rightarrow$(2):** inmediato por el [[Teorema 4.20]] (tomando $\varepsilon=1-\rho(A)>0$, existe una norma inducida con $\|A\|<\rho(A)+\varepsilon=1$).
>
> **(2)$\Rightarrow$(3):** supongamos $\|A\|_*<1$ para alguna norma inducida $\|\cdot\|_*$, y sea $\|\cdot\|$ **cualquier** otra norma inducida. Como el espacio de matrices $n\times n$ tiene dimensión finita, todas las normas son equivalentes: existe $c>0$ tal que
> $$\|A^k\|\leq c\,\|A^k\|_*, \qquad \forall k\in\mathbb N_0.$$
> Por consistencia/submultiplicatividad de la norma inducida ([[4.4 - Submultiplicatividad norma matricial|Problema 4.4]]), $\|A^k\|_*\leq\|A\|_*^k$, así que
> $$\|A^k\|\leq c\,\|A^k\|_*\leq c\,\|A\|_*^k\to0 \quad\text{cuando } k\to\infty,$$
> porque $\|A\|_*<1$.
>
> **(3)$\Rightarrow$(4):** sea $x\in\mathbb R^n$ y $\|\cdot\|$ una norma en $\mathbb R^n$, con $\|\cdot\|$ también su norma inducida en $\mathbb R^{n\times n}$. Por consistencia,
> $$\|A^kx-0\| = \|A^kx\| \leq \|A^k\|\,\|x\| \to 0 \quad\text{cuando } k\to\infty,$$
> usando (3).
>
> **(4)$\Rightarrow$(1):** sea $\lambda\in\sigma(A)$ (posiblemente complejo) con autovector $x\neq0$ (posiblemente complejo; si $A$ es real y $\lambda$ complejo, se trabaja igual en $\mathbb C^n$ - (4) se extiende a vectores complejos aplicando la hipótesis a la parte real e imaginaria, o directamente usando que $\|A^kx\|=\|\lambda\|^k\|x\|$ tiene sentido en $\mathbb C^n$ con cualquier norma). Entonces $Ax=\lambda x\implies A^kx=\lambda^kx$, y para cualquier norma,
> $$\|A^kx\| = |\lambda|^k\|x\|.$$
> Por (4), $\|A^kx\|\to0$; como $x\neq0$ (así que $\|x\|>0$), esto fuerza $|\lambda|^k\to0$, lo cual ocurre si y solo si $|\lambda|<1$. Como $\lambda\in\sigma(A)$ era arbitrario, $\rho(A)=\max_{\lambda\in\sigma(A)}|\lambda|<1$.

> [!idea]
> Esta cadena es la base de **todo** el capítulo: la convergencia de un método iterativo estacionario $x_{k+1}=Mx_k+c$ depende del error $e_k=M^ke_0$, así que "$e_k\to0$ para todo $e_0$" es exactamente la condición (4) de este teorema con $A=M$ - y por la equivalencia, eso pasa si y solo si $\rho(M)<1$ ([[Teorema de convergencia de metodos iterativos estacionarios]]). La equivalencia (1)$\iff$(2) también es la que permite usar "existe una norma con $\|M\|<1$" como condición *suficiente* más fácil de verificar en la práctica (por ejemplo, vía e.d.d. y $\|\cdot\|_\infty$ para Jacobi/Gauss-Seidel) sin tener que calcular $\rho(M)$ exactamente.

## Ver también
- [[Teorema 4.20]]
- [[Teorema de convergencia de metodos iterativos estacionarios]] (aplicación directa)
- [[Teorema 4.24]] (fórmula de Gelfand, consecuencia de este teorema)
- [[Espectro y radio espectral]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
