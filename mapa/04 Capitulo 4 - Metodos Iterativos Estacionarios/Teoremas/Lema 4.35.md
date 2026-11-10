---
tags: [capitulo4, teorema]
---

# Lema 4.35 ($\alpha$ óptimo de Richardson)

El [[Metodo de Richardson|método de Richardson]] converge para $A$ sdp, a partir de cualquier $x_0\in\mathbb R^n$, si y solo si
$$
0<\alpha<\frac{2}{\lambda_{\max}}.
$$
Además, el valor de $\alpha$ que minimiza $\rho(I-\alpha A)$ es
$$
\alpha_{\text{opt}} = \frac{2}{\lambda_{\min}+\lambda_{\max}}.
$$

> [!demostracion]-
> **Condición de convergencia.** Como $A$ es sdp, es diagonalizable ortogonalmente con autovalores reales $0<\lambda_{\min}\leq\cdots\leq\lambda_{\max}$, y $M_\alpha:=I-\alpha A$ tiene los mismos autovectores que $A$ con autovalores $1-\alpha\lambda_i$. Por el [[Teorema de convergencia de metodos iterativos estacionarios|teorema de convergencia de métodos estacionarios]], Richardson converge para todo $x_0$ si y solo si $\rho(M_\alpha)<1$, es decir
> $$|1-\alpha\lambda_i|<1, \qquad \forall i.$$
> Para $\alpha=0$ esto falla ($|1|=1$). Para $\alpha\neq0$, $|1-\alpha\lambda_i|<1\iff0<\alpha\lambda_i<2\iff0<\alpha<2/\lambda_i$ (usando $\lambda_i>0$; si $\alpha<0$, $\alpha\lambda_i<0$ y la desigualdad $0<\alpha\lambda_i$ falla directamente). Esto debe valer para **todos** los $\lambda_i$ simultáneamente, y la condición más restrictiva es la del $\lambda_i$ más grande: $0<\alpha<2/\lambda_{\max}$.
>
> **$\alpha$ óptimo.** Hay que minimizar, sobre $\alpha$,
> $$\rho(M_\alpha) = \max_i|1-\alpha\lambda_i| = \max\{|1-\alpha\lambda_{\min}|,\ |1-\alpha\lambda_{\max}|\}$$
> (el máximo sobre $i$ de $|1-\alpha\lambda_i|$ se alcanza en uno de los dos extremos del espectro, porque $\alpha\mapsto|1-\alpha\lambda|$ es una función convexa de $\lambda$ para $\alpha$ fijo - de hecho lineal en $\lambda$ a cada lado de $\lambda=1/\alpha$, así que entre $\lambda_{\min}$ y $\lambda_{\max}$ el máximo de $|1-\alpha\lambda|$ está en uno de los extremos).
>
> Como función de $\alpha$: $g_1(\alpha):=|1-\alpha\lambda_{\min}|$ es decreciente-luego-creciente con mínimo (valor $0$) en $\alpha=1/\lambda_{\min}$, y $g_2(\alpha):=|1-\alpha\lambda_{\max}|$ tiene su mínimo (valor $0$) en $\alpha=1/\lambda_{\max}<1/\lambda_{\min}$. Para $\alpha$ pequeño ($\alpha<1/\lambda_{\max}$), ambos $1-\alpha\lambda_{\min}$ y $1-\alpha\lambda_{\max}$ son positivos, así que $g_2(\alpha)>g_1(\alpha)$ (restar más de $1$ con $\lambda_{\max}$ más grande); para $\alpha$ grande ($\alpha>1/\lambda_{\min}$), ambos son negativos y $g_1(\alpha)>g_2(\alpha)$. Por continuidad, $\max\{g_1,g_2\}$ se minimiza exactamente donde $g_1=g_2$, es decir donde las dos rectas se cruzan (una creciendo, otra decreciendo, en la región donde $1-\alpha\lambda_{\min}\geq0\geq1-\alpha\lambda_{\max}$ no ocurre para $\lambda_{\min}<\lambda_{\max}$ salvo en el cruce):
> $$\alpha\lambda_{\max}-1 = 1-\alpha\lambda_{\min} \quad\Longleftrightarrow\quad \alpha(\lambda_{\min}+\lambda_{\max})=2 \quad\Longleftrightarrow\quad \alpha_{\text{opt}}=\frac{2}{\lambda_{\min}+\lambda_{\max}}.$$
> En ese punto, $\rho(M_{\alpha_{\text{opt}}})=1-\alpha_{\text{opt}}\lambda_{\min}=1-\dfrac{2\lambda_{\min}}{\lambda_{\min}+\lambda_{\max}}=\dfrac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\max}+\lambda_{\min}}=\dfrac{\kappa_2(A)-1}{\kappa_2(A)+1}$ (dividiendo numerador y denominador por $\lambda_{\min}$, y usando $\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$ del [[Lema 4.14]]).

> [!idea]
> El argumento geométrico (Figura: dos rectas $|1-\alpha\lambda_{\min}|$ y $|1-\alpha\lambda_{\max}|$ como función de $\alpha$, el máximo de ambas se minimiza en el cruce) es el mismo tipo de razonamiento "minimax" que reaparece, mucho más elaborado, en la elección del polinomio de Chebyshev para las cotas de GC ([[Elegir el polinomio para las cotas de GC]]) - ahí también se busca minimizar el máximo de algo sobre un intervalo, y la solución óptima también "equioscila" entre los extremos.

## Ver también
- [[Metodo de Richardson]] (nota de tema completa)
- [[Teorema 4.36]], [[Teorema 4.37]] (cotas de convergencia con este $\alpha_{\text{opt}}$)
- [[Demostraciones - Capitulo 4|Demostraciones]]
