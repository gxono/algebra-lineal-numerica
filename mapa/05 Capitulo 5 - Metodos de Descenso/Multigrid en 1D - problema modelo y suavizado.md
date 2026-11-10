---
tags: [capitulo5, teorema, definicion]
---

# Multigrid en 1D: problema modelo y propiedad de suavizado

> [!idea]
> Dos ideas que, combinadas, dan un precondicionador óptimo: (1) un barrido de [[Metodo de Jacobi|Jacobi]] amortiguado **no resuelve** el sistema, pero **suaviza** el error (elimina las componentes oscilatorias); (2) lo que queda suave en la malla fina se ve **oscilatorio** en una malla más gruesa, donde sí se puede atacar barata y eficazmente. Se puede saltear la deducción (que sigue) y tomar el sistema $(5.16)$ como dato.

## De dónde sale la matriz

Discretizando $-u''(x)=f(x)$, $u(0)=u(1)=0$ por diferencias finitas ($h=1/(n+1)$, diferencia centrada de segundo orden):
$$
A_hu=f, \qquad A_h=\frac1{h^2}\text{tridiag}(-1,2,-1)\in\mathbb{R}^{n\times n}. \tag{5.16}
$$
La misma matriz tridiagonal de [[2.5 - Factorizacion LU tridiagonal|Problema 2.5]] y de los splines del Capítulo 2 - el ejemplo canónico de matriz sparse, sdp y mal condicionada.

## Autovalores y autovectores exactos

$$
\varphi_i^{(k)}=\text{sen}(k\pi x_i), \qquad \lambda_k=\frac4{h^2}\text{sen}^2\Big(\frac{k\pi h}2\Big), \quad k=1,\dots,n.
$$
$k$ chico = modo **suave** (varía poco de nodo a nodo); $k$ grande = modo **oscilatorio**. $\kappa_2(A_h)=O(h^{-2})=O(n^2)$: por la [[Cotas de error del gradiente conjugado|Observación 5.20]], GC sin precondicionar necesita $O(n)$ iteraciones, costo total $O(n^2)$.

> [!demostracion]- Por qué $\varphi^{(k)}$ es autovector de $A_h$ con ese autovalor
> Con $x_i=ih$, $\varphi_i^{(k)}=\text{sen}(k\pi x_i)=\text{sen}(ik\pi h)$, y $\varphi_0^{(k)}=\text{sen}(0)=0$, $\varphi_{n+1}^{(k)}=\text{sen}(k\pi(n+1)h)=\text{sen}(k\pi)=0$ (usando $h=1/(n+1)$) - consistente con las condiciones de borde $u_0=u_{n+1}=0$.
>
> La fila $i$ de $A_h\varphi^{(k)}$ es $\dfrac{1}{h^2}\big(-\varphi_{i-1}^{(k)}+2\varphi_i^{(k)}-\varphi_{i+1}^{(k)}\big)$. Usando la identidad trigonométrica $\text{sen}(\theta-\psi)+\text{sen}(\theta+\psi)=2\,\text{sen}\,\theta\cos\psi$ con $\theta=ik\pi h$, $\psi=k\pi h$:
> $$\varphi_{i-1}^{(k)}+\varphi_{i+1}^{(k)} = \text{sen}\big((i-1)k\pi h\big)+\text{sen}\big((i+1)k\pi h\big) = 2\,\text{sen}(ik\pi h)\cos(k\pi h) = 2\cos(k\pi h)\,\varphi_i^{(k)}.$$
> Por lo tanto,
> $$-\varphi_{i-1}^{(k)}+2\varphi_i^{(k)}-\varphi_{i+1}^{(k)} = 2\varphi_i^{(k)}-2\cos(k\pi h)\varphi_i^{(k)} = 2\big(1-\cos(k\pi h)\big)\varphi_i^{(k)},$$
> así que $(A_h\varphi^{(k)})_i=\dfrac{2}{h^2}\big(1-\cos(k\pi h)\big)\varphi_i^{(k)}$ para todo $i$, es decir $A_h\varphi^{(k)}=\lambda_k\varphi^{(k)}$ con $\lambda_k=\dfrac{2}{h^2}(1-\cos(k\pi h))$. Usando la identidad del ángulo medio $1-\cos\alpha=2\,\text{sen}^2(\alpha/2)$ con $\alpha=k\pi h$:
> $$\lambda_k = \frac{2}{h^2}\cdot2\,\text{sen}^2\Big(\frac{k\pi h}2\Big) = \frac4{h^2}\text{sen}^2\Big(\frac{k\pi h}2\Big).$$
> Como $A_h$ es simétrica (Capítulo 2) y los $\lambda_k$ son todos distintos para $k=1,\dots,n$ (la función $\sin^2$ es inyectiva en $(0,\pi/2)$ y $k\pi h/2\in(0,\pi/2)$ en ese rango), los $\varphi^{(k)}$ son automáticamente ortogonales entre sí y forman una base de $\mathbb R^n$.

## La propiedad de suavizado

Jacobi amortiguado: $u_{m+1}=u_m+\omega D^{-1}(f-A_hu_m)$, $\omega\in(0,1]$ ($\omega=1$ recupera Jacobi del Capítulo 4). El error se amortigua modo a modo: $e_m=\sum_kc_k\mu_k^m\varphi^{(k)}$, con $\mu_k=1-2\omega\,\text{sen}^2(k\pi h/2)$ (mismo cálculo que arriba, ya que $R_\omega=I-\omega\frac{h^2}2A_h$ comparte los autovectores de $A_h$).

$$
\mu_{\text{suav}}(\omega) := \max_{(n+1)/2\leq k\leq n}|\mu_k|, \qquad \text{óptimo: } \omega=\tfrac23,\ \ \mu_{\text{suav}}(2/3)=\tfrac13.
$$

> [!demostracion]- De dónde sale $\omega=2/3$
> Para $k$ en el rango oscilatorio ($\frac{n+1}2\leq k\leq n$), $\theta:=\frac{k\pi h}2\in\big[\frac\pi4,\frac\pi2\big)$ (ya que $k\pi h$ recorre aproximadamente $[\pi/2,\pi)$, usando $h\approx1/n$). Con $s:=\text{sen}^2\theta\in\big[\frac12,1\big)$ (pues $\text{sen}^2$ es creciente en $[0,\pi/2]$), el factor a maximizar es
> $$|\mu_k| = |1-2\omega s|,$$
> una función **lineal** de $s$ (para $\omega$ fijo). El máximo de una función lineal sobre un intervalo se alcanza en uno de los dos extremos, así que
> $$\mu_{\text{suav}}(\omega) = \max_{s\in[1/2,1]}|1-2\omega s| = \max\big\{|1-\omega|,\,|1-2\omega|\big\}.$$
> Para $\omega\in(0,1]$: si $\omega<2/3$, en particular $\omega<1$ y $\omega<1/2$ no necesariamente, pero puede verse que $|1-\omega|>|1-2\omega|$ falla cerca de $\omega$ chico - en cambio, conviene razonar con el cruce directamente. Para $\omega\leq1/2$, ambos $1-\omega\geq1/2>0$ y $1-2\omega\geq0$, con $1-\omega>1-2\omega$: domina $|1-\omega|=1-\omega$, **decreciente** en $\omega$. Para $\omega>1/2$, $1-2\omega<0$, así que $|1-2\omega|=2\omega-1$, **creciente** en $\omega$; y en ese rango $1-\omega$ sigue siendo positivo y decreciente. El máximo de las dos, $\max\{1-\omega,\ 2\omega-1\}$ para $\omega>1/2$ (ambas ya con signo correcto), se minimiza exactamente donde se cruzan las dos rectas (una decreciente, otra creciente):
> $$1-\omega = 2\omega-1 \quad\Longleftrightarrow\quad 2=3\omega \quad\Longleftrightarrow\quad \omega=\frac23.$$
> En ese punto, $\mu_{\text{suav}}(2/3)=1-\tfrac23=\tfrac13$ - y es **independiente de $h$**, porque todo el argumento depende solo de $s\in[1/2,1]$, no de $h$ (la dependencia en $h$ queda enteramente en $\theta=k\pi h/2$, pero el rango $s\in[1/2,1]$ es el mismo para cualquier $h$).

> [!idea]
> Mismo tipo de argumento minimax que la elección del $\alpha$ óptimo de Richardson ([[Lema 4.35]], Capítulo 4): maximizar/minimizar una función lineal sobre un intervalo lleva a comparar los dos extremos, y el óptimo global se alcanza donde esos dos valores de borde se cruzan.

**Con $ω=2/3$, cada barrido reduce *todos* los modos oscilatorios en un factor $1/3$, sin importar $h$.** Los modos suaves ($\mu_k\approx1$) casi no se reducen - por eso Jacoby *solo* (Capítulo 4) es tan lento: está gobernado por esos modos suaves. Esto explica además por qué los métodos estacionarios del Capítulo 4 son lentos en general: $\rho(R_\omega)=\max_k|\mu_k|\approx1$, dominado por los modos suaves.

## La idea de la corrección en malla gruesa

Un modo suave en la malla fina ($\varphi^{(k)}$, $k$ chico respecto de $n$) se ve **oscilatorio** en la malla gruesa (paso $2h$, la mitad de nodos: $n_c\approx n/2$), porque la frecuencia *relativa* $k/n_c$ se duplicó. Trasladando el **error** (no la incógnita directa - el error satisface la *ecuación residual* $A_he=r:=f-A_hu$, con la misma matriz) a esa malla, queda accesible a la relajación, y además la malla gruesa tiene la mitad de incógnitas: más barata.

## Ver también
- [[Multigrid en 1D - transferencia e implementacion]]
- [[Metodo de Jacobi]] (Capítulo 4)
