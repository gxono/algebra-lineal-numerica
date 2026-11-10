---
tags: [capitulo7, algoritmo, teorema, lema, observacion]
---

# Método de las potencias (7.2)

> [!idea]
> Multiplicar repetidamente por $A$ (y normalizar) amplifica la componente en la dirección del autovalor de **módulo máximo**: $u_k=A^ku_0/\|A^ku_0\|$ converge **en dirección** a $x_1$ con tasa $|\lambda_2/\lambda_1|$. El autovalor se estima con el **cociente de Rayleigh** $\mu_k=u_k^HAu_k$. Es simple y barato (1 producto matriz-vector por paso) pero solo da el autovalor dominante, y es lento si $|\lambda_2|\approx|\lambda_1|$.

## La recursión (7.1)

$$
u_k=\sigma_{k-1}Au_{k-1},\qquad\sigma_{k-1}=\frac{1}{\|Au_{k-1}\|},\qquad k=1,2,\dots
$$
Los $\sigma$ solo evitan overflow/underflow.

**Lema 7.1.** Si $\|u_0\|=1$ y $A^ku_0\neq0$ para todo $k$, la recursión está bien definida y $u_k=\dfrac{A^ku_0}{\|A^ku_0\|}$. Normalizar en cada paso o solo al final da lo mismo. Ver [[Lema 7.1]].

## Convergencia: Proposición 7.2

> Si $|\lambda_1|>|\lambda_2|\geq\cdots\geq|\lambda_n|$ ($\lambda_1$ **dominante**, necesariamente simple), $x_1$ autovector a derecha con $\|x_1\|_2=1$, $y_1^H$ autovector a izquierda, y **$y_1^Hu_0\neq0$**, entonces existen $\omega_k\in\mathbb{C}$, $|\omega_k|=1$, con
> $$\omega_ku_k\to x_1,$$
> con tasa asintótica $(|\lambda_2/\lambda_1|+\varepsilon)^k$ para todo $\varepsilon>0$.

Idea: con $X^{-1}AX=\begin{pmatrix}\lambda_1&0^T\\0&B\end{pmatrix}$ y $X^{-1}u_0=\begin{pmatrix}\nu_0\\w_0\end{pmatrix}$,
$$
A^ku_0=\lambda_1^k\big(\nu_0x_1+z_k\big),\qquad z_k=X_2(B/\lambda_1)^kw_0\to0\quad(\rho(B/\lambda_1)=|\lambda_2/\lambda_1|<1).
$$
Demostración completa: [[Proposicion 7.2]].

## Observación 7.3 (comentarios)

- "Converge en dirección" = $\text{dist}(u_k,\{\text{autovectores de }\lambda_1\})\to0$. El signo/fase $\omega_k$ puede oscilar (por ejemplo si $\lambda_1<0$, $u_k$ alterna de signo).
- La tasa $(|\lambda_2/\lambda_1|+\varepsilon)^k$ sale de $\|(B/\lambda_1)^k\|\leq C(\rho(B/\lambda_1)+\varepsilon)^k$ (fórmula de Gelfand, [[Teorema 4.24]]).
- **$y_1^Hu_0\neq0$** no se puede verificar (no conocemos $y_1$). Se toma $u_0$ **al azar**. Además no alcanza con que sea no nulo: si es chico relativo a $\|u_0\|$, hacen falta muchas iteraciones.
- Ventaja: cada paso es 1 producto matriz-vector + normalización (sirve para matrices enormes y ralas).
- Desventajas: solo garantiza el autovalor **dominante**, que puede **no existir** (matriz real con $\lambda_1$ complejo: $\bar\lambda_1$ tiene el mismo módulo). Si $|\lambda_2|\approx|\lambda_1|$, es lentísimo.
- **Conclusión del apunte:** herramienta de propósito especial (autovalor dominante con $|\lambda_2|/|\lambda_1|$ chico), no un algoritmo general.

## El algoritmo

Con el [[Cociente de Rayleigh y criterio de parada|cociente de Rayleigh]] como estimador del autovalor, y reutilizando $v_{k+1}=Au_k$ para tener un solo producto por paso:

> [!algoritmo] Método de las potencias
> ```julia
> function potencias(A, u0; tol=1e-10, kmax=10_000)
>     u = u0 / norm(u0)
>     v = A*u
>     k = 0
>     while true
>         k += 1
>         u = v / norm(v)            # u_k
>         v = A*u                    # v_{k+1} = A u_k (el único producto del paso)
>         μ = dot(u, v)              # μ_k = u_kᴴ A u_k   (dot conjuga el primero)
>         r = v - μ*u                # residuo r_k = A u_k - μ_k u_k, gratis
>         (norm(r) ≤ tol*abs(μ) || k == kmax) && return μ, u, k
>     end
> end
> ```

## Qué tasa se ve en la práctica

| | error del autovector | error del autovalor $\mu_k$ |
|---|---|---|
| $A$ cualquiera | $\sim\vert \lambda_2/\lambda_1\vert ^k$ | $\sim\vert \lambda_2/\lambda_1\vert ^k$ |
| $A$ **hermítica** | $\sim\vert \lambda_2/\lambda_1\vert ^k$ | $\sim\vert \lambda_2/\lambda_1\vert ^{2k}$ (**el doble de rápido**) |

El cuadrado en el caso hermítico viene de que el cociente de Rayleigh tiene error **cuadrático** en el error del autovector ([[Observacion 7.4]]). Verificado en el [[7.1 - Programar el metodo de las potencias|Problema 7.1]] (tasas $0.8$ y $0.64$ con $\lambda_2/\lambda_1=0.8$).

> [!pendiente] El método de las potencias en $A_N=\text{tridiag}(-1,2,-1)$ es pésimo
> $|\lambda_{N-1}/\lambda_N|=\frac{\cos^2(\pi/(N+1))}{\cos^2(\pi/(2(N+1)))}\approx1-\frac{3\pi^2}{4(N+1)^2}$: con $N=50$ la tasa es $0.99716$ y hacen falta **5756** iteraciones para $10^{-10}$ ([[7.2 - Potencias en la matriz A_N|Problema 7.2]]). Las iteraciones crecen como $O(N^2)$.

## Ver también
- [[Cociente de Rayleigh y criterio de parada]]
- [[Metodo de las potencias inversas]] (la mejora: shift + inversión)
- [[Algoritmo QR - definicion y convergencia]] (la primera columna de QR **es** el método de las potencias)
- [[Capitulo 7 - MOC]]
