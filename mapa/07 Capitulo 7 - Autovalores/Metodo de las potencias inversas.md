---
tags: [capitulo7, algoritmo, teorema, observacion]
---

# Método de las potencias inversas (7.3)

> [!idea]
> Aplicar el método de las potencias a $B=(A-\mu I)^{-1}$: sus autovalores son $\frac{1}{\lambda_i-\mu}$, y el dominante corresponde al $\lambda_i$ **más cercano al shift $\mu$**. Tasa $\rho_\mu=\max_{i\neq1}\frac{|\lambda_1-\mu|}{|\lambda_i-\mu|}$: **la elegimos nosotros** a través de $\mu$. Se factoriza $A-\mu I$ una vez y cada paso son dos triangulares. Si además se actualiza $\mu$ con el cociente de Rayleigh, la convergencia pasa a ser **cúbica** para hermíticas (iteración del cociente de Rayleigh).

## La transformación

$\sigma(A-\mu I)=\sigma(A)-\mu$ y $\sigma\big((A-\mu I)^{-1}\big)=\big\{\frac{1}{\lambda-\mu}:\lambda\in\sigma(A)\big\}$, **con los mismos autovectores**. Si $\mu\to\lambda_1$, $\frac{1}{\lambda_1-\mu}\to\infty$ y los demás quedan acotados: $\lambda_1$ pasa a ser dominante para $B$.

## Proposición 7.7

> $\lambda_1$ autovalor simple de $A$ con autovectores $x_1$ ($\|x_1\|_2=1$) e $y_1^H$; $\mu\notin\sigma(A)$ con $|\lambda_1-\mu|<|\lambda_i-\mu|$ para todo $\lambda_i\neq\lambda_1$. Si $y_1^Hu_0\neq0$, las potencias inversas convergen en dirección a $x_1$ con tasa $(\rho_\mu+\varepsilon)^k$,
> $$\rho_\mu:=\max_{\lambda_i\neq\lambda_1}\left|\frac{\lambda_1-\mu}{\lambda_i-\mu}\right|,$$
> y $\mu_k\to\lambda_1$.

Es la [[Proposicion 7.2]] aplicada a $B$. Demostración: [[Proposicion 7.7]].

**Lectura del apunte:** $\rho_\mu$ es **proporcional a $|\lambda_1-\mu|$**. Con $|\lambda_1-\mu|\approx10^{-3}$ y los demás autovalores a distancia $O(1)$, $\rho_\mu\approx10^{-3}$: tres dígitos por iteración. La calidad de $\mu$ afecta la **velocidad**, no el límite (mientras $\lambda_1$ siga siendo el más cercano).

> [!algoritmo] Potencias inversas con shift fijo $\mu$
> ```julia
> function potencias_inversas(A, u0, μ; tol=1e-10, kmax=10_000)
>     F = lu(A - μ*I)                   # se factoriza UNA vez: O(n³)
>     u = u0 / norm(u0); v = F \ u
>     k = 0
>     while true
>         k += 1
>         u = v / norm(v)
>         v = F \ u                     # v_{k+1} = (A-μI)⁻¹ u_k : dos triangulares, O(n²)
>         θ = dot(u, v)                 # cociente de Rayleigh de B = (A-μI)⁻¹
>         λ = 1/θ + μ                   # μ_k del apunte: autovalor de A
>         (norm(v - θ*u) ≤ tol*abs(θ) || k == kmax) && return λ, u, k
>     end
> end
> ```

## Observación 7.8: el mal condicionamiento juega a favor

Si $\mu\approx\lambda_1$, $A-\mu I$ es casi singular y $\kappa_2(A-\mu I)$ es enorme: el vector calculado $\hat v$ tiene un error hacia adelante grande. Pero:
- el solver es **estable hacia atrás** ([[Estabilidad y numero de condicion|Sección 3.5]]): $\hat v$ es solución exacta de $(A-\mu I+E)\hat v=u_k$ con $\|E\|_2\approx\varepsilon_{\text{maq}}\|A\|_2$;
- el error se amplifica justo en la dirección donde $(A-\mu I)^{-1}$ es grande, **que es la de $x_1$**: lo que se agranda es la componente buscada;
- el paso siguiente **normaliza**, y el factor espurio desaparece.

Si $A-\mu I$ resulta numéricamente singular, alcanza con perturbar $\mu$ un poco.

**Por qué no actualizar $\mu=\mu_k$ en cada paso** (en general): habría que refactorizar $A-\mu_kI$ cada vez ($O(n^3)$ por paso), y se pierde el control sobre a qué autovalor converge. Si se acepta ese costo, $\mu_k$ = cociente de Rayleigh da la **iteración del cociente de Rayleigh**, con convergencia **cúbica** para hermíticas.

> [!algoritmo] Iteración del cociente de Rayleigh (shift actualizado, [[7.5 - Potencias inversas con shift actualizado|Problema 7.5]])
> ```julia
> function rayleigh_iter(A, u0; μ0=0.0, tol=1e-10, kmax=100)
>     u = u0 / norm(u0); μ = μ0; k = 0
>     while k < kmax
>         k += 1
>         v = (A - μ*I) \ u             # refactoriza en cada paso
>         u = v / norm(v)
>         μ = dot(u, A*u)               # nuevo shift = cociente de Rayleigh
>         norm(A*u - μ*u) ≤ tol*abs(μ) && break
>     end
>     return μ, u, k
> end
> ```
> Orden medido ([[7.5 - Potencias inversas con shift actualizado|7.5]], en `BigFloat`): **3** para $A$ simétrica, **2** para no simétrica. Errores típicos (simétrica): $8\times10^{-4}\to6\times10^{-10}\to2\times10^{-28}\to5\times10^{-84}$.

> [!pendiente] Con shift actualizado no se elige a qué autovalor converge
> En $A_N$ con $\mu_0=0$ (el más cercano es $\lambda_1$), la iteración con shift actualizado convergió a $\lambda_1$, $\lambda_3$, $\lambda_4$ o $\lambda_2$ según $N$ ([[7.6 - Potencias inversas con y sin shift en A_N|Problema 7.6]]): muy rápido (4–7 iteraciones), pero **no necesariamente al autovalor buscado**. Con shift fijo se garantiza el más cercano a $\mu$.

## Observación 7.9 (ventajas y desventajas)

- **Ventajas:** convergencia muy rápida y se puede obtener **cualquier** autovalor (no solo el dominante), si se tiene una buena aproximación de cada uno.
- **Desventajas:** hace falta una buena aproximación del autovalor y resolver un sistema por paso. Pero:
  - matriz llena/chica: se factoriza una vez y cada paso cuesta $O(n^2)$;
  - matriz grande y rala: método iterativo (Cap. 5) para cada sistema, arrancando de $u_k$ del paso anterior.

## Potencias inversas en $A_N$: la tasa no depende de $N$

Con $\mu=0$ (autovalor mínimo), $\rho_0=\frac{\lambda_1}{\lambda_2}=\frac{1}{4\cos^2\frac{\pi}{2(N+1)}}\to\frac14$: **≈17 iteraciones** para $10^{-10}$ con $N=5,15,50,100$, contra 96, 641 y 5756 del método de las potencias ([[7.4 - Potencias inversas en la matriz A_N|Problema 7.4]]).

## Ver también
- [[Metodo de las potencias]]
- [[Cociente de Rayleigh y criterio de parada]]
- [[Algoritmo QR - Hessenberg, shifts y deflacion]] (la esquina $(n,n)$ de QR **es** una iteración inversa)
- [[Tablas Comparativas - Metodos para autovalores]]
- [[Capitulo 7 - MOC]]
