---
tags: [capitulo7, algoritmo, teorema]
---

# Algoritmo QR - Hessenberg, shifts y deflación (7.4, segunda parte)

> [!idea]
> El QR "crudo" cuesta $O(n^3)$ por paso y converge linealmente. La versión práctica tiene tres mejoras: (1) **reducir a Hessenberg** una vez ($O(n^3)$); QR conserva la forma y cada paso cuesta $O(n^2)$. (2) **Shift** $\mu_k=(A_k)_{nn}$: la esquina $(n,n)$ hace iteración del cociente de Rayleigh (cúbica en hermíticas). (3) **Deflación**: cuando $a_{n,n-1}^{(k)}\approx0$ se acepta $\lambda_n\approx(A_k)_{nn}$ y se sigue con la submatriz $(n-1)\times(n-1)$. Total $O(n^3)$. El shift de Rayleigh puede **estancarse** (espectro simétrico respecto de $a_{nn}$); el shift de **Wilkinson** lo arregla.

## Por qué Hessenberg

**Matriz upper Hessenberg:** $h_{ij}=0$ si $i>j+1$ (triangular superior + una subdiagonal; ver [[Matriz de Hessenberg]]).
1. Si $A_k$ es Hessenberg, $A_{k+1}=R_kQ_k$ también lo es ([[7.14 - QR preserva Hessenberg|Problema 7.14]]).
2. El paso QR de una Hessenberg cuesta $O(n^2)$ (con [[Rotaciones de Givens|Givens]]: $n-1$ rotaciones, cada una toca 2 filas; ver [[3.7 - Conteo Givens en Hessenberg|Problema 3.7]]).
3. Toda $A$ es unitariamente semejante a una Hessenberg, en $O(n^3)$.

QR conserva la estructura **Hessenberg** pero no otras. Una tridiagonal no simétrica se convierte en Hessenberg llena por arriba; una tridiagonal **simétrica** sigue siendo tridiagonal simétrica (Hessenberg + hermítica = tridiagonal), verificado en el [[7.7 - Implementar el algoritmo QR|Problema 7.7]].

## Por qué no se puede llegar a triangular en finitos pasos

Si un procedimiento finito llevara cualquier matriz a una triangular **semejante**, la diagonal daría los autovalores con un método directo, imposible por [[Autovalores - repaso y por que los metodos son iterativos|Abel]]. Concretamente: para anular **toda** la primera columna con $Q$ y conservar semejanza hay que multiplicar también por $Q^H$ a derecha, y eso mezcla las columnas y destruye los ceros. Una $Q$ que fije $e_1$ ($Qe_1=e_1$) no puede anular las posiciones $2..n$ de la primera columna. En cambio, anular las posiciones **$3..n$** (dejando viva la subdiagonal) sí se puede con una $Q=\begin{pmatrix}1&0^T\\0&\tilde Q\end{pmatrix}$, y multiplicar por $Q^H$ a derecha no toca la primera columna. Por eso la forma alcanzable en finitos pasos es la **Hessenberg**.

## Reducción a Hessenberg con Householder

Paso $j$ ($j=1,\dots,n-2$): $x=A^{(j-1)}_{j+1:n,\,j}$, $u=\frac{x/\|x\|_2\pm e_1}{\sqrt{|x_1|/\|x\|_2+1}}$ ($\|u\|=\sqrt2$), $H_j=I-uu^H$ (con $H_jx=\mp\|x\|_2e_1$), $Q_j=\begin{pmatrix}I_j&0\\0&H_j\end{pmatrix}$, y $A^{(j)}=Q_jA^{(j-1)}Q_j^H$. Tras $n-2$ pasos, $A^{(n-2)}$ es Hessenberg y semejante a $A$.

> [!algoritmo] Reducción a Hessenberg (verificada: entradas bajo la subdiagonal $\sim10^{-16}$, autovalores iguales a $10^{-14}$)
> ```julia
> function hessenberg_householder(A)
>     A = copy(A); n = size(A, 1)
>     for j in 1:n-2
>         x = A[j+1:n, j]
>         nx = norm(x); nx == 0 && continue
>         s = x[1] == 0 ? 1.0 : sign(x[1])           # signo que evita cancelación
>         u = (x/nx + s*[1; zeros(n-j-1)]) / sqrt(abs(x[1])/nx + 1)
>         A[j+1:n, :] -= u * (u' * A[j+1:n, :])     # Q_j A     (solo filas j+1..n)
>         A[:, j+1:n] -= (A[:, j+1:n] * u) * u'     # (Q_j A) Q_jᴴ (solo columnas j+1..n)
>     end
>     return A
> end
> ```
> Usa las [[Matrices de Householder]] del Capítulo 3. En Julia: `hessenberg(A)` (LAPACK), con `.H` y `.Q`.

Costo $O(n^3)$: cada paso aplica un reflector a $O(n)$ filas y $O(n)$ columnas de longitud $O(n)$.

## QR con shifts

Dada una sucesión $\mu_0,\mu_1,\dots$:
$$
A_k-\mu_kI=Q_kR_k,\qquad A_{k+1}=R_kQ_k+\mu_kI=Q_k^HA_kQ_k.
$$
Sigue siendo semejanza unitaria para cualquier shift, y el Teorema 7.10(b) vale igual.

**Proposición 7.14.** $(A-\mu_kI)(A-\mu_{k-1}I)\cdots(A-\mu_0I)=P_kU_k$. Ver [[Proposicion 7.14]].

Repitiendo la cuenta de la Obs. 7.12, la última columna de $P_k$ es **potencias inversas con shifts** aplicadas a $A^H$: en el paso $k$ el error se contrae en $\approx\max_{i\neq n}\left|\frac{\lambda_n-\mu_k}{\lambda_i-\mu_k}\right|$, diminuto si $\mu_k\approx\lambda_n$. Y $(A_k)_{nn}$ **ya es** el cociente de Rayleigh de la última columna de $P_{k-1}$: tomar
$$
\mu_k=(A_k)_{nn}\qquad\text{(shift de Rayleigh)}
$$
es la iteración del cociente de Rayleigh (cúbica en hermíticas) **sin pagar refactorización**.

## Deflación

Cuando $|a^{(k)}_{n,n-1}|<\text{tol}$: aceptar $\lambda_n\approx(A_k)_{nn}$, descartar última fila y columna, seguir con la $(n-1)\times(n-1)$. Unas pocas iteraciones $O(n^2)$ por autovalor ⟹ **$O(n^3)$ en total**.

> [!algoritmo] QR con shift y deflación ([[7.9 - QR con shift y deflacion|Problema 7.9]])
> ```julia
> function qr_shift_deflacion(A; tol=1e-12, kmax=1000, shift=:rayleigh)
>     Ak = copy(A); n = size(A, 1); λ = zeros(n); its = Int[]
>     m = n
>     while m > 1
>         B = Ak[1:m, 1:m]; it = 0
>         while abs(B[m, m-1]) ≥ tol && it < kmax
>             it += 1
>             μ = shift == :rayleigh ? B[m, m] : wilkinson(B[m-1:m, m-1:m])
>             F = qr(B - μ*I)
>             B = F.R * Matrix(F.Q) + μ*I
>         end
>         λ[m] = B[m, m]; push!(its, it)
>         Ak[1:m, 1:m] = B
>         m -= 1                                  # deflación
>     end
>     λ[1] = Ak[1, 1]
>     return λ, its
> end
>
> # autovalor del bloque 2x2 final más cercano a su esquina (n,n)
> function wilkinson(B)
>     a, b, c, d = B[1,1], B[1,2], B[2,1], B[2,2]
>     t = (a + d)/2; disc = sqrt(complex(t^2 - (a*d - b*c)))
>     e1, e2 = t + disc, t - disc
>     real(abs(e1 - d) ≤ abs(e2 - d) ? e1 : e2)
> end
> ```

## Cuándo falla el shift de Rayleigh

Ejemplo del apunte: $A=\begin{pmatrix}0&1\\1&0\end{pmatrix}$, autovalores $\pm1$ (mismo módulo). $a_{22}=0$ ⟹ $\mu=0$ ⟹ QR de $A$ es $Q=A$, $R=I$ (salvo signos) ⟹ $A_1=RQ=A$: **quieto para siempre**.

> [!pendiente] Lo mismo pasa con $A_N=\text{tridiag}(-1,2,-1)$, la matriz de los problemas
> La diagonal de $A_N$ es constante $=2$ y su espectro es **simétrico respecto de 2** ($\lambda_i+\lambda_{N+1-i}=4$). $T=A_N-2I$ tiene diagonal nula, y el paso QR **conserva la diagonal nula** (se prueba con la unicidad del QR, [[7.12 - Unicidad de QR salvo fases|Problema 7.12]]). Entonces $(A_k)_{nn}=2$ siempre: el shift de Rayleigh queda clavado en $\mu=2$ y lo que se hace es QR **sin shift** sobre $T$, cuyos autovalores vienen en pares $\pm$ de **igual módulo**: no converge. En aritmética exacta se estanca para siempre; en punto flotante, el redondeo rompe la simetría tras un número impredecible de pasos (41 con $N=5$; más de 500 con $N=15$ y $N=50$, con autovalores **incorrectos**). Ver [[7.9 - QR con shift y deflacion|Problema 7.9]].

**Shift de Wilkinson:** $\mu_k$ = el autovalor del bloque $2\times2$ inferior derecho $\begin{pmatrix}a_{n-1,n-1}&a_{n-1,n}\\a_{n,n-1}&a_{nn}\end{pmatrix}$ más cercano a $a_{nn}$. En el ejemplo da $\mu_0=1$ (empate, se elige uno) y converge en **un** paso: $A-I$ es singular y $A_1=\text{diag}(-1,1)$. En $A_N$: 8, 33 y 102 iteraciones totales para $N=5,15,50$, contra 64, 412 y 3401 sin shift.

Si $A$ es real y el bloque tiene autovalores complejos conjugados, el shift de Wilkinson es complejo; para no salir de la aritmética real se usa el **doble shift** (aplica el par conjugado en un paso).

## Ver también
- [[Algoritmo QR - definicion y convergencia]]
- [[Metodo de las potencias inversas]]
- [[Matriz de Hessenberg]], [[Matrices de Householder]], [[Rotaciones de Givens]]
- [[GMRes - implementacion y rotaciones de Givens]] (otra QR incremental de una Hessenberg, Cap. 5)
- [[Tablas Comparativas - Metodos para autovalores]]
- [[Capitulo 7 - MOC]]
