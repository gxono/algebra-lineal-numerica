---
tags: [capitulo7, algoritmo, teorema, observacion]
---

# Algoritmo QR - definición y convergencia (7.4, primera parte)

> [!idea]
> Factorizar $A_k=Q_kR_k$ y multiplicar al revés: $A_{k+1}=R_kQ_k=Q_k^HA_kQ_k$. Todas las $A_k$ son **unitariamente semejantes** a $A$, y (bajo hipótesis) convergen a una forma de **Schur**: triangular superior con los autovalores en la diagonal, ordenados por módulo. Por qué funciona: $P_kU_k$ es una factorización QR de $A^{k+1}$, así que **QR hace el método de las potencias sobre todos los subespacios $\langle e_1,\dots,e_j\rangle$ a la vez** (esquina $(1,1)$: potencias; esquina $(n,n)$: potencias inversas).

## El algoritmo

$A_0=A$; para $k\geq0$: $A_k=Q_kR_k$ (QR), $A_{k+1}=R_kQ_k$. Como $R_k=Q_k^HA_k$,
$$
A_{k+1}=Q_k^HA_kQ_k.
$$

> [!algoritmo] QR básico ([[7.7 - Implementar el algoritmo QR|Problema 7.7]])
> ```julia
> function qr_basico(A; tol=1e-10, kmax=100_000)
>     Ak = copy(A); d = diag(Ak); k = 0
>     while k < kmax
>         k += 1
>         F = qr(Ak)
>         Ak = F.R * Matrix(F.Q)          # A_{k+1} = R_k Q_k
>         dn = diag(Ak)
>         maximum(abs.(dn - d)) < tol && break
>         d = dn
>     end
>     return Ak, k
> end
> ```

## Teorema 7.10

Con $P_k:=Q_0Q_1\cdots Q_k$ (unitaria) y $U_k:=R_kR_{k-1}\cdots R_0$ (triangular superior):
- (a) $A_{k+1}=Q_k^HA_kQ_k$;
- (b) $A_{k+1}=P_k^HAP_k$;
- (c) $A^{k+1}=P_kU_k$. **Atención: es una potencia de $A$.**

(c) dice que **$P_kU_k$ es una factorización QR de $A^{k+1}$**: el algoritmo QR calcula las QR de $A,A^2,A^3,\dots$ sin formar nunca esas potencias. Demostración: [[Teorema 7.10]].

## Teorema 7.11 (convergencia)

> Hipótesis: (a) $|\lambda_1|>|\lambda_2|>\cdots>|\lambda_n|$ (módulos **distintos**); (b) $A=XDX^{-1}$ con $D=\text{diag}(\lambda_i)$ y **$X^{-1}$ admite factorización LU**. Entonces existen matrices de fase $S_k=\text{diag}(\sigma_i^{(k)})$, $|\sigma_i^{(k)}|=1$, con $S_{k-1}^HQ_kS_k\to I$ y
> $$S_k^HR_kS_{k-1}\to T,\qquad S_{k-1}^HA_kS_{k-1}\to T=\begin{pmatrix}\lambda_1&\times&\cdots&\times\\0&\lambda_2&\ddots&\vdots\\\vdots&\ddots&\ddots&\times\\0&\cdots&0&\lambda_n\end{pmatrix}.$$

**No se demuestra en el curso** (referencia: [SB, Thm. 6.6.4.12]). Ver [[Teorema 7.11]].

## Observación 7.12: la conexión con las potencias

- **Esquina $(1,1)$ = método de las potencias.** De $A^{k+1}e_1=P_kU_ke_1=(U_k)_{11}P_ke_1$: la primera columna de $P_k$ es $A^{k+1}e_1$ normalizada, o sea el iterado de las potencias con $u_0=e_1$. Y $(A_{k+1})_{11}=(P_ke_1)^HA(P_ke_1)$ es su **cociente de Rayleigh**. Luego $(A_k)_{11}\to\lambda_1$ con tasa $|\lambda_2/\lambda_1|$. La hipótesis (b) es exactamente $y_1^He_1=(LU)_{11}=u_{11}\neq0$.
- **Esquina $(n,n)$ = potencias inversas.** Invirtiendo y conjugando: $(A^{-H})^{k+1}e_n=(U_k^{-H})_{nn}P_ke_n$, así que la última columna de $P_k$ es el método de las potencias para $B=A^{-H}$, cuyo dominante es $\overline{\lambda_n}^{-1}$. Luego $(A_k)_{nn}\to\lambda_n$ con tasa $|\lambda_n/\lambda_{n-1}|$. Por eso **los shifts se toman de la esquina $(n,n)$**.
- **El resto:** las primeras $j$ columnas de $P_k$ generan $A^{k+1}\langle e_1,\dots,e_j\rangle$: QR itera simultáneamente sobre $\langle e_1\rangle\subset\langle e_1,e_2\rangle\subset\cdots$, ortonormalizando. Cada subespacio converge al invariante dominante de su dimensión con tasa $|\lambda_{j+1}/\lambda_j|$.

Detalle: [[Observacion 7.12]].

## Observación 7.13: velocidad y módulos iguales

$$
\left\|S_{k-1}^HA_kS_{k-1}-T\right\|\leq C\Big(\max_{1\leq i\leq n-1}\frac{|\lambda_{i+1}|}{|\lambda_i|}\Big)^k,\qquad(A_k)_{ii}\to\lambda_i.
$$
La hipótesis (b) solo sirve para que la diagonal salga **ordenada** por módulo.

**Si $|\lambda_r|=|\lambda_{r+1}|$** (y los demás distintos): $a^{(k)}_{ij}\to0$ para $i>j$ salvo $a^{(k)}_{r+1,r}$; $a_{ii}^{(k)}\to\lambda_i$ para $i\neq r,r+1$; y el bloque $2\times2$ $\begin{pmatrix}a_{rr}&a_{r,r+1}\\a_{r+1,r}&a_{r+1,r+1}\end{pmatrix}$ en general **no converge**, pero **sus autovalores sí** convergen a $\lambda_r,\lambda_{r+1}$. (Es lo que pasa con pares complejos conjugados de una matriz real.)

> [!idea] En hermíticas la diagonal converge el doble de rápido
> Para $A$ hermítica, $A_k$ es hermítica ([[7.15 - QR preserva hermiticidad|Problema 7.15]]) y converge a $\text{diag}(\lambda_i)$; los elementos diagonales son cocientes de Rayleigh y su error es el **cuadrado** del de las columnas de $P_k$. En $A_N$ ([[7.8 - QR en la matriz A_N|Problema 7.8]]), $(A_k)_{11}$ converge con tasa $(\lambda_2/\lambda_1)^2$, verificado a 5 decimales para $N=5,15,50$.

> [!pendiente] Parar por "la diagonal no cambia" deja la parte de abajo grande
> En una tridiagonal simétrica, el QR básico con criterio "$\max|\Delta\,\text{diag}|<10^{-12}$" se detuvo con elementos fuera de la diagonal de tamaño $\sim10^{-6}$: por lo anterior, la diagonal se mueve como el **cuadrado** de la subdiagonal. Los autovalores están bien a $10^{-12}$, pero $A_k$ todavía no es diagonal a esa precisión ([[7.7 - Implementar el algoritmo QR|Problema 7.7]]).

## Ver también
- [[Algoritmo QR - Hessenberg, shifts y deflacion]] (cómo se hace eficiente)
- [[Metodo de las potencias]], [[Metodo de las potencias inversas]]
- [[Teorema 4.21]] (Schur: lo que QR construye)
- [[Capitulo 7 - MOC]]
