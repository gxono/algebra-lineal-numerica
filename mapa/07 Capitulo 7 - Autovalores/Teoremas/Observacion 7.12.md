---
tags: [capitulo7, observacion]
---

# Observación 7.12 (QR = potencias en la esquina $(1,1)$ y potencias inversas en la $(n,n)$)

Bajo las hipótesis del [[Teorema 7.11]] (y $A$ invertible para la segunda parte):
1. $(A_{k+1})_{11}$ es el cociente de Rayleigh del iterado del método de las potencias con $u_0=e_1$; luego $(A_k)_{11}\to\lambda_1$ con tasa $|\lambda_2/\lambda_1|$.
2. $(A_{k+1})_{nn}$ es el cociente de Rayleigh del iterado del método de las potencias para $B=A^{-H}$ con $u_0=e_n$; luego $(A_k)_{nn}\to\lambda_n$ con tasa $|\lambda_n/\lambda_{n-1}|$.

> [!demostracion]-
> **1. Esquina $(1,1)$.** Por el [[Teorema 7.10]](c), $A^{k+1}=P_kU_k$. Como $U_k$ es triangular superior, $U_ke_1=(U_k)_{11}e_1$, así que
> $$A^{k+1}e_1=P_kU_ke_1=(U_k)_{11}P_ke_1.$$
> $P_ke_1$ tiene norma 1 (columna de una unitaria), luego es $A^{k+1}e_1$ normalizado: por el [[Lema 7.1]], salvo un factor de módulo uno, el iterado de las potencias con $u_0=e_1$. Por el Teorema 7.10(b), $(A_{k+1})_{11}=(P_ke_1)^HA(P_ke_1)$: su cociente de Rayleigh (que no cambia al multiplicar el vector por un escalar de módulo uno, [[Observacion 7.4]]). Por la [[Proposicion 7.2]], $(A_{k+1})_{11}\to\lambda_1$ con tasa $|\lambda_2/\lambda_1|$. La hipótesis $y_1^Hu_0\neq0$ es $y_1^He_1=(X^{-1})_{11}=(LU)_{11}=u_{11}\neq0$, que vale porque $X^{-1}$ tiene LU (con $U$ no singular): **esto explica la hipótesis (b)**.
>
> **2. Esquina $(n,n)$.** Invirtiendo $A^{k+1}=P_kU_k$ ($P_k$ unitaria): $(A^{k+1})^{-1}=U_k^{-1}P_k^H$; conjugando y trasponiendo:
> $$(A^{-H})^{k+1}=P_k(U_k^{-1})^H=P_kU_k^{-H}.$$
> $U_k^{-1}$ es triangular superior, así que $U_k^{-H}$ es triangular **inferior** y su última columna es múltiplo de $e_n$: $U_k^{-H}e_n=(U_k^{-H})_{nn}e_n$, con $(U_k^{-H})_{nn}=\overline{1/(U_k)_{nn}}\neq0$. Entonces
> $$(A^{-H})^{k+1}e_n=(U_k^{-H})_{nn}P_ke_n:$$
> la última columna de $P_k$ es el iterado de las potencias para $B=A^{-H}$ con $u_0=e_n$.
>
> *Autopares de $B$.* Si $y_i^HA=\lambda_iy_i^H$, entonces $A^Hy_i=\bar\lambda_iy_i$ y $By_i=\bar\lambda_i^{-1}y_i$. Como $|\bar\lambda_i^{-1}|=1/|\lambda_i|$, la hipótesis (a) da $\frac{1}{|\lambda_n|}>\cdots>\frac{1}{|\lambda_1|}$: el dominante de $B$ es $\bar\lambda_n^{-1}$, con autovector $y_n$.
>
> *No degeneración.* El autovector a izquierda de $B$ para $\bar\lambda_n^{-1}$: $z^HB=\bar\lambda_n^{-1}z^H\iff A^{-1}z=\lambda_n^{-1}z\iff Az=\lambda_nz$, o sea $z=x_n$. Hace falta $x_n^He_n\neq0$, es decir $X_{nn}\neq0$. De $X^{-1}=LU$: $X=U^{-1}L^{-1}$ y $X_{nn}=(U^{-1})_{nn}(L^{-1})_{nn}=1/u_{nn}\neq0$. **La misma hipótesis (b)** garantiza la no degeneración en los dos extremos.
>
> *Conclusión.* Por la Proposición 7.2 aplicada a $B$, $P_ke_n\to y_n$ en dirección con tasa $\left|\bar\lambda_{n-1}^{-1}/\bar\lambda_n^{-1}\right|=|\lambda_n|/|\lambda_{n-1}|$. Con $y_n$ unitario, $y_n^HAy_n=\lambda_ny_n^Hy_n=\lambda_n$, y por la Observación 7.4:
> $$(A_{k+1})_{nn}=(P_ke_n)^HA(P_ke_n)\to y_n^HAy_n=\lambda_n.$$

**El resto de la matriz.** Como $U_k$ es triangular superior, las primeras $j$ columnas de $P_k$ generan $A^{k+1}\langle e_1,\dots,e_j\rangle$: QR es el método de las potencias **sobre subespacios encajados**, ortonormalizando en cada paso; cada uno converge al subespacio invariante dominante de su dimensión con tasa $|\lambda_{j+1}/\lambda_j|$. Que todos converjan es que $A_k$ se vuelva triangular.

## Ver también
- [[Algoritmo QR - definicion y convergencia]]
- [[Algoritmo QR - Hessenberg, shifts y deflacion]] (por eso los shifts salen de la esquina $(n,n)$)
- [[Demostraciones - Capitulo 7|Demostraciones]]
