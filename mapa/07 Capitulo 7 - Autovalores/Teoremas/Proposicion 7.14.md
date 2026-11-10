---
tags: [capitulo7, teorema]
---

# Proposición 7.14 (QR con shifts)

Con $A_k-\mu_kI=Q_kR_k$, $A_{k+1}=R_kQ_k+\mu_kI$, $P_k=Q_0\cdots Q_k$, $U_k=R_k\cdots R_0$:
$$
(A-\mu_kI)(A-\mu_{k-1}I)\cdots(A-\mu_0I)=P_kU_k.
$$

> [!demostracion]-
> **Semejanza (previo).** De $R_k=Q_k^H(A_k-\mu_kI)$: $A_{k+1}=Q_k^H(A_k-\mu_kI)Q_k+\mu_kI=Q_k^HA_kQ_k$. Luego, igual que en el [[Teorema 7.10]](b), $A_k=P_{k-1}^HAP_{k-1}$.
>
> **Inducción.** $k=0$: $A-\mu_0I=A_0-\mu_0I=Q_0R_0=P_0U_0$. Sea $M_k:=(A-\mu_kI)\cdots(A-\mu_0I)$ y supongamos $M_{k-1}=P_{k-1}U_{k-1}$. De $A_k=P_{k-1}^HAP_{k-1}$ sale $(A-\mu_kI)P_{k-1}=P_{k-1}(A_k-\mu_kI)$, y
> $$M_k=(A-\mu_kI)M_{k-1}=(A-\mu_kI)P_{k-1}U_{k-1}=P_{k-1}(A_k-\mu_kI)U_{k-1}=P_{k-1}Q_kR_kU_{k-1}=P_kU_k.$$

**Consecuencia.** Repitiendo la cuenta de la [[Observacion 7.12]] con $M_k$ en lugar de $A^{k+1}$:
$$
M_k^{-H}e_n=(U_k^{-H})_{nn}P_ke_n,\qquad M_k^{-H}=(A^H-\bar\mu_0I)^{-1}\cdots(A^H-\bar\mu_kI)^{-1}.
$$
La última columna de $P_k$ es **potencias inversas con shifts** sobre $A^H$; en el paso $k$ el error se contrae en $\approx\max_{i\neq n}\left|\frac{\lambda_n-\mu_k}{\lambda_i-\mu_k}\right|$. Con $\mu_k=(A_k)_{nn}$ (cociente de Rayleigh de la última columna de $P_{k-1}$) es la iteración del cociente de Rayleigh, y el paso de QR con shift cuesta lo mismo que sin shift: no hay que refactorizar.

## Ver también
- [[Algoritmo QR - Hessenberg, shifts y deflacion]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
