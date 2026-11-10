---
tags: [capitulo7, teorema]
---

# Teorema 7.10 (QR calcula las QR de las potencias de $A$)

Sean $A_0=A$, $A_k=Q_kR_k$, $A_{k+1}=R_kQ_k$, y $P_k=Q_0Q_1\cdots Q_k$, $U_k=R_kR_{k-1}\cdots R_0$. Entonces:
- (a) $A_{k+1}=Q_k^HA_kQ_k$ (unitariamente semejante a $A_k$);
- (b) $A_{k+1}=P_k^HAP_k$;
- (c) $A^{k+1}=P_kU_k$.

> [!demostracion]-
> **(a)** $Q_k$ unitaria y $A_k=Q_kR_k$ ⟹ $R_k=Q_k^HA_k$ ⟹ $A_{k+1}=R_kQ_k=Q_k^HA_kQ_k$.
>
> **(b)** Inducción. $k=0$: es (a), con $P_0=Q_0$, $A_0=A$. Si $A_k=P_{k-1}^HAP_{k-1}$, por (a):
> $$A_{k+1}=Q_k^HA_kQ_k=Q_k^HP_{k-1}^HAP_{k-1}Q_k=(P_{k-1}Q_k)^HA(P_{k-1}Q_k)=P_k^HAP_k.$$
>
> **(c)** Inducción. $k=0$: $A^1=A_0=Q_0R_0=P_0U_0$. Supongamos $A^k=P_{k-1}U_{k-1}$. Multiplicando $A_k=P_{k-1}^HAP_{k-1}$ por $P_{k-1}$ a izquierda: $AP_{k-1}=P_{k-1}A_k$. Entonces
> $$A^{k+1}=A\,A^k=AP_{k-1}U_{k-1}=P_{k-1}A_kU_{k-1}=P_{k-1}Q_kR_kU_{k-1}=(P_{k-1}Q_k)(R_kU_{k-1})=P_kU_k.$$

$P_k$ es unitaria (producto de unitarias) y $U_k$ triangular superior (producto de triangulares superiores): **(c) es una factorización QR de $A^{k+1}$**. Es la raíz de la conexión con los métodos de las potencias ([[Observacion 7.12]]).

## Ver también
- [[Algoritmo QR - definicion y convergencia]]
- [[Proposicion 7.14]] (la versión con shifts)
- [[Demostraciones - Capitulo 7|Demostraciones]]
