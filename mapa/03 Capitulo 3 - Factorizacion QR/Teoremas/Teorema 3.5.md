---
tags: [capitulo3, teorema]
---

# Teorema 3.5 (existencia de QR vía Gram-Schmidt)

Si las columnas de $A\in\mathbb{R}^{n\times m}$ son li, existe $Q_1\in\mathbb{R}^{n\times m}$ ortogonal y $R_1\in\mathbb{R}^{m\times m}$ triangular superior con $A=Q_1R_1$. *(Completando $Q_1$ a una base ortonormal de $\mathbb{R}^n$ se obtiene la versión cuadrada $A=QR$ de la Observación 3.6.)*

> [!demostracion]- Construcción explícita (el procedimiento de Gram-Schmidt), por inducción en $k$
> Mostramos que el procedimiento de arriba está bien definido y produce $u_1,\dots,u_m$ ortonormales con $v_k\in\text{span}\{u_1,\dots,u_k\}$ para cada $k$ - de ahí sale todo el teorema.
>
> **Hipótesis inductiva:** $u_1,\dots,u_{k-1}$ ya construidos, ortonormales, con $\text{span}\{u_1,\dots,u_{k-1}\}=\text{span}\{v_1,\dots,v_{k-1}\}$.
>
> **$w_k\neq0$ (buena definición de $u_k$).** Si $w_k=0$, entonces $v_k=\sum_{i<k}(v_k\cdot u_i)u_i\in\text{span}\{u_1,\dots,u_{k-1}\}\overset{\text{H.I.}}{=}\text{span}\{v_1,\dots,v_{k-1}\}$, contradiciendo que $v_1,\dots,v_k$ son li. Luego $w_k\neq0$ y $u_k:=w_k/\|w_k\|$ está bien definido, con $\|u_k\|=1$.
>
> **Ortogonalidad.** Para $i<k$, usando que $u_1,\dots,u_{k-1}$ son ortonormales (H.I.):
> $$u_k\cdot u_i = \frac1{\|w_k\|}\Big(v_k\cdot u_i-\sum_{j<k}(v_k\cdot u_j)(u_j\cdot u_i)\Big) = \frac1{\|w_k\|}\big(v_k\cdot u_i-(v_k\cdot u_i)\cdot1\big)=0,$$
> (el término $j=i$ sobrevive con coeficiente $1$, los demás $u_j\cdot u_i=0$). Así, $u_1,\dots,u_k$ son ortonormales.
>
> **Triangularidad.** Despejando, $v_k=w_k+\sum_{i<k}(v_k\cdot u_i)u_i=\|w_k\|\,u_k+\sum_{i<k}(v_k\cdot u_i)u_i\in\text{span}\{u_1,\dots,u_k\}$, con coeficiente de $u_k$ igual a $\|w_k\|>0$.
>
> **Conclusión.** Repitiendo para $k=1,\dots,m$: cada columna $v_k$ de $A$ es combinación de $u_1,\dots,u_k$ únicamente, con coeficiente de $u_k$ no nulo - esto es exactamente $A=Q_1R_1$ con $Q_1=[u_1\,\cdots\,u_m]$ (columnas ortonormales, i.e. $Q_1^TQ_1=I$) y $R_1$ triangular superior con diagonal $\|w_1\|,\dots,\|w_m\|>0$ (en particular $R_1$ invertible).

## Ver también
- [[Ortogonalizacion de Gram-Schmidt]] (nota de tema completa)
- [[Demostraciones - Capitulo 3|Demostraciones]]
