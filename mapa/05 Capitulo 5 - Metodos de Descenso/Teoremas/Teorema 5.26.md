---
tags: [capitulo5, teorema]
---

# Teorema 5.26

Sea $k\in\mathbb N$. Si $r_{k-1}\neq0$, entonces
$$
d_{k+1}=r_k+\beta_{k+1}d_k, \qquad \beta_{k+1}=-\frac{r_k^TAd_k}{d_k^TAd_k}.
$$

> [!demostracion]-
> Sabemos (de la definición de $d_{k+1}$, [[Corolario 5.24]]) que existe un único $w_k\in\mathcal K_k$ tal que $d_{k+1}=r_k+w_k$. Por el [[Lema 5.25]], $\mathcal K_k=\text{span}\{d_1,\dots,d_k\}$, así que $w_k=\sum_{i=1}^k\gamma_id_i$ para ciertos $\gamma_1,\dots,\gamma_k\in\mathbb R$, y
> $$d_{k+1}=r_k+\sum_{i=1}^k\gamma_id_i.$$
> Por definición, $d_{k+1}$ es $A$-ortogonal a $\mathcal K_k=\text{span}\{d_1,\dots,d_k\}$, así que $d_{k+1}^TAd_j=0$ para $j=1,\dots,k$. Entonces, para cada $j\in\{1,\dots,k\}$:
> $$0 = d_{k+1}^TAd_j = \Big(r_k+\sum_{i=1}^k\gamma_id_i\Big)^TAd_j = r_k^TAd_j+\sum_{i=1}^k\gamma_i\, d_i^TAd_j.$$
> Como $\{d_1,\dots,d_k\}$ es $A$-ortogonal (por el [[Lema 5.25]]), $d_i^TAd_j=0$ salvo cuando $i=j$, así que la suma colapsa al término $i=j$:
> $$0 = r_k^TAd_j+\gamma_jd_j^TAd_j \quad\Longrightarrow\quad \gamma_j=-\frac{r_k^TAd_j}{d_j^TAd_j}, \qquad j=1,\dots,k.$$
> **"El milagro":** falta ver que $\gamma_j=0$ para $j=1,\dots,k-1$ (solo sobrevive $\gamma_k$). Equivalentemente, basta ver que $r_k^TAd_j=0$ para $j=1,\dots,k-1$. Esto vale porque $Ad_j\in\mathcal K_{j+1}\subset\mathcal K_k$ para $j=1,\dots,k-1$ (ya que $d_j\in\mathcal K_j$ por el [[Lema 5.25]], y $A\mathcal K_j\subset\mathcal K_{j+1}$ por la definición del espacio de Krylov), mientras que $r_k\perp\mathcal K_k$ (por [[Lema 5.10]], $r_k\perp\{r_0,\dots,r_{k-1}\}=\mathcal K_k$). Así, $r_k^TAd_j=0$ para $j<k$, y queda
> $$d_{k+1}=r_k+\gamma_kd_k, \qquad \gamma_k=-\frac{r_k^TAd_k}{d_k^TAd_k}=:\beta_{k+1}.$$

> [!idea]
> Esto es lo que hace viable al gradiente conjugado: **sin** este resultado, $d_{k+1}$ necesitaría los $k$ coeficientes $\gamma_1,\dots,\gamma_k$ (y por lo tanto, guardar $d_1,\dots,d_k$ enteros - como en GMRes). Con él, alcanza con el $d_k$ más reciente: **recurrencia de dos términos**, memoria $O(n)$ en vez de $O(nk)$.

## Ver también
- [[Implementacion del gradiente conjugado]] (nota de tema completa)
- [[Lema 5.25]]
- [[Lema 5.27]] (fórmula de $\beta_{k+1}$ solo en términos de residuos)
- [[Demostraciones - Capitulo 5|Demostraciones]]
