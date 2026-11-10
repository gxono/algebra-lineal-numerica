---
tags: [capitulo5, teorema]
---

# Lema 5.27

Si $r_k\neq0$, entonces
$$
\alpha_{k+1}=\frac{\|r_k\|_2^2}{\|d_{k+1}\|_A^2}=\frac{r_k^Tr_k}{d_{k+1}^TAd_{k+1}}, \qquad \beta_{k+1}=\frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2}=\frac{r_k^Tr_k}{r_{k-1}^Tr_{k-1}}.
$$

> [!demostracion]-
> **Fórmula de $\alpha_{k+1}$.** Por el [[Teorema 5.26]], $d_{k+1}=r_k+\beta_{k+1}d_k$ con $d_k\in\mathcal K_k$. Entonces
> $$d_{k+1}^Tr_k = (r_k+\beta_{k+1}d_k)^Tr_k = r_k^Tr_k+\beta_{k+1}\underbrace{d_k^Tr_k}_{=0} = r_k^Tr_k,$$
> donde $d_k^Tr_k=0$ porque $d_k\in\mathcal K_k$ (por [[Lema 5.25]]) y $r_k\perp\mathcal K_k$ (por [[Lema 5.10]]). Reemplazando en la fórmula del [[Corolario 5.24]], $\alpha_{k+1}=\dfrac{d_{k+1}^Tr_k}{d_{k+1}^TAd_{k+1}}=\dfrac{r_k^Tr_k}{d_{k+1}^TAd_{k+1}}$.
>
> **Fórmula de $\beta_{k+1}$.** Como $r_k=b-Ax_k=b-Ax_{k-1}+A(x_{k-1}-x_k)=r_{k-1}-\alpha_kAd_k$ (usando $x_k=x_{k-1}+\alpha_kd_k$ del [[Corolario 5.24]]), se despeja
> $$Ad_k=\frac{1}{\alpha_k}(r_{k-1}-r_k).$$
> Usando esto y la fórmula ya probada $\alpha_k=\dfrac{r_{k-1}^Tr_{k-1}}{d_k^TAd_k}$ (aplicada un paso antes) junto con $r_k\perp r_{k-1}$ (caso particular de [[Lema 5.10]]):
> $$\beta_{k+1} = -\frac{r_k^TAd_k}{d_k^TAd_k} = -\frac{r_k^T\frac{1}{\alpha_k}(r_{k-1}-r_k)}{d_k^TAd_k} = -\frac{\frac{1}{\alpha_k}\big(r_k^Tr_{k-1}-r_k^Tr_k\big)}{d_k^TAd_k} = \frac{r_k^Tr_k}{\alpha_k\,d_k^TAd_k} = \frac{r_k^Tr_k}{r_{k-1}^Tr_{k-1}},$$
> donde se usó $r_k^Tr_{k-1}=0$ y, en el último paso, $\alpha_k\,d_k^TAd_k=r_{k-1}^Tr_{k-1}$ (despejando de $\alpha_k=r_{k-1}^Tr_{k-1}/d_k^TAd_k$).

> [!idea]
> Esta es la versión **implementable**: ni $\alpha_{k+1}$ ni $\beta_{k+1}$ necesitan guardar $d$'s viejos ni calcular $Ad_j$ para $j<k$ - solo $\|r_k\|_2^2$ en el paso actual y el anterior (y, para $\alpha$, el producto $d_{k+1}^TAd_{k+1}$ que de todas formas hay que calcular para avanzar). Esto es exactamente lo que se implementa en [[Implementacion del gradiente conjugado]].

## Ver también
- [[Implementacion del gradiente conjugado]] (nota de tema completa)
- [[Teorema 5.26]]
- [[Corolario 5.24]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
