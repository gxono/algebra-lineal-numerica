---
tags: [capitulo7, teorema]
---

# Teorema 7.11 (convergencia del algoritmo QR)

Sea $A\in\mathbb{C}^{n\times n}$ tal que:
- (a) $|\lambda_1|>|\lambda_2|>\cdots>|\lambda_n|$;
- (b) $A=XDX^{-1}$ con $D=\text{diag}(\lambda_1,\dots,\lambda_n)$, y $X^{-1}$ tiene factorización LU ($X^{-1}=LU$, $L$ triangular inferior unitaria, $U$ triangular superior).

Entonces existen matrices de fase $S_k=\text{diag}(\sigma_1^{(k)},\dots,\sigma_n^{(k)})$, $|\sigma_i^{(k)}|=1$, tales que $\lim_kS_{k-1}^HQ_kS_k=I$ y
$$
\lim_{k\to\infty}S_k^HR_kS_{k-1}=\lim_{k\to\infty}S_{k-1}^HA_kS_{k-1}=\begin{pmatrix}\lambda_1&\times&\cdots&\times\\0&\lambda_2&\ddots&\vdots\\\vdots&\ddots&\ddots&\times\\0&\cdots&0&\lambda_n\end{pmatrix}.
$$

> [!pendiente] No se demuestra en el curso
> El apunte remite a [SB, Theorem 6.6.4.12] (Stoer–Bulirsch). Lo que **sí** hay que saber: la cota de velocidad $C\big(\max_i|\lambda_{i+1}/\lambda_i|\big)^k$, qué pasa si dos autovalores tienen igual módulo (bloque $2\times2$ que no converge pero cuyos autovalores sí), para qué sirve (b) (diagonal ordenada por módulo), y la **explicación** de por qué funciona vía potencias / potencias inversas ([[Observacion 7.12]]), que sí está completa.

Qué pasa cuando fallan las hipótesis, experimentalmente:
- (a) falla por un autovalor doble **semisimple** ([[7.10 - QR con autovalor semisimple doble|Problema 7.10]]): converge igual, rápido.
- (a) falla por un autovalor **defectuoso** ([[7.11 - QR con autovalor defectuoso|Problema 7.11]]): converge, pero **sublinealmente** ($\sim1/k$).
- (a) falla por $\pm$ pares de igual módulo ($A_N-2I$, [[7.9 - QR con shift y deflacion|Problema 7.9]]): no converge.

## Ver también
- [[Algoritmo QR - definicion y convergencia]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
