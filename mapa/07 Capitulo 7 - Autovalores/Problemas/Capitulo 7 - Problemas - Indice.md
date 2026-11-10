---
tags: [capitulo7, moc]
---

# Capítulo 7 - Índice de problemas

| # | Tema |
|---|---|
| [[7.1 - Programar el metodo de las potencias\|7.1]] | Potencias en $5\times5$: tasa $0.8$ (autovector), $0.64$ (autovalor, simétrica) |
| [[7.2 - Potencias en la matriz A_N\|7.2]] | Autovectores de $A_N$; potencias: tasa $\to1$, 5756 iteraciones con $N=50$ |
| [[7.3 - Programar potencias inversas\|7.3]] | Tasa $\rho_\mu$ verificada para tres shifts |
| [[7.4 - Potencias inversas en la matriz A_N\|7.4]] | Tasa $\to\frac14$ independiente de $N$ (≈17 iteraciones) |
| [[7.5 - Potencias inversas con shift actualizado\|7.5]] | Iteración del cociente de Rayleigh: orden 3 (simétrica), 2 (no simétrica) |
| [[7.6 - Potencias inversas con y sin shift en A_N\|7.6]] | Shift actualizado: 4–7 iteraciones pero converge a **otro** autovalor |
| [[7.7 - Implementar el algoritmo QR\|7.7]] | QR conserva Hessenberg; tridiagonal simétrica sigue tridiagonal simétrica |
| [[7.8 - QR en la matriz A_N\|7.8]] | **Demostración** de los autopares de $A_N$; tasas de las esquinas |
| [[7.9 - QR con shift y deflacion\|7.9]] | Rayleigh **se estanca** en $A_N$ (diagonal clavada en 2); Wilkinson: 8/33/102 iteraciones; cúbica |
| [[7.10 - QR con autovalor semisimple doble\|7.10]] | Autovalor doble semisimple: QR converge igual |
| [[7.11 - QR con autovalor defectuoso\|7.11]] | Autovalor defectuoso: convergencia sublineal $\sim1/k$ |
| [[7.12 - Unicidad de QR salvo fases\|7.12]] | $\tilde Q=QS$, $\tilde R=S^HR$ |
| [[7.13 - Producto de matrices Hessenberg\|7.13]] | **Enunciado falso**: contraejemplo; lo cierto es Hessenberg × triangular |
| [[7.14 - QR preserva Hessenberg\|7.14]] | Sí con Householder / $A$ no singular; **no** en general si $A$ es singular (contraejemplo) |
| [[7.15 - QR preserva hermiticidad\|7.15]] | $A_k\to\text{diag}(\lambda_i)$, $P_k\to$ autovectores ortonormales |

Volver a [[Capitulo 7 - MOC]].
