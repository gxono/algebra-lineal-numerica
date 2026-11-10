---
tags: [glosario, capitulo7]
aliases: [desplazamiento, traslación, shift de Wilkinson]
---

# Shift (desplazamiento)

## Definición

Reemplazar $A$ por $A-\mu I$: mismos autovectores, autovalores trasladados $\lambda_i-\mu$. Combinado con inversión ($(A-\mu I)^{-1}$) vuelve dominante al autovalor más cercano a $\mu$ (potencias inversas). En QR: $A_k-\mu_kI=Q_kR_k$, $A_{k+1}=R_kQ_k+\mu_kI$. **Shift de Rayleigh:** $\mu_k=a_{nn}^{(k)}$. **Shift de Wilkinson:** autovalor del bloque $2\times2$ inferior derecho más cercano a $a_{nn}^{(k)}$, que no se estanca con espectros simétricos.

## Ver también
- [[Metodo de las potencias inversas]]
- [[Algoritmo QR - Hessenberg, shifts y deflacion]]
- [[Glosario - Indice|Glosario]]
