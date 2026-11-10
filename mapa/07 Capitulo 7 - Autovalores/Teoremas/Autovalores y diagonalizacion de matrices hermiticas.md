---
tags: [capitulo7, teorema]
---

# Matrices hermíticas: autovalores reales y diagonalización unitaria

Sea $A\in\mathbb{C}^{n\times n}$ con $A^H=A$. Entonces:
1. $\sigma(A^H)=\overline{\sigma(A)}$ para cualquier $A$; en particular el espectro de una hermítica es invariante por conjugación.
2. Los autovalores de $A$ son **reales**.
3. Existe $U$ unitaria con $U^HAU=\text{diag}(\lambda_1,\dots,\lambda_n)$, $\lambda_i\in\mathbb{R}$: **base ortonormal de autovectores**.

> [!demostracion]-
> **1.** $\det(\bar\lambda I-A^H)=\det\big((\lambda I-A)^H\big)=\overline{\det(\lambda I-A)}$, así que $\lambda\in\sigma(A)\iff\bar\lambda\in\sigma(A^H)$. (Que el espectro sea invariante por conjugación **no** alcanza para concluir que es real: el de cualquier matriz real también lo es.)
>
> **2.** Si $Ax=\lambda x$, $x\neq0$:
> $$x^HAx=\lambda x^Hx=\lambda\|x\|_2^2,\qquad x^HAx=(A^Hx)^Hx=(Ax)^Hx=(\lambda x)^Hx=\bar\lambda\|x\|_2^2.$$
> Como $\|x\|_2\neq0$, $\lambda=\bar\lambda$.
>
> **3.** Por el [[Teorema 4.21]] (Schur) existe $U$ unitaria con $U^HAU=T$ triangular superior y $\text{diag}(T)=(\lambda_1,\dots,\lambda_n)$. Pero $T^H=U^HA^HU=U^HAU=T$: $T$ es hermítica. Triangular superior y hermítica (su traspuesta conjugada, triangular inferior, es ella misma) ⟹ **diagonal**, con diagonal real. Las columnas de $U$ son autovectores ortonormales.

**Autovectores a izquierda.** Si $A^Hy=\bar\lambda y$, conjugando y trasponiendo, $y^HA=\lambda y^H$: los autovectores a izquierda de $A$ son los conjugados traspuestos de los autovectores de $A^H$. Para $A$ hermítica, a izquierda y a derecha coinciden.

## Ver también
- [[Autovalores - repaso y por que los metodos son iterativos]]
- [[7.15 - QR preserva hermiticidad]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
