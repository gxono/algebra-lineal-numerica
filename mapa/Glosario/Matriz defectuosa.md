---
tags: [glosario, capitulo7]
aliases: [defectuosa, multiplicidad geométrica, multiplicidad algebraica]
---

# Matriz defectuosa (multiplicidades)

## Definición

**Multiplicidad algebraica** de $\lambda$: multiplicidad como raíz de $\det(\lambda I-A)$. **Geométrica:** $\dim\ker(A-\lambda I)=n-\text{rank}(A-\lambda I)$. Siempre geométrica $\leq$ algebraica; si es $<$, $\lambda$ (y $A$) es **defectuoso**, y $A$ no es diagonalizable (tiene un bloque de Jordan de tamaño $\geq2$). Ejemplo: $\begin{pmatrix}2&0\\1&2\end{pmatrix}$. El algoritmo QR converge solo como $1/k$ en un autovalor defectuoso.

## Ver también
- [[Autovalores - repaso y por que los metodos son iterativos]]
- [[7.11 - QR con autovalor defectuoso]]
- [[Glosario - Indice|Glosario]]
