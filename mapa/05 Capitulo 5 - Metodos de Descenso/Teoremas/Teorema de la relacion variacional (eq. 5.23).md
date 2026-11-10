---
tags: [capitulo5, teorema]
---

# Teorema (eq. 5.23): relación variacional entre restricción y prolongación

$$
I_h^{2h} = \frac12\big(I_{2h}^h\big)^T.
$$
Garantiza que el operador de dos mallas $I_h^{2h}A_hI_{2h}^h$ sea simétrico - esencial para que el ciclo en V sea sdp y se pueda usar como precondicionador del GC.

> [!demostracion]- Demostración (a partir de las fórmulas explícitas)
> Sea $P:=I_{2h}^h\in\mathbb{R}^{n\times n_c}$. Por la definición de prolongación, la columna $k$ de $P$ (es decir, $Pe_k$, la imagen del $k$-ésimo vector canónico de la malla gruesa) tiene a lo sumo tres entradas no nulas, en las filas $2k-1$, $2k$, $2k+1$ de la malla fina:
> $$P_{2k,k}=1 \quad(\text{fila par, copia directa}), \qquad P_{2k-1,k}=P_{2k+1,k}=\tfrac12 \quad(\text{filas impares vecinas, promedio}).$$
> (Para $k=1$ la fila $2k-1=1$ usa $v_0=0$ en la fórmula de prolongación, pero el coeficiente que multiplica a $v_1$ en esa fila sigue siendo $\tfrac12$; análogamente para $k=n_c$ en la fila $2k+1=n$. Los casos de borde no cambian los coeficientes, solo anulan al vecino inexistente.)
>
> Transponiendo, la fila $k$ de $P^T$ tiene esas mismas tres entradas en las columnas $2k-1,2k,2k+1$:
> $$(P^T)_{k,2k-1}=\tfrac12, \qquad (P^T)_{k,2k}=1, \qquad (P^T)_{k,2k+1}=\tfrac12.$$
> Multiplicando por $\tfrac12$:
> $$\Big(\tfrac12P^T\Big)_{k,2k-1}=\tfrac14, \qquad \Big(\tfrac12P^T\Big)_{k,2k}=\tfrac12, \qquad \Big(\tfrac12P^T\Big)_{k,2k+1}=\tfrac14.$$
> Esto es exactamente la fila $j=k$ de $I_h^{2h}$: $(I_h^{2h})_{k,2k-1}=\tfrac14$, $(I_h^{2h})_{k,2k}=\tfrac12\ (=\tfrac24)$, $(I_h^{2h})_{k,2k+1}=\tfrac14$, y cero en el resto. Como ambas matrices tienen las mismas entradas en todas las posiciones, $I_h^{2h}=\tfrac12P^T=\tfrac12(I_{2h}^h)^T$.

Verificación numérica (reconstruyendo las matrices y comparando) en [[5.12 - Multigrid como metodo iterativo|Problema 5.12]](b).

## Ver también
- [[Multigrid en 1D - transferencia e implementacion]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
