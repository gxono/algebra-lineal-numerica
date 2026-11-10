---
tags: [capitulo5, teorema]
---

# El operador grueso $A_{2h}$ coincide con el operador de Galerkin

## Enunciado

La matriz $A_{2h}$ del problema modelo en la malla gruesa (la misma familia tridiagonal $\frac1{(2h)^2}\text{tridiag}(-1,2,-1)$, pero con paso $2h$) coincide **exactamente** con el *operador de Galerkin*, construido en abstracto a partir de $A_h$ y los operadores de transferencia:
$$
A_{2h} = I_h^{2h}\,A_h\,I_{2h}^h.
$$
Es decir: definir el problema grueso "a mano" (discretizando de nuevo, con paso $2h$) y definirlo "sandwicheando" $A_h$ entre restricción y prolongación dan **la misma matriz** - no es una coincidencia del problema modelo 1D en particular, es una propiedad general de por qué esta elección de $A_{2h}$ es la "correcta" (el operador de Galerkin es la construcción estándar en multigrid cuando no hay una discretización directa obvia en la malla gruesa, por ejemplo con coeficientes variables o en geometrías complicadas).

> [!demostracion]-
> Alcanza con verificarlo columna por columna: mostrar que $A_h$ aplicado a la $k$-ésima columna de $I_{2h}^h$ (el $k$-ésimo vector canónico prolongado), restringido de nuevo con $I_h^{2h}$, da la $k$-ésima columna de $A_{2h}$.
>
> **Prolongar el $k$-ésimo canónico.** Sea $v:=I_{2h}^he_k$ (con $e_k$ el $k$-ésimo canónico de $\mathbb{R}^{n_c}$). Por la fórmula de prolongación ([[Multigrid en 1D - transferencia e implementacion|eq. 5.22]]), $v$ tiene soporte en tres índices finos consecutivos alrededor de $2k$:
> $$v_{2k}=1, \qquad v_{2k-1}=v_{2k+1}=\tfrac12, \qquad v_i=0 \text{ en el resto}.$$
>
> **Aplicar $A_h$.** Como $A_h=\frac1{h^2}\text{tridiag}(-1,2,-1)$, $(A_hv)_i=\frac1{h^2}(2v_i-v_{i-1}-v_{i+1})$. Evaluando en cada índice donde puede ser no nulo (el tridiagonal "ensancha" el soporte de $v$ en uno para cada lado):
> $$(A_hv)_{2k} = \frac1{h^2}(2\cdot1-\tfrac12-\tfrac12) = \frac1{h^2}, \qquad (A_hv)_{2k\pm1} = \frac1{h^2}(2\cdot\tfrac12-0-1) = 0,$$
> $$(A_hv)_{2k\pm2} = \frac1{h^2}(2\cdot0-0-\tfrac12) = -\frac1{2h^2}, \qquad (A_hv)_i=0 \text{ fuera de } [2k-2,2k+2].$$
>
> **Restringir.** Por la fórmula de restricción (ponderación completa, $(I_h^{2h}w)_j=\frac14(w_{2j-1}+2w_{2j}+w_{2j+1})$), calculamos $(I_h^{2h}A_hv)_j$ para $j=k$, $j=k\pm1$ y $|j-k|\geq2$:
> $$(I_h^{2h}A_hv)_k = \tfrac14\big((A_hv)_{2k-1}+2(A_hv)_{2k}+(A_hv)_{2k+1}\big) = \tfrac14\Big(0+\tfrac2{h^2}+0\Big) = \frac1{2h^2},$$
> $$(I_h^{2h}A_hv)_{k+1} = \tfrac14\big((A_hv)_{2k+1}+2(A_hv)_{2k+2}+(A_hv)_{2k+3}\big) = \tfrac14\Big(0+2\big(-\tfrac1{2h^2}\big)+0\Big) = -\frac1{4h^2},$$
> y por simetría $(I_h^{2h}A_hv)_{k-1}=-\frac1{4h^2}$ también. Para $|j-k|\geq2$, los índices finos $2j-1,2j,2j+1$ caen fuera de $[2k-2,2k+2]$ (el soporte de $A_hv$), así que $(I_h^{2h}A_hv)_j=0$.
>
> **Conclusión.** La $k$-ésima columna de $I_h^{2h}A_hI_{2h}^h$ tiene diagonal $\frac1{2h^2}=\frac2{(2h)^2}$ y sub/sobre-diagonal $-\frac1{4h^2}=-\frac1{(2h)^2}$ - exactamente la $k$-ésima columna de $A_{2h}=\frac1{(2h)^2}\text{tridiag}(-1,2,-1)$. Como $k$ era arbitrario, $I_h^{2h}A_hI_{2h}^h=A_{2h}$.

Verificado numéricamente (`norm(Matrix(Ah(nc)) - Ih2h*Matrix(Ah(N))*I2hh)`, para $N=7,15,31$): diferencia $0.0$ exacta en los tres casos.

## Ver también
- [[Multigrid en 1D - transferencia e implementacion]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
