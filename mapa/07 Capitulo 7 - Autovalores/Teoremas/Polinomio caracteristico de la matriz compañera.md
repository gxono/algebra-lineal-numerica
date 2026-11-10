---
tags: [capitulo7, teorema]
---

# El polinomio característico de la matriz compañera $C_q$ es $q$

Sea $q(t)=c_0+c_1t+\cdots+c_{n-1}t^{n-1}+t^n$ y $C_q$ su matriz compañera (primera fila $-c_{n-1},\dots,-c_0$; unos en la subdiagonal; ceros en el resto). Entonces $\det(tI-C_q)=q(t)$.

> [!demostracion]- (desarrollo por la primera fila)
> $$tI-C_q=\begin{pmatrix} t+c_{n-1}&c_{n-2}&\cdots&c_1&c_0\\ -1&t&&&\\ &-1&t&&\\ &&\ddots&\ddots&\\ &&&-1&t \end{pmatrix}.$$
> Desarrollamos por la primera fila. Al eliminar la fila $1$ y la columna $j$, el menor es **diagonal por bloques**:
> - filas $2..j$, columnas $1..j-1$: bidiagonal superior con $-1$ en la diagonal ⟹ determinante $(-1)^{j-1}$;
> - filas $j+1..n$, columnas $j+1..n$: bidiagonal inferior con $t$ en la diagonal ⟹ determinante $t^{n-j}$;
> - los bloques fuera de la diagonal son nulos (la fila $i\leq j$ solo tiene entradas en las columnas $i-1,i\leq j$; la fila $i>j$ solo en $i-1,i\geq j$, y la columna $j$ fue eliminada).
>
> Con el signo del cofactor $(-1)^{1+j}$, el término $j$ es $a_{1j}(-1)^{1+j}(-1)^{j-1}t^{n-j}=a_{1j}\,t^{n-j}$. Sumando:
> $$\det(tI-C_q)=(t+c_{n-1})t^{n-1}+\sum_{j=2}^{n}c_{n-j}t^{n-j}=t^n+c_{n-1}t^{n-1}+\cdots+c_1t+c_0=q(t).$$

Consecuencia: calcular autovalores de matrices generales es al menos tan difícil como hallar raíces de polinomios generales, y por Abel no hay método directo para $n\geq5$.

```julia
c = [-6.0, 11.0, -6.0]                       # q(t) = t³ - 6t² + 11t - 6 = (t-1)(t-2)(t-3)
n = length(c)
Cq = [reshape(-reverse(c), 1, n); I(n-1) zeros(n-1)]
eigvals(Cq)                                  # [1.0, 2.0, 3.0]
```

## Ver también
- [[Autovalores - repaso y por que los metodos son iterativos]]
- [[Demostraciones - Capitulo 7|Demostraciones]]
