---
tags: [capitulo5, teorema]
---

# Teorema 5.22

Si $b$ es combinación lineal de $k$ autovectores de $A$, el GC partiendo de $x_0=0$ converge en a lo sumo $k$ pasos.

> [!demostracion]-
> Sea $b=\sum_{j\in J}c_jv_j$ con $|J|=k$, $Av_j=\lambda_jv_j$ (los $\lambda_j$ no necesariamente distintos entre sí). Con $x_0=0$: $x_*=A^{-1}b=\sum_{j\in J}\frac{c_j}{\lambda_j}v_j$, así que $x_0-x_*=-x_*$ es combinación lineal de **los mismos** $k$ autovectores $\{v_j\}_{j\in J}$.
>
> Tomamos $p(z):=\prod_{j\in J}\big(1-\frac z{\lambda_j}\big)\in\overline{\mathcal P}_k$ (grado $k$, $p(0)=1$; si algún $\lambda_j$ se repite el producto tiene factores iguales, pero sigue siendo un polinomio válido de grado $k$ con esas raíces). Como $p(\lambda_j)=0$ para cada $j\in J$:
> $$p(A)(x_0-x_*) = -\sum_{j\in J}\frac{c_j}{\lambda_j}p(A)v_j = -\sum_{j\in J}\frac{c_j}{\lambda_j}p(\lambda_j)v_j = 0.$$
> Por el [[Teorema 5.13]], $\|x_k-x_*\|_A=\min_{q\in\overline{\mathcal P}_k}\|q(A)(x_0-x_*)\|_A\leq\|p(A)(x_0-x_*)\|_A=0$, así que $x_k=x_*$.
>
> *(Sugerido como ejercicio en el apunte, usando el Teorema 5.13; es el mismo tipo de argumento que el [[Teorema 5.21]], pero aprovechando que $x_0-x_*$ vive en un subespacio invariante chico en vez de usar todos los autovalores distintos de $A$ - comparar con [[5.20 - Demostrar Corolario 5.41|el resultado análogo para GMRes, Corolario 5.41]].)*

> [!idea]
> Más fuerte que el anterior: si la solución "vive" en un subespacio invariante chico (aunque $A$ tenga muchos autovalores distintos), alcanza con $k$ pasos.

## Ver también
- [[GC como metodo directo]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
