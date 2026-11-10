---
tags: [capitulo7, definicion, teorema]
---

# Autovalores - repaso y por qué los métodos tienen que ser iterativos (7.1)

> [!idea]
> Todo el capítulo trabaja en $\mathbb{C}^n$, con $A^H=\bar A^T$. Dos ideas para llevarse: (1) **no existe** un método directo (finito) para autovalores de matrices $n\geq5$: calcular autovalores es tan difícil como hallar raíces de polinomios (Abel), así que todo método es **iterativo**. (2) Las matrices **hermíticas** ($A^H=A$, en particular las reales simétricas) son el caso bueno: autovalores reales y base **ortonormal** de autovectores.

## Definiciones básicas

- $\lambda\in\mathbb{C}$ es **autovalor** de $A\in\mathbb{C}^{n\times n}$ si existe $x\neq0$ con $Ax=\lambda x$; $(\lambda,x)$ es un **autopar**.
- $\lambda$ autovalor $\iff A-\lambda I$ singular $\iff p(\lambda):=\det(\lambda I-A)=0$ (**ecuación característica**).
- $p(\lambda)=(\lambda-\lambda_1)\cdots(\lambda-\lambda_n)$: **una matriz $n\times n$ tiene exactamente $n$ autovalores contando multiplicidades**.
- Producto escalar en $\mathbb{C}^n$: $x\cdot y=y^Hx$.

## Por qué los métodos son iterativos

Dado un polinomio mónico $q(t)=c_0+c_1t+\cdots+c_{n-1}t^{n-1}+t^n$, su **matriz compañera**
$$
C_q=\begin{pmatrix}
-c_{n-1}&-c_{n-2}&\cdots&-c_1&-c_0\\
1&0&\cdots&0&0\\
0&1&\ddots&&\vdots\\
\vdots&&\ddots&0&0\\
0&\cdots&0&1&0
\end{pmatrix}
$$
tiene polinomio característico **$q$** ([[Polinomio caracteristico de la matriz compañera|demostración]]). Entonces un algoritmo que calcule autovalores de cualquier matriz resuelve cualquier ecuación polinomial. Por **Abel** (no hay fórmula por radicales para grado $\geq5$), no puede existir un método **directo** general: **todo método de autovalores es iterativo**. De hecho `roots` de MATLAB/Octave calcula los autovalores de $C_q$ (en Julia, `Polynomials.roots` hace lo mismo).

## Autovalores de $A^T$, $\bar A$, $A^H$; autovectores a izquierda

$\det(\lambda I-A)=0\iff\det(\lambda I-A^T)=0\iff\det(\bar\lambda I-A^H)=0$. Luego:
- $\sigma(A^T)=\sigma(A)$;
- $\sigma(A^H)=\overline{\sigma(A)}$ (los conjugados).

Si $A^Hy=\bar\lambda y$, conjugando y trasponiendo: $y^HA=\lambda y^H$. Ese $y^H$ es un **autovector a izquierda** de $A$ asociado a $\lambda$. Aparece en la hipótesis de convergencia del [[Metodo de las potencias|método de las potencias]] ($y_1^Hu_0\neq0$).

**Matrices reales:** el polinomio característico tiene coeficientes reales, así que los autovalores complejos vienen de a **pares conjugados**, y si $(\lambda,x)$ es autopar también lo es $(\bar\lambda,\bar x)$. Consecuencia: si el autovalor de módulo máximo de una matriz real es complejo, **no hay autovalor dominante** (su conjugado tiene el mismo módulo).

## Multiplicidades y matrices defectuosas

- **Multiplicidad algebraica:** multiplicidad como raíz de $p$.
- **Multiplicidad geométrica:** $\dim\{v:Av=\lambda v\}=n-\text{rank}(A-\lambda I)$.
- Siempre: geométrica $\leq$ algebraica. Si es $<$, el autovalor (y la matriz) es **defectuoso**.
- Defectuosa $\iff$ **no diagonalizable**. Ejemplo: $A=\begin{pmatrix}2&0\\1&2\end{pmatrix}$, $\lambda=2$ con multiplicidad algebraica 2 y geométrica 1.

**Forma de Jordan:** para toda $A$ existe $V$ no singular con $V^{-1}AV=\text{diag}(J_{m_1}(\lambda_1),\dots,J_{m_k}(\lambda_k))$, con bloques $J_m(\lambda)$ que tienen $\lambda$ en la diagonal y $1$ en la superdiagonal. Diagonalizable $\iff$ todos los bloques son $1\times1$.

Los problemas [[7.10 - QR con autovalor semisimple doble|7.10]] y [[7.11 - QR con autovalor defectuoso|7.11]] muestran qué le pasa al algoritmo QR en cada caso.

## Matrices hermíticas

$A^H=A$ (reales simétricas incluidas). Hechos (demostraciones en [[Autovalores y diagonalizacion de matrices hermiticas]]):
1. $a_{ii}\in\mathbb{R}$.
2. **Los autovalores son reales.** Si $Ax=\lambda x$: $\lambda\|x\|^2=x^HAx=(Ax)^Hx=\bar\lambda\|x\|^2$.
3. $U$ unitaria ⟹ $U^HAU$ hermítica.
4. **Teorema espectral:** existe $U$ unitaria con $U^HAU=\text{diag}(\lambda_1,\dots,\lambda_n)$, $\lambda_i\in\mathbb{R}$. Sale de [[Teorema 4.21|Schur]]: $T=U^HAU$ es triangular y hermítica, luego diagonal. Equivalente: **base ortonormal de autovectores**.

Más general: $A$ es **normal** si $AA^H=A^HA$; para normales la $T$ de Schur también es diagonal (pero los autovalores pueden ser complejos).

## En Julia

```julia
eigvals(A)                 # autovalores (orden: ascendente si A es Symmetric/Hermitian)
eigen(A)                   # .values, .vectors
eigvals(Symmetric(A))      # usa el algoritmo para simétricas: reales y más rápido
schur(A)                   # .T (triangular/cuasi-triangular), .Z (unitaria): A = Z T Zᴴ
rank(A - λ*I)              # → multiplicidad geométrica = n - rank
```

## Ver también
- [[Metodo de las potencias]] (siguiente tema)
- [[Espectro y radio espectral]] (Capítulo 4)
- [[Teorema 4.21]] (Schur)
- [[Capitulo 7 - MOC]]
