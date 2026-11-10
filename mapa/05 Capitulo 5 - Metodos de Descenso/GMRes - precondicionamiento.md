---
tags: [capitulo5, teorema, algoritmo]
---

# [[GMRes - definicion y propiedad de minimizacion|GMRes]]: precondicionamiento (izquierda y derecha)

> [!idea]
> A diferencia de GCP (que necesitaba $S=M^{1/2}$ para mantener simetría, y por eso había que reescribir el algoritmo), GMRes no tiene esa restricción: $M$ puede ser **cualquier** aproximación de $A^{-1}$, y no hace falta modificar el algoritmo - alcanza con darle la rutina $ξ\to M(A\xi)$ (izquierda) o $ξ\to A(M\xi)$ (derecha).

## Teorema 5.42 (si $A$ está cerca de la identidad)

Si $\|I-A\|_2\leq\rho<1$: $\|r_k\|_2\leq\rho^k\|r_0\|_2$. Combinado con el [[Numero de condicion y error relativo|Lema 4.10]]: $\|e_k\|_2\leq\kappa_2(A)\rho^k\|e_0\|_2$. Demostración completa (eligiendo $p_k(z)=(1-z)^k$ y aplicando el Teorema 5.36): [[Teorema 5.42]].

## Del residuo al error, para una matriz cualquiera $S$ invertible

> [!demostracion]- Generalización directa de [[Lema 4.10]]
> Sea $S$ invertible y $e,\tilde e\neq0$ vectores con $Se=:\rho$, $S\tilde e=:\tilde\rho$ (en nuestro caso, $e=e_k=x_*-x_k$, $\tilde e=e_0$, $S=$ la matriz del sistema que efectivamente resuelve GMRes, y $\rho,\tilde\rho$ los residuos correspondientes - exactamente la misma situación que en [[Lema 4.10]], salvo que ahí se comparaba con una perturbación de $b$ y acá se comparan dos iterados). Igual que en esa demostración:
> $$\|e\| = \|S^{-1}\rho\| \leq \|S^{-1}\|\|\rho\|, \qquad \|\tilde\rho\| = \|S\tilde e\| \leq \|S\|\|\tilde e\| \implies \|\tilde e\| \geq \frac{\|\tilde\rho\|}{\|S\|}.$$
> Combinando:
> $$\frac{\|e\|}{\|\tilde e\|} \leq \frac{\|S^{-1}\|\|\rho\|}{\|\tilde\rho\|/\|S\|} = \|S\|\|S^{-1}\|\frac{\|\rho\|}{\|\tilde\rho\|} = \kappa(S)\frac{\|\rho\|}{\|\tilde\rho\|}.$$

## Precondicionador a izquierda: $MAx=Mb$

Si $\|I-MA\|_2=\rho_I<1$: el residuo que **ve** GMRes es el precondicionado, $M(b-Ax_k)=MAe_k$. Aplicando la desigualdad de arriba con $S=MA$, $e=e_k$, $\tilde e=e_0$ (y $\rho=MAe_k$, $\tilde\rho=MAe_0$, los residuos precondicionados en los pasos $k$ y $0$):
$$
\frac{\|e_k\|_2}{\|e_0\|_2} \leq \kappa_2(MA)\rho_I^k.
$$
**Ventaja:** si $κ₂(MA)≪κ₂(A)$, el residuo precondicionado se parece más al error relativo verdadero - el criterio de parada es más confiable.

## Precondicionador a derecha: $(AM)y=b$, $x=My$

Si $\|I-AM\|_2=\rho_D<1$: el residuo que ve GMRes es el **verdadero**, $b-Ax_k=Ae_k$. Aplicando la misma desigualdad con $S=A$:
$$
\frac{\|r_k\|_2}{\|r_0\|_2}\leq\rho_D^k, \qquad \frac{\|e_k\|_2}{\|e_0\|_2}\leq\kappa_2(A)\rho_D^k.
$$
**Ventaja:** se minimiza el residuo verdadero (no uno "distorsionado" por $M$), pero la constante de la cota de error sigue siendo $κ₂(A)$ (no mejora como a izquierda).

## Tabla comparativa

> [!comparacion]
> | | A izquierda: $MAx=Mb$ | A derecha: $AMy=b$, $x=My$ |
> |---|---|---|
> | Condición | $\|I-MA\|_2=\rho_I<1$ | $\|I-AM\|_2=\rho_D<1$ |
> | Qué minimiza GMRes | $\|M(b-Ax)\|_2$ (precondicionado) | $\|b-Ax\|_2$ (verdadero) |
> | Cota del error relativo | $\leq\kappa_2(MA)\rho_I^k$ | $\leq\kappa_2(A)\rho_D^k$ |
> | Paso final | ninguno | $x=My$ |

En ambos casos, **no hace falta modificar el algoritmo de GMRes**, a diferencia del GCP donde hubo que reescribirlo para que $S=M^{1/2}$ desapareciera.

## Cómo aplicarlo con el `gmres` ya implementado

`gmres` ([[GMRes - implementacion y rotaciones de Givens]]) acepta tanto una matriz como un *function handle* para $A$ - precondicionar es simplemente envolver esa acción:

```julia
# A, b, x0, M son los del problema concreto; ejemplo con el problema modelo 1D (ver [[Multigrid en 1D - problema modelo y suavizado]])
n = 50
A = Ah(n)
xstar = ones(n)
b = A*xstar
x0 = zeros(n)
M = Diagonal(1 ./ diag(A))     # precondicionador de ejemplo (diagonal); M puede ser matriz o function handle

Av(v) = A isa Function ? A(v) : A*v
Mv(v) = M isa Function ? M(v) : M*v     # acción de M, matrix-free si hace falta

# izquierda: resolver MAx=Mb
x, k = gmres(v -> Mv(Av(v)), Mv(b), x0)

# derecha: resolver (AM)y=b, con x=My
y, k = gmres(v -> Av(Mv(v)), b, x0)
x = Mv(y)
```
(Con un precondicionador diagonal, que para $A_h$ no cambia mucho el acondicionamiento, la cantidad de iteraciones queda igual que sin precondicionar - para ver el efecto real de un buen precondicionador sobre la convergencia, ver los datos de [[Precondicionadores genericos]] y [[5.11 - GCP y Cholesky incompleta|Problema 5.11]].)

## Ver también
- [[GMRes - definicion y propiedad de minimizacion]]
- [[Precondicionamiento y GCP]] (el contraste)
