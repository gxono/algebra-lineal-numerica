---
tags: [moc, julia]
---

# Referencia de Julia

Qué función de `LinearAlgebra` (o de otro paquete estándar) llamar para cada tarea, organizado por lo que querés calcular - no alfabético. Para cada una: qué hace Julia "gratis" (ya implementado, optimizado, usar en la práctica) vs. qué implementa el vault a mano (para entender el algoritmo por dentro, ver la nota de teoría correspondiente). Todo lo de acá es `using LinearAlgebra` salvo que se indique otro paquete.

## Normas

```julia
norm(v)        # norma 2 de un vector (default)
norm(v, 1)     # norma 1
norm(v, Inf)   # norma infinito
opnorm(A, 1)   # norma matricial inducida 1 (máxima suma de columnas)
opnorm(A, 2)   # norma matricial inducida 2 (mayor valor singular) - la única sin fórmula cerrada simple
opnorm(A, Inf) # norma matricial inducida infinito (máxima suma de filas)
```
`norm` es para **vectores**; `opnorm` es para la norma **matricial inducida** - son funciones distintas, no al revés (`norm(A)` existe pero da la norma de Frobenius, *no* ninguna de las tres inducidas del curso). Ver [[Normas de vectores y matrices]], [[Glosario - Indice|Glosario#Norma inducida]].

## Autovalores, autovectores, radio espectral

```julia
eigvals(A)      # todos los autovalores (vector)
eigvecs(A)      # todos los autovectores (matriz, columnas)
eigen(A)        # ambos juntos, como un objeto (eigen(A).values, eigen(A).vectors)
eigmin(A)       # autovalor mínimo (solo matrices simétricas/hermíticas - más rápido que eigvals(A)[1])
eigmax(A)       # autovalor máximo (ídem)
maximum(abs.(eigvals(A)))   # radio espectral ρ(A) - Julia no tiene una función `rho`/`spectralradius` directa
```
Ver [[Espectro y radio espectral]]. `eigmin`/`eigmax` funcionan en cualquier matriz **mientras sus autovalores sean reales** (lanzan `DomainError` si alguno sale complejo) - en la práctica, usarlas para matrices simétricas/sdp (donde los autovalores son reales por teoría) y no para matrices genéricas no simétricas.

## Número de condición

```julia
cond(A)         # κ₂(A) = σ₁/σₘ (norma 2 por default, sirve para A rectangular también)
cond(A, 1)      # κ₁(A)
cond(A, Inf)    # κ∞(A)
```
Ver [[Estabilidad y numero de condicion]], [[Numero de condicion y error relativo]]. Para $A$ sdp, $\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$ (más barato que pasar por valores singulares): `eigmax(A)/eigmin(A)`.

## Resolver sistemas lineales

```julia
x = A \ b             # el método a usar SIEMPRE que se pueda - Julia elige el algoritmo según el tipo de A
x = inv(A) * b        # NUNCA - más lento y menos estable, ver [[Resolver sistemas lineales en Julia]]
```
`A \ b` despacha automáticamente: triangular → sustitución, tridiagonal/banda → algoritmo especializado $O(n)$, simétrica/sdp → Cholesky, densa genérica → LU con pivoteo, rectangular → mínimos cuadrados vía QR. Ver [[2.10 - Funcion lu de Julia|Problema 2.10]] para la comparación directa con la implementación a mano (`eleupiv`).

## Factorizaciones: built-in vs. implementación propia del vault

| Factorización | Julia (LAPACK, usar en la práctica) | Implementación del vault (para entender) |
|---|---|---|
| LU sin pivoteo | - (Julia siempre pivotea) | [[Factorizacion LU\|`eleu`]] |
| LU con pivoteo | `lu(A)` → `.L`, `.U`, `.p`, `.P` | [[Pivoteo y matrices de permutacion\|`eleupiv`]] |
| Cholesky | `cholesky(A)` (sobre `Matrix` denso usa el triángulo superior por default; envolver en `Symmetric(A)`/`Hermitian(A)` si hace falta ser explícito) | [[Factorizacion de Cholesky\|`mi_cholesky!`]] |
| QR | `qr(A)` → `.Q`, `.R` (Householder por detrás) | [[Matrices de Householder\|`qrhouse`/`qrsolve`]] |
| QR (Gram-Schmidt) | *(Julia no expone Gram-Schmidt directamente; `qr` siempre usa Householder)* | [[Ortogonalizacion de Gram-Schmidt\|`gram_schmidt!`]] |

`lu(A)`/`cholesky(A)`/`qr(A)` devuelven un **objeto factorización** (no una tupla de matrices sueltas): `F \ b` resuelve usando esa factorización directamente, sin rehacerla.

## Tipos de matrices estructuradas (y por qué importan)

Usar el tipo correcto no es solo prolijidad - cambia **qué algoritmo corre por detrás** (dispatch de Julia):

```julia
Diagonal(v)                              # matriz diagonal a partir de un vector - O(n) en todo
SymTridiagonal(diag, offdiag)            # tridiagonal simétrica - eigvals/eigmin/eigmax/\ en O(n) o casi
Tridiagonal(sub, diag, sup)              # tridiagonal general (no necesariamente simétrica)
UpperTriangular(A)  # / LowerTriangular(A)  # fuerza sustitución hacia atrás/adelante en vez de LU completo
Symmetric(A)                             # le dice a Julia "usá el algoritmo para simétricas" (ej. habilita cholesky(A))
```
Si se arma `A` como matriz densa genérica cuando en realidad es tridiagonal, `A \ b` hace LU denso $O(n^3)$ en vez del $O(n)$ que el tipo `SymTridiagonal`/`Tridiagonal` habilita - la diferencia de performance es enorme para $n$ grande (ver [[Costo computacional metodos directos]]). Ejemplo usado en todo el Capítulo 5: `A = SymTridiagonal(fill(2/h^2, N), fill(-1/h^2, N-1))`.

### `UpperTriangular`/`LowerTriangular`/`Diagonal` son *wrappers*, no copias con ceros

`UpperTriangular(A)` **no** copia `A` ni pone ceros en la memoria de abajo de la diagonal - `parent(UpperTriangular(A)) === A` sigue siendo la matriz original, con lo que sea que tenga guardado en esa mitad. Lo que SÍ garantiza el wrapper es que las operaciones especializadas (`\`, `*`, `det`, `inv`, etc.) **nunca leen** las entradas fuera del triángulo correspondiente - no es que las traten como cero y por accidente alguna se cuele: el bucle del algoritmo de sustitución directamente no itera sobre esos índices. `Diagonal(A)` es distinto: ahí sí se **extrae** (copia) solo la diagonal a un vector nuevo, así que cualquier otra entrada de `A` queda completamente descartada.

**Verificado** con una matriz a propósito no simétrica y con basura en la mitad que debería ignorarse:
```julia
A = [4.0 99.0 99.0;
     1.0 5.0  99.0;
     2.0 3.0  6.0]
v = [1.0, 2.0, 3.0]

LowerTriangular(A) \ v   # da exactamente lo mismo que resolver con [4 0 0; 1 5 0; 2 3 6] \ v
UpperTriangular(A) \ v   # da exactamente lo mismo que resolver con [4 99 99; 0 5 99; 0 0 6] \ v
```
Los `99.0` nunca entran en la cuenta en ninguno de los dos casos - confirmado comparando contra la versión con ceros explícitos a mano, coinciden a precisión de máquina.

**Por qué importa:** esto es lo que hace correcto a `sgs_precon` ([[Precondicionadores genericos]]) sin tener que separar `A` en dos matrices explícitas - `LowerTriangular(A)` ya es exactamente $L+D$ y `UpperTriangular(A)` ya es exactamente $L^T+D$, leyendo directamente de la misma `A` simétrica. Y es la razón por la que estos wrappers son baratos: no hay copia de memoria, solo cambia qué algoritmo ($O(n^2)$ de sustitución en vez de $O(n^3)$ de LU genérico) se despacha.

## Construcción de matrices y vectores

```julia
zeros(n)               # vector de ceros
zeros(n,m)             # matriz de ceros
ones(n)                # vector de unos
ones(n,m)              # matriz de unos
fill(v, n)             # vector con el valor v repetido n veces
range(a, b, length=n)  # n puntos equiespaciados entre a y b (no es un Vector hasta que se usa collect o se indexa)
I                      # identidad "abstracta" (UniformScaling) - se adapta al tamaño que necesite la operación
diagm(0 => v)                       # matriz diagonal cuadrada a partir del vector v
diagm(0 => d, 1 => s, -1 => s)      # con sobre/sub-diagonales (offsets)
diagm(m, n, k => v)                 # versión RECTANGULAR m×n, con offset k (usada en el parcial modelo para A(ε))
```

## Operaciones básicas

```julia
A'          # transpuesta (o conjugada-transpuesta si A es compleja)
dot(x, y)   # x ⋅ y, producto escalar (equivalente a x'y para vectores reales)
det(A)      # determinante - caro (O(n³)) y casi nunca hace falta calcularlo para resolver sistemas
rank(A)     # rango numérico (vía SVD internamente)
isposdef(A) # true si A es simétrica definida positiva (intenta Cholesky por detrás)
```

## Benchmarking

```julia
using BenchmarkTools
@benchmark A \ b        # corre muchas veces, da estadísticas (min/mediana/media, memoria)
@time A \ b             # una sola corrida - más rápido de escribir, menos confiable (incluye compilación la primera vez)
```
Ver [[2.2 - Benchmark sustitucion|Problema 2.2]].

## Gráficos

Ver [[Graficos con CairoMakie]] (curvas, leyendas, escala log, curvas de nivel).

## IterativeSolvers.jl y Preconditioners.jl: métodos y precondicionadores de fábrica

Paquetes externos (no vienen con `LinearAlgebra`) con implementaciones de todo lo del Capítulo 5 - útiles para verificar la implementación propia contra algo independiente, o para resolver sistemas grandes sin reescribir `gradconj`/`gmres` a mano.

### `IterativeSolvers.jl`

```julia
using IterativeSolvers

x, hist = cg(A, b; Pl=precon, reltol=1e-10, log=true)       # gradiente conjugado (A sdp)
x, hist = gmres(A, b; restart=n, reltol=1e-10, log=true)    # GMRes
hist.iters    # cantidad de iteraciones (solo si log=true)
```
`cg`/`gmres` iteran **hasta cumplir la tolerancia** (`reltol`, `abstol`), igual que `gradconj`/`gmres` del vault ([[Implementacion del gradiente conjugado]], [[GMRes - implementacion y rotaciones de Givens]]). `Pl` es el precondicionador a izquierda (ver más abajo); para GMRes sin precondicionar y comparar contra la teoría del curso, dejar `Pl` por default.

> [!pendiente]
> **`restart` de `gmres` tiene default `min(20, n)`** - es decir, por default `gmres(A,b)` es en realidad `GMRes(20)` **con restart**, no el GMRes "completo" del curso ([[GMRes - implementacion y rotaciones de Givens|Restart: GMRes(m)]]). Sin pasar `restart=n` explícitamente, puede no converger a la tolerancia pedida en el número de iteraciones esperado - *verificado*: en el problema modelo 1D ($n=50$), `gmres(Ah,b)` default se queda con residuo relativo $\approx0.2$ tras 50 "iteraciones" (en realidad, $\sim2$ ciclos de reinicio de $20$), mientras que `gmres(Ah,b;restart=n)` converge en 25 iteraciones, igual que `cg`. Coincide exactamente con lo que dice la teoría: GMRes(m) no tiene garantía de convergencia.

```julia
x = jacobi(A, b; maxiter=500)        # NO iteran hasta converger -
x = gauss_seidel(A, b; maxiter=500)  # ejecutan EXACTAMENTE `maxiter` iteraciones (default 10)
x = sor(A, b, ω; maxiter=500)        # bajo su propio riesgo: si el método diverge
x = ssor(A, b, ω; maxiter=500)       # con ese ω, `x` sale sin ningún aviso
```

> [!pendiente]
> **Pitfall real, verificado numéricamente:** `jacobi`/`gauss_seidel`/`sor`/`ssor` de este paquete no tienen criterio de parada por tolerancia - hacen la cantidad de iteraciones que se les pida, converjan o no, y no lanzan ningún error ni warning si diverge. Además, **`ssor` de este paquete parece tener un bug**: en el problema modelo 1D ($A_h$ sdp tridiagonal, donde SOR/SSOR con $\omega\in(0,2)$ debería converger por el análogo de Ostrowski-Reich), se probó `ssor(Ah,b,ω)` para $\omega\in\{0.5,1.0,1.2,1.5,1.8,1.9,1.95\}$: para $\omega\gtrsim1.3$ el error decrece lentamente las primeras ~20 iteraciones y luego **diverge catastróficamente** (de $\|x-x_*\|\approx1$ a $10^{12}$ en 200 iteraciones) - un comportamiento no lineal incompatible con una iteración estacionaria genuina (que debería decrecer o crecer monótonamente como $\rho^k$). `sor` (solo el barrido hacia adelante, sin el paso hacia atrás) con los mismos $\omega$ se comportó correctamente. **Conclusión: no confiar en `ssor` de este paquete sin verificar contra la solución exacta**; preferir `sor`, o el precondicionador `sgs_precon` propio del vault ([[Precondicionadores genericos]]), que fue verificado y da la identidad correcta ([[Identidad del precondicionador Gauss-Seidel simetrico]]).

### `Preconditioners.jl`

Construyen un objeto precondicionador para pasarle a `Pl=` en `cg`/`gmres` de `IterativeSolvers.jl` (acá la convención es $M\approx A^{-1}$, igual que en [[Precondicionamiento y GCP]] - no la convención de `pcg` de matlab/octave mencionada en las diapositivas, donde $M\approx A$).

```julia
using Preconditioners

DiagonalPreconditioner(A)                    # M = D⁻¹ (equivalente al Jacobi de Precondicionadores genericos)
CholeskyPreconditioner(A, memory=2)          # Cholesky incompleta; `memory` controla cuánto fill-in se permite retener
AMGPreconditioner(RugeStuben, A)             # multigrid algebraico (AMG) como precondicionador
AMGPreconditioner(SmoothedAggregation, A)    # otra variante de AMG, agregación suavizada

x, hist = cg(A, b; Pl=AMGPreconditioner(RugeStuben, A), reltol=1e-10, log=true)
```
`AMGPreconditioner` es la versión "de fábrica" de la idea de [[Multigrid en 1D - transferencia e implementacion|multigrid]], pero **algebraica**: no necesita saber que $A$ viene de discretizar una ecuación diferencial (a diferencia del multigrid geométrico del vault, que usa explícitamente la malla) - construye la jerarquía de mallas gruesas mirando solo los coeficientes de $A$. `RugeStuben` y `SmoothedAggregation` son dos estrategias distintas para elegir esa jerarquía.

**Verificado** (laplaciano 2D, $N=20$, $n=400$ - mismo problema que [[Precondicionadores genericos]]), iteraciones de `cg` hasta $\|r_k\|_2/\|b\|_2<10^{-10}$:

> [!comparacion]
> | Precondicionador | Iteraciones |
> |---|---|
> | Ninguno | 41 |
> | `DiagonalPreconditioner` | 41 (sin cambio - diagonal constante, mismo fenómeno que en [[Precondicionadores genericos]]) |
> | `CholeskyPreconditioner(A, 0)` | 41 (memoria 0: casi no retiene fill-in, casi no ayuda) |
> | `CholeskyPreconditioner(A, 2)` | 30 |
> | `AMGPreconditioner(RugeStuben, A)` | 6 |
> | `AMGPreconditioner(SmoothedAggregation, A)` | 8 |

Confirma, con una implementación independiente de la del vault, el mismo patrón ya documentado en [[Tablas Comparativas - Precondicionadores]]: los precondicionadores genéricos (diagonal, Cholesky incompleta) mejoran poco; el que usa la estructura (AMG, análogo algebraico de multigrid) gana por mucho.

## Sistemas no lineales (Capítulo 6)

```julia
# Newton a mano: el jacobiano a mano o con Symbolics
using Symbolics
@variables x y
J = Symbolics.jacobian([x^2 + y^2 - 4, x^2 - y^2 - 1], [x, y])   # Num[2x 2y; 2x -2y]
DF = eval(build_function(J, x, y)[1])                            # función numérica DF(x, y)

h = DF(x0...) \ F(x0...)          # paso de Newton: SIEMPRE \, nunca inv(...)
F0 = lu(DF(x0...)); F0 \ F(xk...) # Newton simplificado: factorizar una vez

# Minimizar sin derivadas (equivalente de fminsearch = Nelder–Mead)
using Optim                       # no viene con Julia: ] add Optim
r = optimize(h, x0, NelderMead(), Optim.Options(g_tol=1e-14))
Optim.minimizer(r), Optim.minimum(r), Optim.f_calls(r)

# Curvas de nivel cero (sin meshgrid)
contour!(ax, xs, ys, [f(x, y) for x in xs, y in ys], levels=[0])

# Orden de convergencia experimental, y precisión extendida para verlo
orden(e) = [log(e[k+1]/e[k]) / log(e[k]/e[k-1]) for k in 2:length(e)-1]
setprecision(BigFloat, 500); x = BigFloat[16//10, 12//10]
```
Ver [[Metodo de Newton - definicion y buena definicion]], [[6.9 - Programar Newton]], [[6.11 - Cuatro ceros y fminsearch]].

> [!pendiente] `optimize` con `NelderMead()`: el símplex inicial por defecto es grande (≈ $0.5\,x_0$ por coordenada) y puede saltar a otro mínimo. Ver [[6.11 - Cuatro ceros y fminsearch|Problema 6.11]].

## Autovalores (Capítulo 7)

```julia
eigvals(A); eigen(A)              # el QR práctico de LAPACK (Hessenberg + shifts + deflación)
eigvals(Symmetric(A))             # simétrica: reales, ascendentes, más rápido
S = schur(A); S.T, S.Z            # forma de Schur: A = Z T Zᴴ (T cuasi-triangular si A es real)
Hh = hessenberg(A); Hh.H, Hh.Q    # reducción a Hessenberg: A = Q H Qᴴ
F = qr(A); Matrix(F.Q), F.R       # F.Q es compacta: Matrix(F.Q) para la Q explícita
rank(A - λ*I)                     # multiplicidad geométrica = n - rank
nullspace(A - λ*I)                # base del autoespacio
dot(u, v)                         # = uᴴv (conjuga el primero): cociente de Rayleigh dot(u, A*u)
```
Lo que el vault implementa a mano, para entender: [[Metodo de las potencias|potencias]], [[Metodo de las potencias inversas|potencias inversas]] (con `lu` una vez), [[Algoritmo QR - definicion y convergencia|QR básico]], [[Algoritmo QR - Hessenberg, shifts y deflacion|QR con shift y deflación y reducción a Hessenberg]].

## Lo mínimo de Symbolics.jl

Para manipular expresiones simbólicamente en vez de numéricamente (ej. resolver las ecuaciones normales de la [[Parcial Modelo - Problema 2|matriz de Läuchli]] "en función de $\varepsilon$", como pide el parcial modelo).

```julia
using Symbolics

@variables eps          # declara eps como símbolo (tipo Num), no como número
A = [1 1 1; eps 0 0; 0 eps 0; 0 0 eps]     # matriz simbólica: tipo Matrix{Num}
b = [1, 0, 0, 0]

AtA = A'A                        # las operaciones matriciales normales funcionan sobre símbolos
Atb = A'b

x = Symbolics.simplify.(AtA \ Atb)    # SÍ, `\` funciona sobre matrices simbólicas - resuelve el sistema simbólico
# x = Num[1 / (3 + eps^2), 1 / (3 + eps^2), ...]  (verificado - puede no quedar simplificado igual en cada entrada,
# aplicar simplify.() de nuevo o comparar numéricamente si hace falta confirmar que son iguales)

xval = Symbolics.substitute(x, Dict(eps => 0.01))   # sustituir un valor numérico
Symbolics.value.(xval)                              # "desenvolver" a Float64 de verdad
```
Puntos clave:
- `@variables nombre` crea un símbolo de tipo `Num`; se pueden declarar varios juntos: `@variables x y z`.
- Una matriz/vector con entradas `Num` se comporta como cualquier `Matrix`/`Vector` de Julia: `+`, `*`, `'`, `\` funcionan igual, solo que el resultado queda en forma simbólica.
- `Symbolics.simplify.(expr)` (con el punto, porque `expr` suele ser un array) intenta simplificar cada entrada - no siempre llega a la forma más corta posible, pero sí a una matemáticamente correcta.
- `Symbolics.substitute(expr, Dict(var => valor))` evalúa la expresión en un valor concreto, sin "perder" el símbolo - el resultado sigue siendo tipo `Num` hasta que se llama `Symbolics.value.(...)` para pasarlo a `Float64`.
- No hace falta mucho más que esto para los usos del curso: declarar, operar, resolver con `\`, sustituir un valor, listo.

## Ver también
- [[Sintaxis basica de Julia]]
- [[Resolver sistemas lineales en Julia]]
- [[Glosario - Indice|Glosario]]
- [[Formulario - Indice|Formulario]]
- [[Inicio]]
