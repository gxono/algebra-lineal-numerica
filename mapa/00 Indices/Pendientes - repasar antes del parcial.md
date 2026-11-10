---
tags: [moc]
---

# Pendientes - repasar con cuidado antes del parcial

Todos los `> [!pendiente]` del vault, juntos en un solo lugar. No es contenido faltante - es contenido ya completo, señalado porque es un error típico, una sutileza de implementación, o un punto que vale la pena repasar paso a paso.

## Capítulo 1

> [!pendiente] $A^2$ vs. $A.^2$
> $A^2$ es producto matricial ($A*A$); $A.^2$ es cada entrada al cuadrado. Error típico de parcial. Ver [[Sintaxis basica de Julia]].

> [!pendiente] `diag` en Julia vs. Octave
> En Julia, `diag` **solo extrae** la diagonal; para construir una matriz diagonal hace falta `Diagonal(v)` o `diagm(v)` (en Octave, `diag` hace ambas cosas). Ver [[Sintaxis basica de Julia]].

> [!pendiente] Nunca calcular `inv(A)` para resolver un sistema
> `A \ b` es más rápido (no hace falta la inversa completa) y más estable numéricamente. Ver [[Resolver sistemas lineales en Julia]].

> [!pendiente] `lines!` de CairoMakie no acepta `marker`
> Para puntos marcados sobre una curva, usar `scatterlines!` o `scatter!` en su lugar. Ver [[Graficos con CairoMakie]].

## Capítulo 2

> [!pendiente] No siempre existe [[Factorizacion LU|factorización LU]]
> Contraejemplo: $A=\begin{bmatrix}0&1\\1&0\end{bmatrix}$. La hipótesis de menores principales no singulares es necesaria - por eso hace falta pivotear. Ver [[Factorizacion LU]], [[2.7 - Matriz sin factorizacion LU|Problema 2.7]].

> [!pendiente] Ninguna implementación manual le gana a `L \ b`
> `\` despacha a LAPACK/BLAS (multi-hilo, bloqueo por caché). Implementar a mano es para *entender*, no para producción. Ver [[Matrices triangulares - sustitucion]].

> [!pendiente] Julia vs. MATLAB/Octave: bucles explícitos
> En Julia (JIT), un bucle bien tipado con `@inbounds`/`@simd`/`@views` puede competir con código vectorizado - a diferencia de MATLAB/Octave, donde vectorizar es mucho más determinante. Ver [[2.2 - Benchmark sustitucion|Problema 2.2]].

> [!pendiente] `trisup` a partir de `triinf` no es lo más eficiente
> Invertir índices (`U[end:-1:1,end:-1:1]`) genera una vista no siempre amigable con la caché. Buen truco para no duplicar código, no para producción. Ver [[2.4 - trisup a partir de triinf|Problema 2.4]].

## Capítulo 3

> [!pendiente] Ecuaciones normales: más barato, menos estable
> $\cong nm^2$ (mitad que QR) pero $\kappa_2(A^TA)=\kappa_2(A)^2$ - si $A$ está mal condicionada, se pierde el doble de cifras significativas. Ver [[Minimos cuadrados y ecuaciones normales]], [[Estabilidad y numero de condicion]].

> [!pendiente] [[Ortogonalizacion de Gram-Schmidt|Gram-Schmidt]] modificado vs. clásico
> El modificado resta las proyecciones de a una, recalculando con la columna ya corregida - más estable que el clásico (que proyecta sobre la columna original sin corregir). Mismo fenómeno que el Gram-Schmidt modificado dentro de [[GMRes - definicion y propiedad de minimizacion|GMRes]] (Capítulo 5). Ver [[Ortogonalizacion de Gram-Schmidt]].

## Capítulo 4

> [!pendiente] $\rho(A)\leq\|A\|$ puede ser estricta
> $\rho(A)$ puede valer $0$ aunque $A\neq0$ (ej.: triangular con diagonal nula). Ver [[Espectro y radio espectral]], [[4.7 - Radio espectral y norma|Problema 4.7]].

> [!pendiente] La demostración **clásica** de convergencia de [[Metodo de Gauss-Seidel|Gauss-Seidel]] (vía e.d.d.) es más sutil que la de [[Metodo de Jacobi|Jacobi]]
> La prueba de $\rho(M_{GS})<1$ no sale de acotar una norma de $M_{GS}$ directamente - hace falta trabajar con el autovector y la fila $k$ donde se alcanza el máximo ([[Teorema de convergencia de Gauss-Seidel con e.d.d.]]). **Pero** sí existe una cota cuantitativa directa de $\|M_{GS}\|_\infty$ (con un argumento distinto, inducción fila por fila): da la **misma cota** $\alpha$ que Jacobi - ver [[Cota de norma infinito de Gauss-Seidel para matrices e.d.d.]], usada en [[Parcial Modelo - Problema 4]](a)/(b). Repasar ambas: la cualitativa (por qué converge) y la cuantitativa (qué tan rápido).

> [!pendiente] Dos bugs típicos implementando Gauss-Seidel
> (1) `tril(A)`/`triu(A)` en Julia **incluyen** la diagonal - hace falta `tril(A,-1)`/`triu(A,1)`. (2) Si $M_{gs}=-L Dinv*U$ ya tiene el signo puesto, la iteración es $x_n=M_{gs}x+LDinv\cdot b$ **sin** volver a negar. Ver [[Metodo de Gauss-Seidel]].

> [!pendiente] En Jacobi, `A*x` ya incluye el término diagonal
> No debe estar en la suma $\sum_{j\neq i}a_{ij}x_j$ - olvidar sacarlo da una iteración distinta a Jacobi. Ver [[Metodo de Jacobi]].

> [!pendiente] Otros casos donde [[Metodo de Richardson|Richardson]] no tiene $\alpha$ que funcione
> $A$ singular (autovalor $0$: $1-\alpha\cdot0=1$ siempre, nunca $<1$); o $A$ con autovalores puramente imaginarios, ej. $A=\begin{bmatrix}0&-1\\1&0\end{bmatrix}$, donde $|1-\alpha i|\geq1$ para todo $\alpha$ real. Ver [[4.16 - Matriz sin alpha de Richardson que converja|Problema 4.16]].

## Capítulo 5

> [!pendiente] El Teorema 5.34 (convergencia del ciclo en V) no se demuestra
> Es el único resultado del curso que no se demuestra - el propio apunte dice explícitamente que excede los alcances del curso (usa herramientas de la teoría general de multigrid que no se desarrollan acá). Lo que **sí** hay que saber reproducir es la cadena de deducciones que va *desde* ese teorema (tomado como dato) *hasta* la cota de $\kappa_2$ del sistema precondicionado. Ver [[Multigrid en 1D - transferencia e implementacion]], [[Demostraciones - Capitulo 5]].

> [!pendiente] "Varios bloques aislados" en la cota de GC no tiene una receta garantizada
> A diferencia de los demás casos (un intervalo, un intervalo + un bloque aislado, pocos autovalores distintos), la construcción de un polinomio con una raíz por bloque (la del [[Parcial Modelo - Problema 3]]) **no siempre funciona** - depende de la posición relativa exacta de los bloques, de una forma difícil de predecir sin probar. Verificado con 3 corridas numéricas con la misma estructura y resultados muy distintos (un caso funciona bien, otro se estanca, otro explota). Si aparece este caso en el parcial: proponer la construcción, pero verificar graficando/calculando $\max|p|$ en cada bloque por separado antes de reportar una cota. Ver [[Elegir el polinomio para las cotas de GC]].

> [!pendiente] `IterativeSolvers.jl`: `jacobi`/`gauss_seidel`/`sor`/`ssor` no paran por tolerancia, y `ssor` parece tener un bug
> Estas cuatro funciones del paquete ejecutan exactamente `maxiter` iteraciones (converjan o no) sin avisar. Además, `ssor` diverge catastróficamente para $\omega\gtrsim1.3$ en un caso donde debería converger (Ostrowski-Reich) - verificado numéricamente. También: el `restart` default de `gmres` es $\min(20,n)$, o sea que por default es GMRes(20) con reinicio, no el GMRes completo del curso - sin `restart=n` explícito puede no converger en las iteraciones esperadas. Ver [[Referencia de Julia]].

## Capítulo 6

> [!pendiente] Ostrowski es solo suficiente y solo local; con $\rho(G'(x_*))=1$ no decide
> $G(x)=x-x^3$ y $G(x)=x+x^3$ tienen ambas $G'(0)=1$: la primera atrae (lentísimo, $x_k\sim1/\sqrt{2k}$) y la segunda no. Ver [[Iteracion de punto fijo y teorema de Ostrowski]], [[6.1 - Rho igual a 1 no decide|Problema 6.1]].

> [!pendiente] Para aplicar el Teorema 6.12 hay que verificar **las tres** hipótesis
> $D$ cerrado, **$G(D)\subset D$** (la que más se olvida) y $\|G'\|\leq\gamma<1$ en $D$ **convexo** (Corolario 6.19). "Acercar puntos" sin constante uniforme no alcanza: $\sqrt{1+x^2}$ no tiene punto fijo ([[6.3 - Contraccion en acotados|Problema 6.3]]). Receta en [[Criterio de contraccion via la derivada#Receta para los ejercicios]].

> [!pendiente] No hay Teorema del Valor Medio con igualdad en varias variables
> Solo la versión integral y la desigualdad $\|G(y)-G(x)\|\leq\max\|G'\|\,\|y-x\|$. Contraejemplo: $G=(x_1^3,x_1^2)$ ([[6.6 - No hay teorema del valor medio vectorial|Problema 6.6]]).

> [!pendiente] Newton cuadrático necesita **dos** hipótesis: $F'(x_*)$ invertible **y** $F'$ Lipschitz en $x_*$
> Si falla la primera (raíz múltiple), cae a **lineal** ($x^2$: tasa $\frac12$, [[6.14 - Newton para x al cuadrado|6.14]]). Si falla solo la segunda, queda **superlineal** de orden intermedio ($x+x^{1+\alpha}$: orden $1+\alpha$, [[6.15 - Newton para x mas x a la 1 mas alfa|6.15]]). Y puede ser **más** que cuadrático si $F''(x_*)=0$ ($Ax+h^2x^3$: orden 3, [[6.17 - Newton en problema discretizado|6.17]]). Ver [[Metodo de Newton - convergencia cuadratica]].

> [!pendiente] Medir el orden: en `Float64` quedan 2–3 cocientes útiles
> Usar $p\approx\log(e_{k+1}/e_k)/\log(e_k/e_{k-1})$ y repetir en `BigFloat` para verlo bien. Ver [[6.10 - Orden experimental de Newton|Problema 6.10]].

> [!pendiente] Newton-Baluev: sin convexidad no hay convergencia global
> $\arctan x$ diverge si $|x_0|>1.3917$; $2x+\sin x$ queda atrapado en el 2-ciclo $\pi\leftrightarrow-\pi$ (superatractor). Ambas tienen $F'>0$ y un único cero. Ver [[6.18 - Newton-Baluev sin convexidad|Problema 6.18]].

> [!pendiente] Newton sobre $\nabla f=0$ encuentra puntos críticos (no necesariamente mínimos) y puede ciclar
> Chequear la hessiana al final. En el [[6.13 - Minimo con Newton|Problema 6.13]] hay un ciclo atractor de período 5 para $x_0<0.2$.

> [!pendiente] Nelder–Mead (`fminsearch`) también depende del punto inicial
> Desde $(1.4,2.3)$ terminó en **otra** raíz del 6.11, por el símplex inicial grande de Optim.jl. Y la precisión alcanzable es $\sim\sqrt{\varepsilon_{\text{maq}}}$ ($h=\|F\|^2$ es cuadrática en el error). Ver [[6.11 - Cuatro ceros y fminsearch|Problema 6.11]].

## Capítulo 7

> [!pendiente] El autovector converge "en dirección": medir el error salvo signo/fase
> $\min(\|u_k-x_1\|,\|u_k+x_1\|)$; si $\lambda_1<0$, $u_k$ alterna de signo. Ver [[Metodo de las potencias]].

> [!pendiente] En hermíticas el autovalor converge el **doble** de rápido que el autovector
> El cociente de Rayleigh tiene error cuadrático: potencias con tasa $|\lambda_2/\lambda_1|^2$ en el autovalor; QR con la diagonal al cuadrado de la subdiagonal. Ver [[Observacion 7.4]], [[7.1 - Programar el metodo de las potencias|Problema 7.1]].

> [!pendiente] Potencias en $A_N$: $O(N^2)$ iteraciones; potencias inversas: $O(1)$
> $\lambda_{N-1}/\lambda_N\to1$ pero $\lambda_1/\lambda_2\to\frac14$. Ver [[7.2 - Potencias en la matriz A_N|7.2]], [[7.4 - Potencias inversas en la matriz A_N|7.4]].

> [!pendiente] Shift actualizado (cociente de Rayleigh): rapidísimo, pero no se elige el autovalor
> En $A_N$ desde $\mu_0=0$ convergió a $\lambda_1$, $\lambda_3$, $\lambda_4$ o $\lambda_2$ según $N$. Ver [[Metodo de las potencias inversas]], [[7.6 - Potencias inversas con y sin shift en A_N|Problema 7.6]].

> [!pendiente] El shift de Rayleigh se estanca en $A_N$ (no solo en $\begin{bmatrix}0&1\\1&0\end{bmatrix}$)
> Diagonal constante + espectro simétrico respecto de 2 ⟹ $(A_k)_{nn}=2$ para siempre y QR sin shift sobre pares $\pm$ de igual módulo. Usar **Wilkinson**. Ver [[7.9 - QR con shift y deflacion|Problema 7.9]], [[Algoritmo QR - Hessenberg, shifts y deflacion]].

> [!pendiente] El Problema 7.13 tal como está escrito es **falso**
> Hessenberg × Hessenberg tiene dos subdiagonales. Lo cierto (y lo que se usa en el 7.14): Hessenberg × triangular y triangular × Hessenberg son Hessenberg. Ver [[7.13 - Producto de matrices Hessenberg|Problema 7.13]].

> [!pendiente] "$A=QR$ con cualquier QR ⟹ $RQ$ Hessenberg" vale solo si $A$ es no singular
> Con $A$ singular la QR no es única y hay contraejemplo. Ver [[7.14 - QR preserva Hessenberg|Problema 7.14]].

> [!pendiente] Defectuoso ≠ repetido
> Autovalor doble semisimple: QR converge igual ([[7.10 - QR con autovalor semisimple doble|7.10]]). Defectuoso (Jordan): converge como $1/k$ ([[7.11 - QR con autovalor defectuoso|7.11]]).

> [!pendiente] Parar QR por "la diagonal no cambia" deja la subdiagonal en $\sim\sqrt{\text{tol}}$
> Parar por $|a_{i+1,i}|<\text{tol}$ (deflación). Ver [[7.7 - Implementar el algoritmo QR|Problema 7.7]].

> [!pendiente] El Teorema 7.11 (convergencia de QR) no se demuestra
> Sí hay que saber la explicación vía potencias / potencias inversas ([[Observacion 7.12]]) y qué pasa cuando fallan sus hipótesis.

## Ver también
- [[Inicio]]
- [[Tablas Comparativas - Indice]]
- [[Formulario - Indice]]
- [[Glosario - Indice|Glosario]]
