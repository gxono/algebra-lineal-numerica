---
tags: [parcial, capitulo4, capitulo5]
---

# Parcial Modelo - Problema 4 (35 pts): $-u''+u=f$ por diferencias finitas

## Enunciado

> Consideramos el problema de hallar la función $u$ que satisface
$$
\begin{cases}-u''(x)+u(x)=f(x), & x\in(0,1)\\u(0)=u(1)=0.\end{cases}
$$
> Mediante el método de diferencias finitas podemos obtener una aproximación de la solución $u$ en puntos de una partición del intervalo $[0,1]$. Más precisamente, si consideramos una partición de puntos equiespaciados $0=X_0<X_1<\dots<X_n<X_{n+1}=1$, $X_i=i/(n+1)$, el vector $x$ que aproxima $(u(X_1),\dots,u(X_n))^T$ es la solución del sistema lineal $Ax=b$ donde
> $$A = I_n + \frac1{h^2}\,\text{tridiag}(-1,2,-1), \qquad b=\begin{pmatrix}f(X_1)\\\vdots\\f(X_n)\end{pmatrix},$$
> con $I_n$ la identidad de $n\times n$ y $h=\frac1{n+1}$.
>
> **(a)** Demostrar que el método de Gauss-Seidel converge a la solución $x_*$ del sistema partiendo de cualquier valor inicial $x_0$. Más aún, hallar un valor $\mu_n\in(0,1)$ tal que $\|x_{k+1}-x_*\|_\infty\leq\mu_n\|x_k-x_*\|_\infty$, $k=0,1,\dots$
>
> **(b)** Considerar $n=1000$ y estimar la cantidad de iteraciones $k$ necesarias con el método de Gauss-Seidel para lograr que $\|x_k-x_*\|_\infty/\|x_0-x_*\|_\infty\leq10^{-5}$. Interpretar el número obtenido: ¿es un método viable para este problema? ¿Por qué?
>
> **(c)** Calcular el número de condición $\kappa_2(A)$, sabiendo que, para cada $i=1,\dots,n$, $\lambda_i=1+\dfrac4{h^2}\text{sen}^2\!\big(\dfrac{i\pi}{2(n+1)}\big)$ es autovalor con autovector $v_i=(\text{sen}(i\pi X_1),\dots,\text{sen}(i\pi X_n))^T$.
>
> **(d)** Considerar $n=1000$ y estimar la cantidad de iteraciones $k$ necesarias con el método del descenso más pronunciado y con el gradiente conjugado para lograr $\|x_k-x_*\|_A/\|x_0-x_*\|_A\leq10^{-5}$.
>
> **(e)** Considerando $f(x)=e^{-2\pi x^2}$, $n=1000$ y $x_0=0$, calcular las iteraciones del método del gradiente conjugado y graficar el decrecimiento de $\|x_k-x_*\|_A/\|x_0-x_*\|_A$. Verificar si se cumple la estimación que realizó en el ítem anterior. Explicar y justificar lo que se observa.
>
> **(f)** Considerando $f(x)=e^{-2\pi x^2}$, $n=1000$ y $x_0=0$, graficar el cociente entre el error y el residuo (ambos en norma dos) $\|x_k-x_*\|_2/\|r_k\|_2$, junto con las rectas horizontales $1/\lambda_{\min}$ y $1/\lambda_{\max}$. ¿Se mantiene el cociente entre ambas? ¿Se acerca a alguna de ellas? ¿Qué relación tiene la distancia entre las dos rectas con $\kappa_2(A)$? Explicar además qué ocurre con el gráfico para $k$ mayor que $n$ y por qué.
>
> **(g)** Considerar $f(x)=\text{sen}(\pi x)$ y un valor de $n\in\mathbb{N}$ arbitrario. Demostrar analíticamente que el método del gradiente conjugado partiendo de $x_0=0$ converge siempre en una iteración.

Es la misma familia que el [[Multigrid en 1D - problema modelo y suavizado|problema modelo de multigrid]] ($-u''=f$), pero con un **término de orden cero extra** ($+u$): eso suma $I_n$ a la matriz, lo cual **mejora el condicionamiento** respecto del caso puro (ver (c)).

## Resolución

### (a) Gauss-Seidel converge, con tasa explícita

**Paso 1 - $A$ es e.d.d.** La fila $i$ tiene $|a_{ii}|=1+2/h^2$. Las filas interiores tienen exactamente dos vecinos no nulos, cada uno de valor absoluto $1/h^2$, así que suma de off-diagonales $=2/h^2<1+2/h^2=|a_{ii}|$; las filas de borde ($i=1,n$) tienen *un solo* vecino, suma $1/h^2<2/h^2$, dominancia todavía más holgada. Luego $A$ es [[Estrictamente diagonalmente dominante|e.d.d.]], con razón de dominancia (por fila, $\sum_{j\neq i}|a_{ij}|/|a_{ii}|$) **uniformemente acotada** por
$$
\alpha_n := \frac{2/h^2}{1+2/h^2} = \frac{2}{h^2+2} \in(0,1)
$$
(la razón exacta en las filas interiores; estrictamente menor en las de borde). Por el [[Teorema de convergencia de Gauss-Seidel con e.d.d.|teorema general]] (e.d.d. $\implies\rho(M_{GS})<1$), Gauss-Seidel converge para todo $x_0$ - esto ya responde la primera parte del enunciado, pero no da ninguna *tasa*.

**Paso 2 - de dónde sale $\mu_n$: no es una fórmula que se propone "porque sí".** Para Jacobi, $\|M_J\|_\infty\le\alpha_n$ es inmediato (la fila $i$ de $M_J=-D^{-1}(L+U)$ tiene suma de valores absolutos exactamente $\sum_{j\neq i}|a_{ij}|/|a_{ii}|\le\alpha_n$, por definición de e.d.d.). La pregunta natural es: ¿vale la **misma** cota para Gauss-Seidel? $M_{GS}=-(L+D)^{-1}U$ involucra una inversa, así que sus entradas no se leen directamente de $A$ como en Jacobi - pero resulta que **sí**, con un argumento inductivo (fila por fila, aprovechando que GS usa valores *ya actualizados*): $\|M_{GS}\|_\infty\le\alpha_n$ también. Esta es la razón por la que $\mu_n:=\alpha_n$ es el candidato natural: es la misma cantidad que ya se usa para Jacobi (la razón de dominancia de la propia definición de e.d.d.), y resulta que **también** acota a Gauss-Seidel. Demostración completa, general (vale para cualquier matriz e.d.d. con razón uniforme $\alpha$, no solo esta tridiagonal - verificado además numéricamente con matrices aleatorias no tridiagonales): [[Cota de norma infinito de Gauss-Seidel para matrices e.d.d.]].

Aplicando ese teorema a esta $A$ con $\alpha=\alpha_n$:
$$
\|M_{GS}\|_\infty \leq \mu_n := \alpha_n = \frac{2}{h^2+2} = \frac{2(n+1)^2}{2(n+1)^2+1},
$$
y como $x_{k+1}-x_*=M_{GS}(x_k-x_*)$ (restando el punto fijo $x_*=M_{GS}x_*+c$ de la iteración, igual que en el [[Teorema de convergencia de metodos iterativos estacionarios|teorema general de convergencia]]), esto da directamente $\|x_{k+1}-x_*\|_\infty\le\mu_n\|x_k-x_*\|_\infty$ - la cota pedida.

Verificado numéricamente, comparando $\mu_n$ contra $\|M_{GS}\|_\infty$ calculado directamente sobre la matriz (construyendo $M_{GS}=-(L+D)^{-1}U$ y usando `opnorm(Mgs, Inf)`):

| $n$ | $\mu_n$ (cota) | $\|M_{GS}\|_\infty$ (real) |
|---|---|---|
| $5$ | $0.986$ | $0.915$ |
| $20$ | $0.999$ | $0.998$ |
| $100$ | $0.99995$ | $0.9999$ |
| $1000$ (el de (b)) | $0.999999501$ | $0.999999002$ |

La cota es válida en los cuatro casos ($\mu_n\ge\|M_{GS}\|_\infty$ siempre) y razonablemente ajustada - no exacta (la demostración "tira" un poco de margen al acotar $\alpha_n\sum_{j<i}|a_{ij}|$ por $\sum_{j<i}|a_{ij}|$ en un paso intermedio), pero del mismo orden.

### (b) $n=1000$: ¿es viable?

Usando el $\mu_{1000}$ calculado en (a):
$$
\mu_{1000} = \frac{2}{2+h^2} \approx 0.999999501, \qquad \mu_{1000}^k\leq10^{-5} \iff k \ge \frac{\ln(10^{-5})}{\ln\mu_{1000}} \implies k\approx2.31\times10^7.
$$
**No es viable.** $23$ millones de iteraciones de Gauss-Seidel (cada una $O(n)$) son completamente impracticables. El problema de fondo: $\mu_n\to1$ cuando $n\to\infty$ (porque $h\to0$ y $\mu_n=\frac{2}{h^2+2}\to1$), igual que en el [[Multigrid en 1D - problema modelo y suavizado|problema modelo puro]] - los métodos estacionarios clásicos degradan con la resolución de la malla, por eso el Capítulo 5 existe. (Usando la norma real $\|M_{GS}\|_\infty=0.999999002$ en vez de la cota $\mu_{1000}$ el número baja a $\approx1.15\times10^7$ - la mitad, pero sigue siendo igual de inviable en la práctica; la conclusión cualitativa no cambia.)

### (c) Número de condición vía los autovalores dados

$\lambda_i=1+\dfrac4{h^2}\text{sen}^2\!\Big(\dfrac{i\pi}{2(n+1)}\Big)$ es creciente en $i\in\{1,\dots,n\}$ ($\text{sen}^2$ creciente en $(0,\pi/2)$), así que $\lambda_{\min}=\lambda_1$, $\lambda_{\max}=\lambda_n$:
```julia
lambda(i,n) = 1 + 4/(1/(n+1))^2 * sin(i*pi/(2*(n+1)))^2
n = 1000
κ = lambda(n,n) / lambda(1,n)   # -> 368734.4976
```
$\kappa_2(A)=368734.4976$ (verificado contra `cond(A)` calculado directamente sobre la matriz: coincide). Comparar con el problema modelo puro $-u''=f$ (sin el $+u$), donde $\kappa_2=O(n^2)=O(4\times10^6)$ para $n=1000$ (ver [[Multigrid en 1D - problema modelo y suavizado]]) - acá sale **~10 veces mejor** porque el término $+u$ suma $1$ a todos los autovalores, lo cual ayuda proporcionalmente mucho más al autovalor chico ($\lambda_{\min}$ pasa de $O(h^2)\approx10^{-5}$ a $O(1)$) que al grande.

### (d) Iteraciones estimadas para SD y GC

**SD.** El [[Teorema 5.9|teorema de convergencia del descenso más pronunciado]] da $\|x_{k}-x_*\|_A\leq r^k\|x_0-x_*\|_A$ con $r=\dfrac{\kappa-1}{\kappa+1}$ - la misma estructura que (b) con Gauss-Seidel, solo que con otro factor. Resolviendo $r^k\leq10^{-5}$ igual que en (b):
$$
r = \frac{368734.5-1}{368734.5+1} = 0.99999458, \qquad k\geq\frac{\ln(10^{-5})}{\ln r} \implies k\approx2.1226\times10^6 \ (k=2122607).
$$

**GC.** La [[Observacion 5.20|Observación 5.20]] ya resuelve de antemano la desigualdad análoga para la cota de Chebyshev del [[Teorema 5.19]] ($\|x_k-x_*\|_A\leq2\big(\frac{\sqrt\kappa-1}{\sqrt\kappa+1}\big)^k\|x_0-x_*\|_A\leq\varepsilon\|x_0-x_*\|_A$), dando directamente la fórmula cerrada
$$
k\geq\frac{\sqrt\kappa+1}2\ln\frac2\varepsilon = \frac{\sqrt{368734.5}+1}{2}\ln\frac2{10^{-5}} \approx 3712.08 \implies k=3713.
$$

GC es $k_{SD}/k_{GC}\approx571.8$ veces más rápido que SD acá, y **miles** de veces más rápido que Gauss-Seidel (b) - la misma jerarquía de siempre ([[Tablas Comparativas - GC vs SD]]), mucho más marcada cuanto peor es $\kappa$ (la ventaja de GC escala como $\sqrt\kappa$ en el exponente, así que cuanto más grande $\kappa$, más se nota).

### (e) Corrida real con $f(x)=e^{-2\pi x^2}$

```julia
using LinearAlgebra, CairoMakie

# Igual que gradconj de [[Implementacion del gradiente conjugado]], pero guardando
# ‖x_k-x*‖_A en cada paso (necesita x* de antemano, por eso no es la versión genérica).
function gradconj_errA(A, b, x0, xstar; tol=1e-12, max_iter=5000)
    x = copy(x0)
    r = b - A*x
    ρ = r'r
    d = copy(r)
    ρold = 0.0
    iter = 0
    errA = [sqrt((x - xstar)'A*(x - xstar))]
    while sqrt(ρ) > tol && iter < max_iter
        if iter > 0
            β = ρ/ρold
            d = r + β*d
        end
        w = A*d
        α = ρ/(d'w)
        x += α*d
        r -= α*w
        ρold = ρ
        ρ = r'r
        iter += 1
        push!(errA, sqrt((x - xstar)'A*(x - xstar)))
    end
    return x, iter, errA
end

n = 1000
h = 1/(n+1)
A = SymTridiagonal(fill(1+2/h^2, n), fill(-1/h^2, n-1))
Xs = [i*h for i in 1:n]
b = exp.(-2π .* Xs.^2)
xstar = A \ b
x0 = zeros(n)

x, k, errA = gradconj_errA(A, b, x0, xstar)
errA_rel = errA ./ errA[1]    # ‖x_k-x*‖_A / ‖x_0-x*‖_A

fig = Figure()
ax = Axis(fig[1,1], yscale=log10, xlabel="iteración k", ylabel="‖xₖ-x*‖_A / ‖x₀-x*‖_A",
          title="GC: decrecimiento del error (f(x)=exp(-2πx²), n=1000)")
lines!(ax, 0:k, errA_rel)
save("parcial-modelo-p4e-decrecimiento-gc.png", fig)
```
![[parcial-modelo-p4e-decrecimiento-gc.png]]

Cruza $10^{-5}$ en $k=992$ (verificado: la corrida completa termina en $k=1000=n$ iteraciones hasta $\|r_k\|_2<10^{-12}$) - **menos** que la estimación de (d) ($3713$), pero la forma de la curva es reveladora: se mantiene casi plana hasta $k\approx900$ y recién ahí cae en picada. Esto **no** es el decaimiento geométrico suave que predice la cota genérica (factor fijo por paso) - es la [[GC como metodo directo|terminación finita]] dominando: con $\kappa$ tan grande, el factor geométrico por paso es $\approx1-2/\sqrt\kappa\approx0.9967$, casi invisible en pocos cientos de pasos, así que lo que realmente tira el error a cero es que GC resuelve **exacto** en a lo sumo $n=1000$ pasos ([[Teorema 5.38|terminación finita]], que para GC es [[Teorema 5.21]]/[[Teorema 5.22]]). La cota de (d) describe el peor caso teórico; acá gana la garantía de terminación finita antes de que el decaimiento geométrico "normal" llegue a notarse.

### (f) Cociente $\|e_k\|_2/\|r_k\|_2$ vs. $1/\lambda_{\min}$, $1/\lambda_{\max}$

```julia
function gradconj_ratio(A, b, x0, xstar; tol=1e-12, max_iter=1200)
    x = copy(x0)
    r = b - A*x
    ρ = r'r
    d = copy(r)
    ρold = 0.0
    iter = 0
    ratios = Float64[]
    while iter < max_iter
        ek = norm(x - xstar)
        rk = norm(r)
        push!(ratios, ek/rk)
        sqrt(ρ) <= tol && break
        if iter > 0
            β = ρ/ρold
            d = r + β*d
        end
        w = A*d
        α = ρ/(d'w)
        x += α*d
        r -= α*w
        ρold = ρ
        ρ = r'r
        iter += 1
    end
    return ratios
end

ratios = gradconj_ratio(A, b, x0, xstar)    # A, b, x0, xstar de (e)
λmin, λmax = eigmin(A), eigmax(A)           # coincide con la fórmula cerrada de (c): lambda(1,n), lambda(n,n)

fig2 = Figure()
ax2 = Axis(fig2[1,1], yscale=log10, xlabel="iteración k", ylabel="‖eₖ‖₂ / ‖rₖ‖₂",
           title="GC: cociente error/residuo vs 1/λmin, 1/λmax (n=1000)")
lines!(ax2, 0:length(ratios)-1, ratios, label="‖eₖ‖₂/‖rₖ‖₂")
hlines!(ax2, [1/λmin], color=:red, linestyle=:dash, label="1/λmin")
hlines!(ax2, [1/λmax], color=:green, linestyle=:dash, label="1/λmax")
axislegend(ax2, position=:rb)
save("parcial-modelo-p4f-cociente-error-residuo.png", fig2)
```
![[parcial-modelo-p4f-cociente-error-residuo.png]]

Verificado: $\lambda_{\min}\approx10.87$, $\lambda_{\max}\approx4.008\times10^6$ (coincide con `eigmin`/`eigmax` calculados sobre $A$ y con la fórmula cerrada $\lambda(1,n)$, $\lambda(n,n)$ de (c)), así que $1/\lambda_{\max}\approx2.50\times10^{-7}$ y $1/\lambda_{\min}\approx0.0920$ - las dos rectas del gráfico.

Como $r_k=b-Ax_k=A(x_*-x_k)=-Ae_k$, y $A$ es sdp: $\lambda_{\min}\|e_k\|_2\le\|Ae_k\|_2\le\lambda_{\max}\|e_k\|_2$.

> [!demostracion]- De dónde sale $\lambda_{\min}\|e\|_2\le\|Ae\|_2\le\lambda_{\max}\|e\|_2$
> Como $A$ es sdp, es diagonalizable ortogonalmente: $A=Q\Lambda Q^T$, con $Q$ ortogonal y $\Lambda=\text{diag}(\lambda_1,\dots,\lambda_n)$, $0<\lambda_{\min}\leq\lambda_i\leq\lambda_{\max}$. Sea $e\in\mathbb R^n$ cualquiera y $c:=Q^Te$ (sus coordenadas en la base de autovectores; $\|c\|_2=\|e\|_2$ porque $Q$ es ortogonal, ver [[Matrices ortogonales]]). Entonces
> $$\|Ae\|_2^2 = e^TA^TAe = e^TA^2e = (Q^Te)^T\Lambda^2(Q^Te) = c^T\Lambda^2c = \sum_i\lambda_i^2c_i^2.$$
> Como $\lambda_{\min}^2\leq\lambda_i^2\leq\lambda_{\max}^2$ para todo $i$ (los $\lambda_i>0$ al ser $A$ sdp, así que elevar al cuadrado preserva el orden):
> $$\lambda_{\min}^2\sum_ic_i^2 \leq \sum_i\lambda_i^2c_i^2 \leq \lambda_{\max}^2\sum_ic_i^2 \quad\Longrightarrow\quad \lambda_{\min}^2\|e\|_2^2 \leq \|Ae\|_2^2 \leq \lambda_{\max}^2\|e\|_2^2,$$
> usando $\sum_ic_i^2=\|c\|_2^2=\|e\|_2^2$. Tomando raíz cuadrada (todo no negativo) se obtiene la desigualdad. *(Es, en esencia, la misma idea que [[4.5(d) - Formula norma 2 matriz simetrica|Problema 4.5(d)]] - $\|A\|_2=\lambda_{\max}$ para $A$ simétrica - aplicada también a la cota inferior, vía $A^{-1}$: $\|A^{-1}\|_2=1/\lambda_{\min}$ por el mismo argumento, ver [[Lema 4.14]].)*

Aplicando esto a $e=e_k$, y usando $\|r_k\|_2=\|Ae_k\|_2$:
$$
\lambda_{\min}\|e_k\|_2\leq\|r_k\|_2\leq\lambda_{\max}\|e_k\|_2 \quad\Longrightarrow\quad \frac1{\lambda_{\max}} \le \frac{\|e_k\|_2}{\|r_k\|_2} \le \frac1{\lambda_{\min}}.
$$
El cociente observado ($\approx10^{-3}$ a $10^{-2.5}$ para $k<n$) se mantiene **cómodamente dentro** de la banda ($1/\lambda_{\max}\approx2.5\times10^{-7}$, $1/\lambda_{\min}\approx0.092$), sin pegarse a ninguno de los dos extremos: el residuo en cualquier iteración intermedia tiene componentes repartidas entre varios modos, no alineado con el autovector extremo. La **distancia** entre las dos rectas (en escala log) es $\log(\lambda_{\max}/\lambda_{\min})=\log\kappa_2(A)$ - cuanto peor condicionada $A$, más ancha la banda en la que puede moverse el cociente, y más lejos puede estar $\|e_k\|$ de ser proporcional a $\|r_k\|$.

**Para $k\ge n$:** el cociente se vuelve errático (pico y caída brusca). Por la terminación finita de (e), para $k$ cercano a $n=1000$ tanto $\|e_k\|$ como $\|r_k\|$ ya están en el piso de precisión de máquina simultáneamente - el cociente de dos cantidades que son básicamente ruido de redondeo deja de tener el significado teórico de arriba (mismo fenómeno de "piso de punto flotante" que apareció en [[5.12 - Multigrid como metodo iterativo|Problema 5.12]](d) con $\rho$ medido).

### (g) $f(x)=\text{sen}(\pi x)$: GC converge en exactamente 1 paso

Con $f(x)=\text{sen}(\pi x)$, $b_i=\text{sen}(\pi X_i)$ para $i=1,\dots,n$ - que es **exactamente** el primer autovector dado en el enunciado, $v_1=(\text{sen}(\pi X_1),\dots,\text{sen}(\pi X_n))^T$, evaluado con $i=1$. Es decir $b=v_1$, autovector puro de $A$ con autovalor $\lambda_1$.

> [!demostracion]-
> Con $x_0=0$: $r_0=b-Ax_0=b=v_1$, combinación lineal de **un solo** autovector. Por el [[Teorema 5.22]] (si $b$ es combinación lineal de $k$ autovectores y $x_0=0$, GC converge en a lo sumo $k$ pasos), con $k=1$: GC converge en a lo sumo $1$ paso. Y no puede ser $0$ pasos salvo que $b=0$, así que converge en **exactamente** $1$ paso (si $x_0$ no es ya la solución).

Verificado: corriendo GC con este $b$, el error $\|x_1-x_*\|_A/\|x_0-x_*\|_A$ ya está en $3.5\times10^{-13}$ (precisión de máquina) después de la primera iteración, para **cualquier** $n$ (se verificó para $n=1000$; el argumento de arriba no depende de $n$, solo de que $v_1$ sea autovector, lo cual vale para todo $n\in\mathbb{N}$ por construcción del enunciado).

## Ver también
- [[Parcial Modelo - Indice|Parcial Modelo]]
- [[Metodo de Gauss-Seidel]]
- [[GC como metodo directo]]
- [[Multigrid en 1D - problema modelo y suavizado]]
- [[Cotas de error del gradiente conjugado]]
