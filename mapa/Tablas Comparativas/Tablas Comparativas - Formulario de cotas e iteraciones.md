---
tags: [moc, comparacion, capitulo4, capitulo5]
---

# Tablas de consulta rápida: $\kappa_2$, iteraciones, error y residuo por método

Página pensada para las consignas que más se repiten en los parciales ("calcular $\kappa_2(A)$", "estimar cuántas iteraciones hacen falta", "¿cuál es el error/residuo después de $k$ iteraciones?"): en vez de ir a buscar cada fórmula a su nota de teoría, están todas acá, con el link a la demostración si hace falta reconstruir de dónde salen.

## Glosario: qué es cada símbolo y cómo se calcula si no está dado

Estos símbolos se repiten en casi todas las tablas de abajo - se definen una sola vez acá para no repetir en cada fila. Los símbolos específicos de una sola tabla ($q(z)$, $c_1,c_2$) se explican junto a esa tabla.

| Símbolo | Qué es | Cómo se calcula si no está dado |
|---|---|---|
| $n$ | tamaño del sistema ($A\in\mathbb R^{n\times n}$) | dato del problema |
| $h$ | paso de malla (en problemas de diferencias finitas) | $h=\frac1{n+1}$ |
| $\lambda_i$ | autovalores de $A$ | dados en el enunciado, o `eigvals(A)` |
| $\lambda_{\min},\lambda_{\max}$ | mínimo y máximo autovalor de $A$ (en valor absoluto si $A$ no es sdp) | si no están dados: `eigmin(A)`, `eigmax(A)` (sirven directo si $A$ es sdp); para $A$ simétrica no necesariamente sdp, tomar mínimo/máximo de $\lvert\lambda_i\rvert$ sobre `eigvals(A)` |
| $\sigma_{\min},\sigma_{\max}$ | mínimo y máximo valor singular de $A$ (solo hace falta si $A$ **no** es simétrica) | `svdvals(A)`, o $\sqrt{\lambda_{\min}(A^TA)}$ / $\sqrt{\lambda_{\max}(A^TA)}$ |
| $\kappa_2(A)$ | número de condición en norma 2 | ver Tabla 1 completa más abajo - depende de qué se conoce de $A$ |
| $\alpha$ (razón de dominancia) | qué tan "e.d.d." es $A$, fila por fila | $\alpha=\max_i\dfrac{\sum_{j\neq i}\lvert a_{ij}\rvert}{\lvert a_{ii}\rvert}$ - se lee directo de las filas de $A$, sin autovalores ni factorizaciones |
| $M$ (matriz de iteración) | la matriz tal que $e_{k+1}=Me_k$ para el método estacionario en cuestión | $M_J=-D^{-1}(L+U)$ (Jacobi), $M_{GS}=-(L+D)^{-1}U$ (Gauss-Seidel), $M_{SOR}=(D/\omega+L)^{-1}\big((1-\omega)/\omega\,D-U\big)$ (SOR) - $D,L,U$ son la diagonal, y las partes estrictamente triangular inferior/superior de $A$ |
| $\rho(M)$ | radio espectral de $M$ (el autovalor de módulo más grande) | si no hay fórmula cerrada: construir $M$ explícitamente y `maximum(abs.(eigvals(M)))` |
| $\varepsilon$ | tolerancia relativa deseada ($\|\text{error}_k\|/\|\text{error}_0\|\leq\varepsilon$, o lo mismo con el residuo) | dato del enunciado |
| $k$ | número de iteraciones | lo que casi siempre se pide despejar (Tabla 3) o se da para evaluar la cota (Tabla 4) |
| $e_k:=x_k-x_*$ | error en la iteración $k$ ($x_*$ = solución exacta) | no se puede calcular en la práctica sin conocer $x_*$ (si se conoce $x_*$, como en un ejercicio de verificación, es una resta directa) |
| $r_k:=b-Ax_k$ | residuo en la iteración $k$ | siempre calculable (solo necesita $x_k$, no $x_*$) - por eso se usa como criterio de parada práctico |
| $\text{tol}$ | tolerancia sobre el residuo **relativo** usada como criterio de corte de un algoritmo ($\|r_k\|_2/\|b\|_2<\text{tol}$) | dato del enunciado o del código (parámetro `tol` de la función); no confundir con $\varepsilon$, que es la tolerancia sobre el **error** relativo - están relacionadas por la fila "criterio de parada práctico" de la Tabla 5 |

## 1. Cómo calcular $\kappa_2(A)$, según lo que se conoce de $A$

| Se conoce... | $\kappa_2(A)$ | Por qué |
|---|---|---|
| $A$ simétrica (autovalores con signo) | $\dfrac{\max_i\lvert\lambda_i\rvert}{\min_i\lvert\lambda_i\rvert}$ | [[Lema 4.14]] ($\sigma_i=\lvert\lambda_i\rvert$ cuando $A$ es simétrica) |
| $A$ sdp (autovalores $>0$) | $\dfrac{\lambda_{\max}}{\lambda_{\min}}$ | caso particular de lo anterior, sin valor absoluto |
| $A$ general (no simétrica) | $\dfrac{\sigma_{\max}(A)}{\sigma_{\min}(A)}$ (valores singulares, no autovalores) | definición; ver [[Numero de condicion y error relativo]] |
| Se conoce $\kappa_2(A^TA)$ en vez de $\kappa_2(A)$ (p. ej. ecuaciones normales) | $\kappa_2(A) = \sqrt{\kappa_2(A^TA)}$, o al revés $\kappa_2(A^TA)=\kappa_2(A)^2$ | [[Lema 4.17]] - **clave** en [[Parcial Modelo - Problema 2]]: las ecuaciones normales elevan $\kappa$ al cuadrado |
| $A$ triangular | autovalores = diagonal, $\kappa_2$ con esos (si además es simétrica/normal; si no, $\kappa_2$ real puede diferir mucho de la razón de autovalores) | [[Autovalores de una matriz triangular]] |
| Espectro dado como unión de intervalos/bloques, $A$ sdp | $\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$ usando directamente los extremos globales dados (no hace falta ir bloque por bloque) | [[Parcial Modelo - Problema 3]](a) |
| Matriz tridiagonal Toeplitz $\alpha I+\beta\,\text{tridiag}(-1,2,-1)$ | autovalores cerrados $\lambda_i=\alpha+\frac{4\beta}{h^2}\text{sen}^2\!\big(\frac{i\pi}{2(n+1)}\big)$, creciente en $i$ $\Rightarrow\lambda_{\min}=\lambda_1,\ \lambda_{\max}=\lambda_n$ | [[Parcial Modelo - Problema 4]](c), [[Multigrid en 1D - problema modelo y suavizado]] |

## 2. Tasa de convergencia $r$ por iteración, según el método

Toda la familia de cotas de abajo tiene la misma forma $\|\text{error}_k\|\leq C\cdot r^k\|\text{error}_0\|$ (con $C=1$ salvo que se diga lo contrario) - lo único que cambia es qué es $r$ y en qué norma.

| Método | Hipótesis | $r$ (factor por iteración) | Norma de la cota | Demostración |
|---|---|---|---|---|
| [[Metodo de Jacobi\|Jacobi]] | $A$ e.d.d., razón $\alpha$ | $r=\alpha=\max_i\sum_{j\neq i}\lvert a_{ij}\rvert/\lvert a_{ii}\rvert$ | $\|\cdot\|_\infty$ | $\|M_J\|_\infty\leq\alpha$, inmediato de la definición |
| [[Metodo de Gauss-Seidel\|Gauss-Seidel]] | $A$ e.d.d., razón $\alpha$ uniforme | $r=\alpha$ (**misma cota que Jacobi**) | $\|\cdot\|_\infty$ | [[Cota de norma infinito de Gauss-Seidel para matrices e.d.d.]] |
| [[Metodo SOR\|SOR]]($\omega$) | $\rho(M_{SOR})<1$ (hay que calcularlo para el $\omega$ elegido, no hay fórmula cerrada simple salvo casos particulares) | $r=\rho(M_{SOR})$ | asintótica (radio espectral) | [[Metodo SOR]] |
| [[Metodo de Richardson\|Richardson]]/SD, $\alpha$ óptimo | $A$ sdp | $r=\dfrac{\kappa_2(A)-1}{\kappa_2(A)+1}$ | $\|\cdot\|_A$ | [[Teorema 5.9]] |
| GC, cota genérica (solo con $\kappa_2$) | $A$ sdp | $r=\dfrac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}$ (con $C=2$ adelante) | $\|\cdot\|_A$ | [[Teorema 5.19]] / [[Observacion 5.20]] |
| GC, con distribución de autovalores (bloques aislados) | $A$ sdp, $\sigma(A)$ en bloques conocidos | $r=\max_{z\in\sigma(A)}\lvert q(z)\rvert$ para el $q$ elegido (polinomio auxiliar) | $\|\cdot\|_A$ | [[Elegir el polinomio para las cotas de GC]], [[Corolario 5.14]] |
| GCP | $M\approx A^{-1}$, $c_1\leq x^TMx/x^TA^{-1}x\leq c_2$ | $r=\dfrac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}$ (reemplaza $\kappa_2(A)$ por $c_2/c_1$) | $\|\cdot\|_A$ | [[5.10 - Cota de kappa2 SAS\|Problema 5.10]] |
| Multigrid, ciclo en V | hipótesis del Teorema 5.34 (suavizado + aproximación) | $r=\rho$ (medido, no hay fórmula cerrada - ver [[5.12 - Multigrid como metodo iterativo\|Problema 5.12]](d)) | $\|r_k\|_2$ (residuo) | [[Multigrid en 1D - transferencia e implementacion]] |
| GMRes | en general, **no hay** fórmula cerrada en términos de $\kappa_2(A)$ solo (depende de todo $\sigma(A)$ vía $\kappa_2(V)$); si $A$ es normal, sí se puede acotar por $\max_{z\in\sigma(A)}\lvert p(z)\rvert$ igual que GC | - | $\|r_k\|_2$ (residuo) | [[GMRes - definicion y propiedad de minimizacion|Teorema 5.39]] |

**Variables específicas de esta tabla** (las demás - $\alpha,\kappa_2(A),\rho$ - están en el glosario de arriba):
- $C$: constante que multiplica a $r^k$ en la cota ($C=1$ salvo donde se indica $C=2$, que aparece solo en la cota genérica de GC - [[Teorema 5.19]]).
- $q(z)$ (GC con bloques aislados): el polinomio auxiliar **elegido** (no dado, hay que construirlo) - ver la receta completa en [[Elegir el polinomio para las cotas de GC]] (un factor con una raíz en cada bloque, elevado a una potencia). $\max_{z\in\sigma(A)}\lvert q(z)\rvert$ se calcula evaluando $q$ en los extremos de cada bloque (ver [[Parcial Modelo - Problema 3]](d) para el argumento de por qué los extremos alcanzan).
- $c_1,c_2$ (GCP): cotas del cociente $x^TMx/x^TA^{-1}x$ sobre todo $x\neq0$. Si no están dadas directamente, se obtienen como $c_1=\lambda_{\min}(MA)$, $c_2=\lambda_{\max}(MA)$ (autovalores del producto $MA$, que son reales porque $MA$ es semejante a la matriz simétrica $A^{1/2}MA^{1/2}$ - ver [[5.10 - Cota de kappa2 SAS|Problema 5.10]] y el cálculo numérico concreto en [[5.13 - Multigrid como precondicionador|Problema 5.13]](b)).

## 3. Iteraciones necesarias para $\|\text{error}_k\|/\|\text{error}_0\|\leq\varepsilon$

Despejando $k$ de $C\cdot r^k\leq\varepsilon$ (y usando $\ln r<0$, por eso se da vuelta la desigualdad al dividir):
$$
k \geq \frac{\ln(\varepsilon/C)}{\ln r}.
$$

| Método | $k\geq\dots$ | Ejemplo verificado |
|---|---|---|
| Jacobi/Gauss-Seidel (e.d.d., razón $\alpha$) | $\dfrac{\ln\varepsilon}{\ln\alpha}$ | [[Parcial Modelo - Problema 4]](b): $\alpha=\mu_{1000}\approx0.9999995$, $k\approx2.31\times10^7$ |
| Richardson/SD óptimo | $\dfrac{\ln\varepsilon}{\ln r}$, $r=\frac{\kappa-1}{\kappa+1}$ | Problema 4(d): $\kappa=368734.5\Rightarrow k=2122607$ |
| GC genérico | $\dfrac{\ln(\varepsilon/2)}{\ln r}$, $r=\frac{\sqrt\kappa-1}{\sqrt\kappa+1}$ - forma cerrada ya resuelta: $k\geq\frac{\sqrt\kappa+1}2\ln\frac2\varepsilon$ | Problema 4(d): $\kappa=368734.5\Rightarrow k=3713$ |
| GC con bloques aislados | $\dfrac{\ln\varepsilon}{\ln(\max\lvert q\rvert)}$ (igual forma, con el $r$ de la fila de arriba) | [[Parcial Modelo - Problema 3]](d): $k=17$ repeticiones ($3k=51$) |
| GCP | $\dfrac{\sqrt{c_2/c_1}+1}2\ln\frac2\varepsilon$ (mismo formato que GC, con $c_2/c_1$ en vez de $\kappa$) | [[5.13 - Multigrid como precondicionador\|Problema 5.13]](d): $k$ predicho $9$ vs. observado $8$ ($\nu=1$) |
| Multigrid / GMRes precondicionado | $\dfrac{\ln\varepsilon}{\ln\rho}$ | [[5.12 - Multigrid como metodo iterativo\|Problema 5.12]](c): no crece con $n$ |

**Receta para resolver cualquier consigna de este tipo:** identificar el método → fila de la Tabla 2 para $r$ (o $C\cdot r$) → aplicar la fórmula de la Tabla 3 → listo. Si la consigna da la distribución de autovalores en vez de solo $\kappa_2$, usar la fila "con distribución de autovalores" en vez de la genérica (da un $k$ mucho más chico - ver la comparación completa en [[Parcial Modelo - Problema 3]]).

## 4. Error (o residuo) dado $k$ fijo - la cuenta inversa

Es la fila de la Tabla 2 **sin despejar $k$**: evaluar directamente $C\cdot r^k$ (o $\rho^k$) con el $k$ dado.
$$
\frac{\|x_k-x_*\|}{\|x_0-x_*\|} \leq C\cdot r^k.
$$
Ejemplo: [[Parcial Modelo - Problema 4]](e) - la cota teórica con $k=3713$ (de la Tabla 3) da $\varepsilon=10^{-5}$ exactamente por construcción, pero la corrida **real** cruza $10^{-5}$ en $k=992$ (mucho antes): la cota genérica es solo una cota superior, no hay que esperar que sea ajustada - en ese caso ganó la [[GC como metodo directo|terminación finita]] (Teorema 5.21/5.22), no el decaimiento geométrico "normal".

> [!idea]
> Cuando la cota teórica y el valor observado difieren mucho, **no es un error** - es información: casi siempre quiere decir que hay algo más fuerte en juego (terminación finita si $b$ tiene pocas componentes espectrales, o que el espectro tiene estructura que la cota genérica no ve). Ver también [[5.8 - Cota pesimista autovalor aislado|Problema 5.8]].

## 5. Relación error $\leftrightarrow$ residuo

$r_k=b-Ax_k=A(x_*-x_k)=-Ae_k$ siempre (por definición), así que acotar uno acota el otro - pero la relación exacta depende de si $A$ es sdp o no:

| Hipótesis sobre $A$ | Relación | Demostración |
|---|---|---|
| $A$ sdp | $\lambda_{\min}\|e_k\|_2\leq\|r_k\|_2\leq\lambda_{\max}\|e_k\|_2 \iff \dfrac1{\lambda_{\max}}\leq\dfrac{\|e_k\|_2}{\|r_k\|_2}\leq\dfrac1{\lambda_{\min}}$ | [[Parcial Modelo - Problema 4]](f), via diagonalización ortogonal |
| $A$ general (no simétrica) | $\sigma_{\min}(A)\|e_k\|_2\leq\|r_k\|_2\leq\sigma_{\max}(A)\|e_k\|_2$ (mismo argumento, con valores singulares en vez de autovalores) | análogo, vía $A=U\Sigma V^T$ |
| Criterio de parada práctico (cualquier $A$ no singular) | $\dfrac{\|r_k\|_2}{\|b\|_2}<\text{tol} \implies \dfrac{\|x_k-x_*\|_2}{\|x_*\|_2}\leq\kappa_2(A)\cdot\text{tol}$ | [[Implementacion del gradiente conjugado]] - "Criterio de parada"; útil para traducir una tolerancia en residuo a una cota de error relativo |
| $\|\cdot\|_A$ vs. $\|\cdot\|_2$ (para GC/SD) | $\sqrt{\lambda_{\min}}\|e_k\|_2\leq\|e_k\|_A\leq\sqrt{\lambda_{\max}}\|e_k\|_2$ | mismo argumento espectral que arriba, aplicado a $\|z\|_A^2=z^TAz$ |

**Por qué importa la distancia entre $1/\lambda_{\max}$ y $1/\lambda_{\min}$:** en escala logarítmica es $\log(\lambda_{\max}/\lambda_{\min})=\log\kappa_2(A)$ - cuanto peor condicionada $A$, más ancha la banda en la que puede moverse el cociente $\|e_k\|/\|r_k\|$, y menos informativo es el residuo solo como indicador del error real (ver [[Parcial Modelo - Problema 4]](f) para el experimento completo).

## 6. Qué hipótesis hace falta para que cada método converja

Distinguiendo **condición necesaria y suficiente exacta** (siempre $\rho(M)<1$ para los estacionarios, por [[Teorema de convergencia de metodos iterativos estacionarios|el teorema general]]) de **condiciones suficientes** más fáciles de chequear a simple vista (que no siempre son necesarias - de ahí los contraejemplos).

| Método | Condición suficiente (fácil de chequear) | Exacta (necesaria y suficiente) | Si la suficiente falla... |
|---|---|---|---|
| Jacobi | $A$ e.d.d. | $\rho(M_J)<1$ | puede converger igual (contraejemplo: [[4.11 - Convergencia Jacobi y Gauss-Seidel caso particular\|Problema 4.11]], $A$ no e.d.d. pero Jacobi converge) o no - hay que calcular $\rho(M_J)$ |
| Gauss-Seidel | $A$ e.d.d., **o** $A$ sdp | $\rho(M_{GS})<1$ | mismo contraejemplo 4.11: ahí $A$ no es e.d.d. (ni sdp) y Gauss-Seidel **no** converge - muestra que ninguna de las dos hipótesis es necesaria, y que Jacobi y Gauss-Seidel no se pueden comparar en general (ahí Jacobi converge y GS no, al revés de la intuición) |
| SOR($\omega$) | $A$ sdp **y** $\omega\in(0,2)$ (Ostrowski-Reich, [[4.13 - Convergencia global de SOR para A sdp\|Problema 4.13]]) | $\rho(M_{SOR})<1$ | $\omega\in(0,2)$ es además **necesaria para cualquier $A$** ([[Metodo SOR]]: si $\omega\notin(0,2)$, $\rho(M_{SOR})\geq\lvert1-\omega\rvert\geq1$, nunca converge); sin $A$ sdp, hay que calcular $\rho(M_{SOR})$ caso por caso |
| Richardson/SD | $A$ sdp (y ahí **cualquier** $\alpha\in(0,2/\lambda_{\max})$ converge) | existe algún $\alpha\in\mathbb R$ con $\rho(M_\alpha)<1$ | puede no existir **ningún** $\alpha$ que funcione (contraejemplo: [[4.16 - Matriz sin alpha de Richardson que converja\|Problema 4.16]], autovalores de signos opuestos) |
| GC / GCP | $A$ sdp (GCP además necesita $M$ sdp) | - (GC no está definido si $A$ no es sdp; para $A$ general hace falta CGNR/CGNE o GMRes) | - |
| CGNR / CGNE | $A$ no singular (cualquiera, no necesita ser simétrica) | $A^TA$ (o $AA^T$) sdp, lo cual vale automáticamente si $A$ tiene rango completo | si $A$ no tiene rango completo, $A^TA$ es solo sdp (no definida positiva), GC sobre $A^TA$ puede no converger a una única solución |
| GMRes | $A$ no singular (nada más - no exige simetría) | $A$ no singular es ya necesaria y suficiente (termina en $\leq n$ pasos siempre) | - |
| Multigrid (ciclo en V) | hipótesis técnicas del Teorema 5.34 (suavizado + aproximación) - no se demuestran en el curso, se toman como dato | - | sin esas hipótesis (p. ej. coeficientes muy oscilantes), el suavizador puede no amortiguar los modos oscilatorios y el ciclo deja de converger a tasa fija |

**Variables:** $M_J,M_{GS},M_{SOR},M_\alpha$ son las matrices de iteración del glosario (con $M_\alpha=I-\alpha A$ para Richardson); $\rho(\cdot)$ es el radio espectral (glosario). $\omega$ es el parámetro de sobre-relajación de SOR (dato que se elige, no se calcula). $\alpha$ en esta tabla es el **paso de Richardson** (un escalar que se elige) - **no confundir** con el $\alpha$ de la Tabla 2/glosario (razón de dominancia e.d.d.): mismo símbolo, dos usos distintos según el método, hay que mirar el contexto.

## 7. Costo por iteración y memoria, todos los métodos juntos

Asumiendo $A$ sparse con $\text{nnz}(A)$ entradas no nulas (para una matriz tridiagonal o pentadiagonal, $\text{nnz}(A)=O(n)$; para una matriz densa genérica, $\text{nnz}(A)=O(n^2)$ y estos costos "por iteración" se multiplican por $n$ en vez de quedar en $O(n)$).

| Método | Costo por iteración | Memoria | Costo acumulado hasta $k$ |
|---|---|---|---|
| Jacobi / Gauss-Seidel / SOR | $O(\text{nnz}(A))$ (un producto matriz-vector, o una resolución triangular equivalente) | $O(n)$ (un vector $x$, fijo) | $O(k\cdot\text{nnz}(A))$ |
| Richardson / SD | $O(\text{nnz}(A))$ (1 mat-vec + 2 productos escalares) | $O(n)$, fija (2-3 vectores) | $O(k\cdot\text{nnz}(A))$ |
| GC | $O(\text{nnz}(A))$ (1 mat-vec, 2 prod. escalares, 3 $axpy$ - ver [[Implementacion del gradiente conjugado]]) | $O(n)$, **fija** (4 vectores: $x,r,d,w$) | $O(k\cdot\text{nnz}(A))$ |
| GCP | $O(\text{nnz}(A))+O(\text{costo de aplicar }M)$ | $O(n)$, fija | $O\big(k\cdot(\text{nnz}(A)+\text{costo de }M)\big)$ |
| Multigrid (1 ciclo en V) | $O(n)$ (suma geométrica sobre todos los niveles: $n+n/2+n/4+\dots\approx2n$) | $O(n)$ (mismo argumento geométrico sobre todas las mallas) | $O(k\cdot n)$ - y $k=O(1)$, independiente de $n$, así que el total es $O(n)$ |
| GMRes (sin restart) | $O(\text{nnz}(A))+O(nk)$ en la iteración $k$ (la ortogonalización contra **toda** la base de Krylov crece con $k$) | $O(nk)$, **creciente** (hay que guardar los $k$ vectores de la base) | $O(nk^2)$ - domina el costo de ortogonalización, no el de los mat-vec |
| GMRes($m$) (con restart) | igual que GMRes, pero $k$ se resetea a $0$ cada $m$ iteraciones | $O(nm)$, acotada (nunca crece más allá de $m$ vectores) | depende de cuántos restarts hagan falta - ver la advertencia de [[5.16 - Programar GMRes con restart\|Problema 5.16]] (sin buen precondicionador, el restart puede necesitar muchos más mat-vec en total que sin restart) |
| CGNR / CGNE | $O(\text{nnz}(A))$, pero **2 mat-vec** por iteración ($A$ y $A^T$) en vez de 1 | $O(n)$, fija | $O(k\cdot\text{nnz}(A))$, con el factor $2\times$ ya contado |

**Variables:** $n$ = tamaño del sistema; $\text{nnz}(A)$ = cantidad de entradas no nulas de $A$ (si no se especifica estructura, asumir denso, $\text{nnz}(A)=n^2$); $k$ = número de iteraciones (de la Tabla 3 si hace falta estimarlo); $m$ = tamaño de reinicio elegido para GMRes($m$) (dato del enunciado, no se calcula).

## 8. Capítulos 6 y 7: tasas, órdenes e iteraciones

Los métodos lineales de los Capítulos 6–7 se tratan igual que los de arriba ($\|\text{error}_k\|\leq C r^k$, $k\geq\ln(\varepsilon/C)/\ln r$). Los superlineales se cuentan en **dígitos**: con orden $p$, los dígitos correctos se multiplican por $\approx p$ en cada paso.

| Método | Hipótesis | $r$ u orden $p$ | Fuente |
|---|---|---|---|
| Punto fijo, contracción | $G(D)\subset D$, $\Vert G'\Vert \leq\gamma<1$ en $D$ convexo cerrado | $r=\gamma$; a priori $k\geq\frac{\ln(\varepsilon(1-\gamma)/\Vert x_1-x_0\Vert )}{\ln\gamma}$ | [[Teorema 6.12]] |
| Punto fijo, local | $\rho(G'(x_*))<1$ | $r\to\rho(G'(x_*))$ | [[Teorema 6.8]] |
| Newton | $F'(x_*)$ invertible, $F'$ Lipschitz | $p=2$ ($c=3L\Vert F'(x_*)^{-1}\Vert $) | [[Teorema 6.24]] |
| Newton, raíz múltiple (1D) | multiplicidad $m$ | $r=1-\frac1m$ | [[6.14 - Newton para x al cuadrado\|6.14]] |
| Newton simplificado | $F'(x_0)$ cerca de $F'(x_*)$ | $r\approx\rho(I-F'(x_0)^{-1}F'(x_*))$ | [[6.16 - Newton simplificado\|6.16]] |
| Potencias | $\vert \lambda_1\vert >\vert \lambda_2\vert $ | $r=\vert \lambda_2/\lambda_1\vert $ (autovector); $r^2$ para el autovalor si $A$ es hermítica | [[Proposicion 7.2]] |
| Potencias inversas, shift $\mu$ | $\mu$ más cerca de $\lambda_*$ que de los demás | $r=\rho_\mu=\max_i\frac{\vert \lambda_*-\mu\vert }{\vert \lambda_i-\mu\vert }$ | [[Proposicion 7.7]] |
| Cociente de Rayleigh (shift actualizado) | cerca del autopar | $p=3$ hermítica, $p=2$ general | [[7.5 - Potencias inversas con shift actualizado\|7.5]] |
| QR sin shift | Teorema 7.11 | $r=\max_i\vert \lambda_{i+1}/\lambda_i\vert $ | [[Teorema 7.11]] |
| QR con shift (Rayleigh/Wilkinson) + deflación | genérico | $p=3$ hermítica (unas 2–4 iteraciones por autovalor) | [[7.9 - QR con shift y deflacion\|7.9]] |

Para $A_N=\text{tridiag}(-1,2,-1)$: potencias $r\approx1-\frac{3\pi^2}{4(N+1)^2}$ ($k=O(N^2)$); potencias inversas con $\mu=0$, $r\to\frac14$ ($k$ independiente de $N$). Ver [[Formulario - Capitulo 7]].

## Ver también
- [[Tablas Comparativas - Indice]]
- [[Formulario - Capitulo 4]] · [[Formulario - Capitulo 5]] (fórmulas de construcción de cada método, sin la parte de "cuántas iteraciones")
- [[Parcial Modelo - Indice]] (las consignas de donde salió la mayoría de esta tabla)
