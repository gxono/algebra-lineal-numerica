---
tags: [parcial, capitulo5]
---

# Parcial Modelo - Problema 3 (30 pts): estimaciones a partir del espectro

## Enunciado

> Se tiene una matriz $A$ sdp de $10^5\times10^5$. Se sabe que sus autovalores están contenidos en
> $$[1,2]\cup[2000,2001]\cup[4999,5000].$$
> El mínimo autovalor es $1$ y el máximo es $5000$.
>
> **(a)** Calcular el número de condición $\kappa_2(A)$.
>
> **(b)** Si sólo se conociera este número de condición, estimar cuántas iteraciones del método de Richardson harían falta para lograr $\dfrac{\|x_k-x_*\|_A}{\|x_0-x_*\|_A}\leq10^{-5}$.
>
> **(c)** Si sólo se conociera este número de condición, estimar cuántas iteraciones del método de gradiente conjugado harían falta para lograr $\dfrac{\|x_k-x_*\|_A}{\|x_0-x_*\|_A}\leq10^{-5}$.
>
> **(d)** Utilizando la información sobre la distribución de los autovalores, obtener una estimación mejor para el número de iteraciones que harían falta con el gradiente conjugado para lograr $\dfrac{\|x_k-x_*\|_A}{\|x_0-x_*\|_A}\leq10^{-5}$.
> *(Quien obtenga el menor número tendrá 10 puntos más en el examen, siempre que el resultado sea correcto.)*

## Resolución

### (a) Número de condición

Para $A$ sdp, $\kappa_2(A)=\lambda_{\max}(A)/\lambda_{\min}(A)$ (caso particular de $\kappa_2(A)=\sigma_1/\sigma_m$ cuando los valores singulares coinciden con los autovalores, que pasa exactamente cuando $A$ es simétrica - ver [[Numero de condicion y error relativo]]). El enunciado ya da ambos extremos directamente ($\lambda_{\min}=1$, $\lambda_{\max}=5000$, no hace falta ir a buscarlos a los tres bloques):
$$
\kappa_2(A) = \lambda_{\max}/\lambda_{\min} = 5000/1 = 5000.
$$

### (b) Richardson, cota genérica

¿Por qué esta fórmula y no otra? El método de Richardson con el $\alpha$ óptimo es, en este problema (y para cualquier $A$ sdp), exactamente el método de [[Descenso mas pronunciado|descenso más pronunciado]] (SD), y su [[Teorema 5.9|cota de convergencia]] es
$$
\|x_{k+1}-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A \quad\Longrightarrow\quad \|x_k-x_*\|_A \leq r^k\|x_0-x_*\|_A, \quad r:=\frac{\kappa-1}{\kappa+1}.
$$
Con $\kappa=5000$:
$$
r=\frac{4999}{5001}=0.99960008\dots
$$
Se pide $r^k\leq10^{-5}$; tomando logaritmo (y como $\ln r<0$, la desigualdad se da vuelta al dividir):
$$
k\ln r \leq \ln(10^{-5}) \iff k \geq \frac{\ln(10^{-5})}{\ln r}.
$$
Verificado numéricamente: $k=28783$ (con $k=28782$ no alcanza todavía).

### (c) GC, cota genérica

A diferencia de Richardson/SD, el GC no tiene una cota "directa" tan simple - sale de construir un polinomio auxiliar $p$ con $p(0)=1$ que sea chico en $\sigma(A)$ y usar el [[Corolario 5.14]] (ver [[Elegir el polinomio para las cotas de GC]] si no queda claro de dónde sale esto - ahí está la construcción paso a paso). Cuando lo único que se sabe es $\kappa_2(A)$, la elección de $p$ es un polinomio de Chebyshev sobre $[\lambda_{\min},\lambda_{\max}]$, y eso da la cota del [[Teorema 5.19]], ya resuelta para "cuántas iteraciones hacen falta" en la [[Observacion 5.20|Observación 5.20]]:
$$
k\geq\frac{\sqrt{\kappa_2(A)}+1}{2}\ln\frac2\varepsilon.
$$
Con $\kappa_2(A)=5000$, $\varepsilon=10^{-5}$:
$$
k \geq \frac{\sqrt{5000}+1}{2}\ln(2\times10^5) \approx 437.65 \implies k=438.
$$
(Chequeo cruzado, para no confiar en una sola fórmula: aplicando directamente la cota exacta $2\big(\frac{\sqrt\kappa-1}{\sqrt\kappa+1}\big)^k\leq10^{-5}$ de la que sale la Observación 5.20, el menor $k$ que la cumple es también $438$ - coincide.)

**Por qué GC le gana tanto a Richardson acá:** Richardson/SD solo puede mejorar su velocidad recalculando $\alpha$ óptimo en cada paso, lo cual (se prueba en el Teorema 5.9) da como mucho la misma cota que la mejor elección de $p$ de grado $1$. GC, en cambio, minimiza sobre **todos** los polinomios de grado $\leq k$ simultáneamente (es la libertad extra que le da trabajar sobre todo el espacio de Krylov, no solo la dirección del residuo actual) - esa libertad extra es la que transforma la cota de $O(\kappa)$ a $O(\sqrt\kappa)$ (ver [[Cotas de error del gradiente conjugado]]). Numéricamente: **28783 vs. 438**, casi dos órdenes de magnitud.

### (d) Usando la distribución de autovalores (cota mejorada)

El $438$ de (c) sale de elegir $p$ = Chebyshev sobre **todo** $[1,5000]$ - una elección que ignora que el espectro está concentrado en tres bloques angostos (de ancho $1$ cada uno) con huecos enormes entre ellos, en vez de esparcido uniformemente. Un polinomio más inteligente puede aprovechar esos huecos: la técnica general (con la intuición de por qué funciona, ilustrada) está en [[Elegir el polinomio para las cotas de GC#Caso 3: varios bloques aislados|Elegir el polinomio - Caso 3]], que es exactamente generalizar [[5.8 - Cota pesimista autovalor aislado|Problema 5.8]] (un bloque aislado) a **tres** bloques.

La idea: un polinomio de grado $3$ con una raíz en cada bloque, repetido $k$ veces (grado total $3k$), balanceado para que ningún factor individual se dispare al evaluarlo en los otros bloques (ver la nota de por qué *no* conviene usar tres factores independientes, en la página linkeada arriba). **Importante:** esta construcción funciona bien para *estos* centros en particular, pero no es una garantía universal - con otras posiciones de bloques puede quedarse estancada o incluso explotar (ver [[Elegir el polinomio para las cotas de GC#Receta: qué hacer según lo que se sabe de $\sigma(A)$|la advertencia completa, con contraejemplos]]), por eso abajo se verifica numéricamente bloque por bloque en vez de solo confiar en la fórmula:
$$
p_{3k}(z) = \left(\frac{(c_1-z)(c_2-z)(c_3-z)}{c_1c_2c_3}\right)^{k}, \qquad c_1=1.5,\ c_2=2000.5,\ c_3=4999.5\ \text{(centros de cada bloque)}.
$$
Por el Corolario 5.14, $\|x_{3k}-x_*\|_A/\|x_0-x_*\|_A \leq \max_{z\in\sigma(A)}|p_{3k}(z)|$, así que alcanza con encontrar el $k$ más chico que hace esa cota $\leq10^{-5}$:
```julia
z1 = range(1,2,length=2000)
z2 = range(2000,2001,length=2000)
z3 = range(4999,5000,length=2000)
q(z) = (1.5-z)*(2000.5-z)*(4999.5-z) / (1.5*2000.5*4999.5)
maxq = max(maximum(abs.(q.(z1))), maximum(abs.(q.(z2))), maximum(abs.(q.(z3))))   # -> 0.4997
k = ceil(Int, log(1e-5)/log(maxq))    # -> 17  =>  grado total 3k = 51
```
Con los centros exactos, $\max|q|\approx0.4997$ (el máximo se alcanza cerca del bloque 3, el más lejano - ver el gráfico en la página linkeada), y alcanza con $k=17$ repeticiones (**grado $3k=51$**) para bajar de $10^{-5}$ - verificado. Moviendo un poco los puntos $c_1,c_2,c_3$ dentro de cada bloque (búsqueda numérica fina en una grilla) se llega a $\max|q|\approx0.428$ con $(c_1,c_2,c_3)=(1.75,\,2000.55,\,4999.5)$, lo que baja a $k=14$ (**grado $42$**) - no necesariamente el óptimo global, pero ya una mejora adicional sobre los centros "obvios".

> [!idea]
> **Dos atajos que usa el script, justificados.**
>
> **1. Por qué $\max_z|p_{3k}(z)|=(\max_z|q(z)|)^k$, sin recalcular el máximo para cada $k$.** Como $p_{3k}(z)=q(z)^k$ y elevar a una potencia fija $k$ es una función creciente de $|q(z)|\geq0$, el punto $z$ donde $|q(z)|$ es máximo es **el mismo** donde $|q(z)|^k$ es máximo - elevar a la $k$ no mueve el máximo, solo cambia el valor. Por eso alcanza con calcular $\max|q(z)|$ **una sola vez** y después despejar $k$ algebraicamente con logaritmo, en vez de recorrer la grilla una vez por cada $k$ candidato:
> $$(\max q)^k\leq10^{-5} \iff k\log(\max q)\leq\log(10^{-5}) \iff k\geq\frac{\log(10^{-5})}{\log(\max q)}$$
> (la desigualdad se da vuelta al dividir porque $\log(\max q)<0$, ya que $0<\max q<1$). Verificado numéricamente para $k=1,5,17$: `maximum(abs.(q.(z).^k))` da exactamente lo mismo que `maxq^k`, a precisión de máquina.
>
> **2. Por qué evaluar en una grilla (y no analíticamente en todo el dominio continuo) ya da el máximo exacto, no una aproximación.** $q$ es un polinomio **cúbico** en $z$, así que $q'$ es cuadrática, con a lo sumo dos raíces reales, calculables en forma cerrada. Para $(c_1,c_2,c_3)=(1.5,\,2000.5,\,4999.5)$:
> $$q'(z)=0 \iff 3z^2-2(c_1+c_2+c_3)z+(c_1c_2+c_1c_3+c_2c_3)=0 \implies z\approx881.44 \ \text{ó}\ z\approx3786.23.$$
> Ninguno de los dos cae **dentro** de ninguno de los tres bloques ($[1,2]$, $[2000,2001]$, $[4999,5000]$) - y se repite para los centros optimizados $(1.75,2000.55,4999.5)$ y para otros centros probados: los puntos críticos de $q$ apenas se mueven (quedan siempre cerca de $881$ y $3786$) porque dependen continuamente de $c_1,c_2,c_3$, que solo se mueven dentro de bloques angostos. Como no hay ningún punto crítico adentro, $q$ es **monótona** en cada bloque, así que su máximo continuo está necesariamente en uno de los dos bordes - y como la grilla `range(a,b,length=2000)` incluye esos bordes exactos como primer y último punto, el máximo de la grilla **coincide exactamente** (no aproximadamente) con el máximo continuo real. Confirmado calculando ambos por separado: dan el mismo valor a precisión de máquina en los tres bloques, para los tres juegos de centros probados en esta página.
>
> Esto es específico de esta construcción (un solo factor cúbico, repetido) - con un polinomio de grado más alto o más factores, los puntos críticos podrían no ser tan fáciles de ubicar a mano, y ahí sí conviene una grilla bien fina (o una optimización numérica) en lugar de confiar en que los bordes alcanzan.

**Conclusión: $42$ (o $51$ sin optimizar la elección de los centros) en vez de $438$** - un orden de magnitud menos, por la misma razón que en el Problema 5.8: la cota genérica solo ve $\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$, ciega a que casi todo ese rango $[1,5000]$ está **vacío** de autovalores; explotar la estructura de bloques aislados cuesta apenas un factor constante en el grado (acá, $\times3$ por repetición en vez de $\times1$) a cambio de una ganancia enorme en la tasa de decrecimiento por iteración.

> [!pendiente]
> No se verificó que $42$ sea el mínimo absoluto posible - el enunciado da puntos extra a quien encuentre el número más chico, lo cual es en sí un problema de optimización min-max que admite mejores soluciones (p. ej. con polinomios de Chebyshev genuinos sobre cada bloque en vez de factores lineales a potencia, en vez de un simple grid search sobre los centros). Lo importante para entender el método es que alcanza con un factor constante en el grado por cada bloque aislado que se anula explícitamente - no hace falta encontrar el óptimo exacto para capturar la mejora de fondo.

## Ver también
- [[Parcial Modelo - Indice|Parcial Modelo]]
- [[Elegir el polinomio para las cotas de GC]]
- [[Cotas de error del gradiente conjugado]]
- [[5.8 - Cota pesimista autovalor aislado|Problema 5.8]]
- [[Numero de condicion y error relativo]]
