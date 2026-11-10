---
tags: [parcial, capitulo3]
---

# Parcial Modelo - Problema 2 (20 pts): la matriz de Läuchli

## Enunciado

> Para $\varepsilon>0$ se considera la matriz de Läuchli
> $$A(\varepsilon) = \begin{bmatrix}1&1&1\\\varepsilon&0&0\\0&\varepsilon&0\\0&0&\varepsilon\end{bmatrix}\in\mathbb{R}^{4\times3}, \qquad b=\begin{bmatrix}1\\0\\0\\0\end{bmatrix}.$$
> **(a)** Plantee las ecuaciones normales $A(\varepsilon)^TA(\varepsilon)\,x=A(\varepsilon)^Tb$ y resuélvalas *simbólicamente* (en función de $\varepsilon$). Verifique que la solución exacta de $\min\|A(\varepsilon)x-b\|_2$ no depende de $\varepsilon$.
>
> **(b)** Calcule $\kappa_2(A(\varepsilon))$ para $\varepsilon=10^{-2},10^{-4},10^{-6},10^{-7},10^{-8},10^{-9}$. ¿Cómo depende de $\varepsilon$?
>
> **(c)** Para cada uno de esos valores de $\varepsilon$, resuelva el problema de cuadrados mínimos de tres maneras: ecuaciones normales (intentando la factorización de Cholesky de $A^TA$), QR por Householder y QR por Gram-Schmidt clásico. Arme una tabla con el error relativo respecto de la solución exacta del ítem (a), e indique en qué casos algún método no llega a producir una respuesta.
>
> **(d)** Para $\varepsilon=10^{-8}$, imprima la matriz $A^TA$ calculada por la computadora y compárela con la que obtuvo simbólicamente en (a). ¿Qué le pasó? Evalúe también en la computadora la expresión `1 + eps^2 == 1`. Relacione lo que observa con lo que ocurre en la tabla del ítem (c).
>
> **(e)** A partir de la tabla, ordene los tres métodos por robustez y explique la diferencia distinguiendo qué parte del error se debe al condicionamiento del problema y qué parte al algoritmo empleado.

Es el ejemplo clásico (Läuchli, 1961) para mostrar *en un caso minúsculo* lo que [[3.6 - Estabilidad experimental minimos cuadrados|Problema 3.6]] muestra con una Vandermonde: las [[Minimos cuadrados y ecuaciones normales|ecuaciones normales]] pueden colapsar por completo, no solo perder un par de cifras.

## Resolución

### (a) Ecuaciones normales simbólicas

El enunciado pide dos cosas encadenadas: primero plantear y resolver $A^TAx=A^Tb$ (las [[Minimos cuadrados y ecuaciones normales|ecuaciones normales]] del problema de cuadrados mínimos $\min\|A(\varepsilon)x-b\|_2$ - la vía estándar para resolverlo, aunque acá el punto del problema es mostrar que esa vía falla numéricamente más adelante), y después verificar la afirmación "no depende de $\varepsilon$". Primero lo algebraico:

$$
A^TA = \begin{bmatrix}1+\varepsilon^2&1&1\\1&1+\varepsilon^2&1\\1&1&1+\varepsilon^2\end{bmatrix}, \qquad A^Tb=\begin{bmatrix}1\\1\\1\end{bmatrix}.
$$
Por simetría, la solución de $A^TA\,x=A^Tb$ es $x=(c,c,c)^T$ con $(1+\varepsilon^2)c+2c=1$, es decir
$$
x(\varepsilon) = \frac{1}{3+\varepsilon^2}(1,1,1)^T.
$$
Esto **sí depende de $\varepsilon$** en aritmética exacta - pero de forma $O(\varepsilon^2)$, así que para los $\varepsilon\le10^{-2}$ de este problema la dependencia es inapreciable: en $\varepsilon=10^{-2}$ ya da $(0.333322,\dots)$, y para $\varepsilon\le10^{-6}$ coincide con $(1/3,1/3,1/3)^T$ hasta precisión de máquina (verificado: `A(ε)\b` da exactamente esto para todo $\varepsilon$ de la lista). Esa es la "solución exacta" de referencia para (c).

### (b) $\kappa_2(A(\varepsilon))$

$A^TA=\varepsilon^2I+\mathbb{1}\mathbb{1}^T$ (con $\mathbb{1}=(1,1,1)^T$), y $\mathbb{1}\mathbb{1}^T$ tiene autovalores $\{3,0,0\}$ (rango 1, traza 3). Como $I$ y $\mathbb{1}\mathbb{1}^T$ conmutan, los autovalores de $A^TA$ son $\{3+\varepsilon^2,\varepsilon^2,\varepsilon^2\}$, así que
$$
\kappa_2(A(\varepsilon)) = \sqrt{\frac{3+\varepsilon^2}{\varepsilon^2}} = \frac{\sqrt{3+\varepsilon^2}}{\varepsilon} \approx \frac{\sqrt3}{\varepsilon}.
$$
Crece como $O(1/\varepsilon)$ - calculado para los 6 valores pedidos, con la fórmula y verificado numéricamente contra `cond(A(ε))`:

| $\varepsilon$ | $\kappa_2(A(\varepsilon))$ (fórmula) | `cond(A(ε))` |
|---|---|---|
| $10^{-2}$ | $173.208$ | $173.208$ |
| $10^{-4}$ | $1.73205\times10^4$ | $1.73205\times10^4$ |
| $10^{-6}$ | $1.73205\times10^6$ | $1.73205\times10^6$ |
| $10^{-7}$ | $1.73205\times10^7$ | $1.73205\times10^7$ |
| $10^{-8}$ | $1.73205\times10^8$ | $1.73205\times10^8$ |
| $10^{-9}$ | $1.73205\times10^9$ | $1.73205\times10^9$ |

(coinciden en los 6 casos, a precisión de máquina)

### (c) Comparación de los tres métodos

Implementando `mincuad_normales` (ecuaciones normales + Cholesky, como en [[3.6 - Estabilidad experimental minimos cuadrados|Problema 3.6]]), `qrhouse`/`qrsolve` ([[Matrices de Householder]]) y una versión de Gram-Schmidt **clásico** (no modificado: todas las proyecciones usan las columnas originales de $A$, a diferencia de `gram_schmidt!` del vault que es la versión modificada de [[Ortogonalizacion de Gram-Schmidt]]):

```julia
using LinearAlgebra

function Aeps(eps)
    A = zeros(4,3)
    A[1,:] .= 1.0
    A[2,1]=eps
    A[3,2]=eps
    A[4,3]=eps
    return A
end
b = [1.0,0,0,0]

mincuad_normales(A,b) = cholesky(Symmetric(A'A)) \ (A'b)

function gram_schmidt_clasico(A)
    n,m = size(A)
    Q = zeros(n,m)
    R = zeros(m,m)
    for j in 1:m
        v = A[:,j]
        for i in 1:j-1
            R[i,j] = Q[:,i]'A[:,j]      # proyección sobre la columna ORIGINAL, no la corregida
            v -= R[i,j]*Q[:,i]
        end
        R[j,j] = norm(v)
        Q[:,j] = v/R[j,j]
    end
    return Q, R
end
function mincuad_gs_clasico(A,b)
    Q, R = gram_schmidt_clasico(A)
    return R \ (Q'b)
end
```

> [!comparacion]
> | $\varepsilon$ | ec. normales | Householder | Gram-Schmidt clásico |
> |---|---|---|---|
> | $10^{-2}$ | $1.5\times10^{-12}$ | $3.0\times10^{-16}$ | $2.6\times10^{-12}$ |
> | $10^{-4}$ | $5.6\times10^{-9}$ | $7.0\times10^{-16}$ | $3.2\times10^{-9}$ |
> | $10^{-6}$ | $5.6\times10^{-13}$ | $2.4\times10^{-16}$ | $7.3\times10^{-5}$ |
> | $10^{-7}$ | $0.0140$ | $7.3\times10^{-16}$ | $0.0188$ |
> | $10^{-8}$ | **falla** (Cholesky) | $6.4\times10^{-16}$ | $1.414$ (respuesta falsa) |
> | $10^{-9}$ | **falla** (Cholesky) | $6.2\times10^{-16}$ | $1.414$ (respuesta falsa) |

(error relativo $\|x_{\text{método}}-x(\varepsilon)\|_2/\|x(\varepsilon)\|_2$, verificado corriendo los tres métodos.)

**Householder se mantiene en precisión de máquina en los seis casos**, sin importar $\kappa_2(A)$. Las ecuaciones normales se degradan y a partir de $\varepsilon=10^{-8}$ la factorización de Cholesky de $A^TA$ **directamente falla** (`PosDefException`, ver (d)). Gram-Schmidt clásico es el peor: para $\varepsilon\le10^{-8}$ no falla con un error visible, **da una respuesta completamente incorrecta sin avisar** ($\|x-x(\varepsilon)\|_2/\|x(\varepsilon)\|_2\approx\sqrt2$, es decir, ni siquiera está cerca).

### (d) Qué pasa en la computadora con $\varepsilon=10^{-8}$

```julia
A = Aeps(1e-8)
A'A    # -> [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]  (exactamente!)
1 + (1e-8)^2 == 1    # -> true
```
$\varepsilon^2=10^{-16}$ es **menor que la mitad del épsilon de máquina** ($\texttt{eps}(\text{Float64})/2\approx1.11\times10^{-16}$), así que $1+\varepsilon^2$ redondea exactamente a $1$ en punto flotante. Por eso la diagonal de $A^TA$ calculada por la computadora da $1$, no $1+\varepsilon^2$: la matriz calculada es **exactamente** $\mathbb{1}\mathbb{1}^T$ (rango 1, singular, ni siquiera semidefinida positiva en el sentido estricto que exige Cholesky) - de ahí el `PosDefException` en (c). El problema, que en aritmética exacta tiene $\kappa_2(A^TA)=\kappa_2(A)^2\approx3/\varepsilon^2=3\times10^{16}$ (ya del orden de $1/\varepsilon_{\text{maq}}$), se vuelve **numéricamente indistinguible de un problema singular** apenas se forma $A^TA$ - esto es exactamente $\kappa_2(A^TA)=\kappa_2(A)^2$ llevado al extremo (ver [[Estabilidad y numero de condicion]]).

### (e) Ranking de robustez

**Householder > ecuaciones normales > Gram-Schmidt clásico.**

- **Householder** usa solo reflexiones (isometrías), nunca forma $A^TA$: el error se mantiene en $\kappa_2(A)\cdot\varepsilon_{\text{maq}}$ en el peor caso, y acá ni siquiera llega a notarse - **estable por el algoritmo**, independiente de qué tan mal condicionado esté el problema en sí.
- **Ecuaciones normales**: el error crece como $\kappa_2(A)^2\varepsilon_{\text{maq}}$ porque el algoritmo mismo (formar $A^TA$) **eleva al cuadrado** el condicionamiento del problema original - es **inestabilidad del algoritmo**, no del problema (el problema de cuadrados mínimos en sí tiene $\kappa_2(A)\approx\sqrt3/\varepsilon$, no $\kappa_2(A)^2$).
- **Gram-Schmidt clásico** es el más inestable de los tres: pierde ortogonalidad entre columnas mucho antes que la versión modificada (ver [[Ortogonalizacion de Gram-Schmidt]]) porque cada proyección usa la columna *original* en vez de la ya corregida, acumulando errores de cancelación - y lo hace **sin señal de alarma** (no lanza una excepción como Cholesky), lo cual es más peligroso en la práctica: una respuesta visiblemente mala es mejor que una respuesta silenciosamente mala.

La distinción pedida en el enunciado: el **condicionamiento del problema** ($\kappa_2(A)\approx\sqrt3/\varepsilon$) es el mismo para los tres métodos - es una propiedad de $A$ y $b$, no del algoritmo. Lo que cambia entre métodos es **cuánto error adicional introduce el algoritmo** encima de ese piso inevitable: Householder prácticamente no agrega nada; ecuaciones normales agrega un factor $\kappa_2(A)$ extra (porque trabaja con $\kappa_2(A)^2$ en vez de $\kappa_2(A)$); Gram-Schmidt clásico agrega aún más por la pérdida de ortogonalidad.

## Ver también
- [[Parcial Modelo - Indice|Parcial Modelo]]
- [[Minimos cuadrados y ecuaciones normales]]
- [[Matrices de Householder]]
- [[Ortogonalizacion de Gram-Schmidt]]
- [[Estabilidad y numero de condicion]]
- [[3.6 - Estabilidad experimental minimos cuadrados|Problema 3.6]] (el mismo fenómeno con una Vandermonde)
