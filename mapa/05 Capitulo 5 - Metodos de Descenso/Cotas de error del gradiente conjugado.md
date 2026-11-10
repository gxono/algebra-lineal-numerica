---
tags: [capitulo5, teorema, definicion]
---

# Cotas de error del [[Gradiente conjugado - definicion y propiedades|gradiente conjugado]]

> [!idea]
> El resultado más importante del capítulo: $x_k$ es, entre **todos** los puntos de $x_0+\mathcal{K}_k$, el que minimiza el error en norma $A$ - y eso se puede escribir como "el mejor polinomio residual evaluado en $A$". Eligiendo un polinomio de Chebyshev se obtiene la cota clásica $O(\sqrt{\kappa_2(A)})$ iteraciones, muy superior al $O(\kappa_2(A))$ de steepest descent.

## Teorema 5.13

$$
\|x_k-x_*\|_A = \min_{p\in\overline{\mathcal P}_k}\|p(A)(x_0-x_*)\|_A, \qquad \overline{\mathcal P}_k:=\{p\in\mathcal P_k:p(0)=1\}.
$$

Demostración completa (reescribiendo $x_0+\mathcal K_k$ vía el [[Funcional de energia y caracterizacion|Lema 5.3]]): [[Teorema 5.13]].

## Corolario 5.14 (cota espectral)

$$
\|x_k-x_*\|_A \leq \Big(\max_{\lambda\in\sigma(A)}|p(\lambda)|\Big)\|x_0-x_*\|_A, \qquad \forall p\in\overline{\mathcal P}_k.
$$

Demostración completa (en dos pasos: $\|p(A)\|_A\leq\|p(A)\|_2$, y $\|p(A)\|_2\leq\max_{\lambda\in\sigma(A)}|p(\lambda)|$ por diagonalización): [[Corolario 5.14]].

> [!idea]
> **Estrategia general** (usada en todos los ejemplos siguientes): elegir **cualquier** $p\in\overline{\mathcal P}_k$ que sea chico sobre el espectro de $A$, sin preocuparse de qué hace fuera de él. El GC, sin "saberlo", siempre hace al menos tan bien como la mejor elección de $p$ - por eso todas estas cotas son cotas *superiores*, nunca el comportamiento exacto.

## Ejemplos de uso (eligiendo $p$ según lo que se sabe de $σ(A)$)

Versión extendida, con la intuición paso a paso de *por qué* cada elección funciona e ilustraciones de cómo se ve $p(z)$ en cada caso: [[Elegir el polinomio para las cotas de GC]].

- **Espectro agrupado** $\sigma(A)\subset[1-\varepsilon,1+\varepsilon]$: con $p(z)=(1-z)^k$, $\|x_k-x_*\|_A\leq\varepsilon^k\|x_0-x_*\|_A$ - convergencia muy rápida, **independiente de $n$**.
- **Pocos autovalores distintos** ($m$ de ellos): con $p(z)=\prod_j(1-z/\mu_j)$, $x_m=x_*$ exacto (ver [[GC como metodo directo]]).
- **Autovalores aislados** ($σ(A)\subset[a,b]\cup\{\mu_1,\dots,\mu_m\}$): combinar un factor que anule los $\mu_j$ con un polinomio de Chebyshev sobre $[a,b]$. La velocidad queda gobernada por $b/a$, **no** por $\kappa_2(A)=\mu_{\max}/a$ - un autovalor aislado "cuesta" solo $m$ iteraciones extra, no arruina la tasa asintótica.

## Teorema 5.19 (cota general vía Chebyshev)

$$
\|x_k-x_*\|_A \leq 2\Big(\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\Big)^k\|x_0-x_*\|_A.
$$

Demostración completa (construyendo $p_k$ con polinomios de Chebyshev y aplicando el Corolario 5.14): [[Teorema 5.19]].

Comparar con steepest descent: $\frac{\kappa_2(A)-1}{\kappa_2(A)+1}$ - la diferencia entre $\kappa$ y $\sqrt\kappa$ en el exponente es la mejora clave (ver [[Tablas Comparativas - GC vs SD]]).

## Observación 5.20 (cuántas iteraciones hacen falta)

Para $\|x_k-x_*\|_A\leq\varepsilon\|x_0-x_*\|_A$ alcanza con
$$
k \geq \frac{\sqrt{\kappa_2(A)}+1}{2}\ln\frac2\varepsilon = O\big(\sqrt{\kappa_2(A)}\big).
$$

Demostración completa: [[Observacion 5.20|Observación 5.20]].

La dependencia en la tolerancia es solo logarítmica: pedir una cifra más de precisión cuesta un número *fijo* de iteraciones; duplicar $\sqrt\kappa$ duplica el total.

## Ver también
- [[Elegir el polinomio para las cotas de GC]]
- [[Gradiente conjugado - definicion y propiedades]]
- [[Descenso mas pronunciado]]
- [[GC como metodo directo]]
