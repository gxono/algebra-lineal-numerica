---
tags: [capitulo5, comparacion]
---

# Precondicionadores: genéricos vs. multigrid

> [!idea]
> Distinción central que hay que tener clara: los precondicionadores **genéricos** (diagonal, Gauss-Seidel simétrico, Cholesky incompleta) solo miran las *entradas* de $A$, y mejoran la **constante** $c_2/c_1$ pero no el **orden** - si $\kappa_2(A)=O(h^{-2})$, siguen dando $c_2/c_1=O(h^{-2})$. [[Multigrid en 1D - transferencia e implementacion|Multigrid]] es el único que usa la *estructura* del problema (que viene de una EDP discretizada), y por eso es el único que logra $c_2/c_1=O(1)$, **independiente de $h$**.

## Tabla comparativa completa

> [!comparacion]
> | Precondicionador | Acción $\xi\to M\xi$ | Preparación | ¿Mejora el orden de $\kappa$? | Ventaja principal | Desventaja principal |
> |---|---|---|---|---|---|
> | Ninguno ($M=I$) | gratis | ninguna | - | nada que preparar | $\kappa_2(A)$ sin cambios |
> | Diagonal (Jacobi) | $n$ divisiones | ninguna | No | el más barato posible | inútil si la diagonal es ~constante (caso típico en discretizaciones uniformes) |
> | Gauss-Seidel simétrico | 2 triangulares + 1 diagonal | ninguna | No | mejor que diagonal, sigue sin preparación | todavía $O(h^{-2})$, solo mejora la constante |
> | Cholesky incompleta | 2 triangulares | factorizar 1 vez (costo no trivial) | No | el mejor de los genéricos en la práctica | preparación cara; degenerado en 1D (coincide con la completa) |
> | Multigrid (ciclo en V) | $O(n)$ | construir jerarquía de mallas | **Sí** | $\kappa$ acotado independiente de $h$ - el único que logra esto | mucho más complejo de implementar; exige que el problema tenga estructura de grilla/EDP |

## Por qué los genéricos no alcanzan (y multigrid sí)

Los precondicionadores genéricos reciclan un método estacionario del Capítulo 4 ($M$ tal que $I-MA$ sea la matriz de iteración de Jacobi/GS - ver [[Precondicionadores genericos]]): están gobernados por los mismos modos **suaves** que hacen lentos a esos métodos en primer lugar. Multigrid ataca exactamente ese punto ciego: lo que es suave en la malla fina se ve oscilatorio en una malla más gruesa, así que la combinación suavizado + corrección en malla gruesa elimina **todos** los modos, no solo los oscilatorios - ver [[Multigrid en 1D - problema modelo y suavizado]].

## Evidencia numérica

Problema modelo 1D, $n$ de 31 a 4095, GC sin precondicionar vs. GC + multigrid (código completo y reproducible en [[Multigrid en 1D - transferencia e implementacion]]):

![[cap5-multigrid-independencia-h.png]]

GC se duplica cada vez que $n$ se duplica (16→2048, lineal). GC+multigrid se queda en 9-10 iteraciones en **todo** el rango - la firma numérica de "independiente de $h$", sin necesidad de creerle a la cota teórica.

Consistente con la tabla real de la cátedra ([[5.14 - Comparacion de los cuatro metodos|Problema 5.14]]), que además compara contra GCP-Jacobi (crece igual que GC puro, confirmando que Jacobi no cambia el orden):

> [!comparacion]
> | $n$ | GC (sin precond.) | GCP-Jacobi | Multigrid solo | GCP-Multigrid |
> |---|---|---|---|---|
> | 31 | 16 | 16 | 15 | 6-8 |
> | 511 | 256 | 256 | 16 | 6-8 |
> | 4095 | 2049 | 2049 | 16 | 5-6 |

En el laplaciano 2D (Problema 5.7/5.11), comparación de los tres precondicionadores genéricos entre sí:

![[cap5-precondicionadores.png]]

Cholesky incompleta gana a Gauss-Seidel simétrico, que gana a diagonal - pero **ninguno** de los tres cambia el orden de convergencia, solo la constante (y en este caso particular, con diagonal constante, Jacobi ni siquiera mejora la constante: queda exactamente igual que no precondicionar).

## Cuándo usar cada uno

- **Sistema chico o $\kappa_2(A)$ ya razonable:** no hace falta precondicionar.
- **Sistema grande, sparse, sin estructura particular conocida:** Cholesky incompleta (el mejor genérico) o Gauss-Seidel simétrico si la preparación de Cholesky incompleta es muy cara.
- **Sistema grande, sparse, viene de discretizar una EDP:** multigrid, sin dudarlo - es la única opción que no degrada con $h$.
- **Diagonal:** prácticamente nunca vale la pena salvo que $A$ tenga escalas muy distintas en la diagonal (poco común en los problemas del curso, donde suele ser constante).

## Ver también
- [[Precondicionadores genericos]] (fórmulas y código de diagonal/SGS/Cholesky incompleta)
- [[Multigrid en 1D - transferencia e implementacion]] (ciclo en V, implementación)
- [[Precondicionamiento y GCP]] (cómo se usa $M$ dentro de GC)
- [[Tablas Comparativas - GC vs GCP]]
- [[Tablas Comparativas - Indice]]
- [[Referencia de Julia#IterativeSolvers.jl y Preconditioners.jl: métodos y precondicionadores de fábrica|Referencia de Julia]] - mismo patrón (genéricos vs. AMG), confirmado con una implementación de paquete independiente de la del vault
