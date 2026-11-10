---
tags: [moc, comparacion]
---

# Índice de tablas comparativas

## Consulta rápida (todos los capítulos)
- [[Tablas Comparativas - Formulario de cotas e iteraciones]] - glosario de variables, cómo calcular $\kappa_2(A)$ según el caso, tasa de convergencia e iteraciones necesarias por método, error/residuo dado $k$, relación error↔residuo, qué hipótesis hace falta para que cada método converja (necesaria vs. suficiente, con contraejemplos), y costo por iteración/memoria de todos los métodos juntos. Pensada para las consignas tipo "estimar cuántas iteraciones hacen falta" / "calcular $\kappa_2$" / "¿cuál es el error después de $k$ pasos?" / "demostrar que converge" que abundan en los parciales.

## Capítulo 2 - Métodos directos
- [[Costo computacional metodos directos]]

## Capítulo 3 - Factorización QR
- [[Estabilidad y numero de condicion]] (ecuaciones normales vs. [[Ortogonalizacion de Gram-Schmidt|Gram-Schmidt]] vs. [[Matrices de Householder|Householder]])

## Capítulo 4 - Métodos iterativos estacionarios
- [[Comparacion metodos iterativos estacionarios]] ([[Metodo de Jacobi|Jacobi]], [[Metodo de Gauss-Seidel|Gauss-Seidel]], [[Metodo SOR|SOR]], [[Metodo de Richardson|Richardson]])

## Capítulo 5 - Métodos de descenso
- [[Tablas Comparativas - GC vs SD]]
- [[Tablas Comparativas - GC vs GCP]]
- [[Tablas Comparativas - GC vs GMRes]] (incluye costo por iteración)
- [[Precondicionadores genericos]] (tabla al final: diagonal, GS simétrico, Cholesky incompleta)
- [[Tablas Comparativas - Precondicionadores]] - los tres genéricos vs. multigrid, con ventajas y desventajas de cada uno
- [[GMRes - precondicionamiento]] (tabla: izquierda vs. derecha)

## Capítulo 6 - Sistemas no lineales
- [[Tablas Comparativas - Punto fijo vs Newton]] (incluye Newton simplificado y Nelder–Mead, y una tabla síntoma → causa)

## Capítulo 7 - Autovalores
- [[Tablas Comparativas - Metodos para autovalores]] (potencias, inversas, Rayleigh, QR; números de los problemas en $A_N$; hermítica vs. general)
- Sección 8 de [[Tablas Comparativas - Formulario de cotas e iteraciones]]: tasas y órdenes de los Capítulos 6–7.

## La tabla de todas las tablas: qué método elegir

| Situación | Método recomendado | Por qué |
|---|---|---|
| Sistema chico, denso, $A$ genérica | LU con [[Pivoteo y matrices de permutacion\|pivoteo]] ([[Factorizacion LU]]) | directo, $O(n^3)$, robusto |
| Sistema chico, denso, $A$ sdp | [[Factorizacion de Cholesky\|Cholesky]] | la mitad de operaciones que [[Factorizacion LU\|LU]] |
| Mínimos cuadrados, $A$ bien condicionada | Ecuaciones normales ([[Minimos cuadrados y ecuaciones normales]]) | más barato |
| Mínimos cuadrados, $A$ mal condicionada | QR por Householder ([[Matrices de Householder]]) | estable ($κ$ vs. $κ²$) |
| Sistema grande y sparse, $A$ sdp, sin estructura especial | Gradiente conjugado, quizás con Jacobi/Cholesky incompleta ([[Implementacion del gradiente conjugado]]) | $O(\sqrt\kappa)$ iteraciones |
| Sistema grande y sparse, $A$ sdp, viene de una EDP | GC + multigrid como precondicionador ([[Multigrid en 1D - transferencia e implementacion]]) | $O(1)$ iteraciones, independiente de $h$ |
| Sistema grande, $A$ no simétrica | [[GMRes - definicion y propiedad de minimizacion\|GMRes]], con precondicionamiento si hace falta | minimiza el residuo, no exige simetría |
| Sistema **no lineal**, probar existencia/unicidad en una región | punto fijo + [[Criterio de contraccion via la derivada\|contracción]] (Teo. 6.12) | da el resultado teórico y una cota a priori |
| Sistema **no lineal**, calcular la raíz con precisión | [[Metodo de Newton - definicion y buena definicion\|Newton]] desde un $x_0$ cercano | cuadrático; cada paso es un sistema lineal |
| Todos los autovalores, matriz densa | [[Algoritmo QR - Hessenberg, shifts y deflacion\|QR]] con Hessenberg + shift de Wilkinson + deflación (`eigvals`) | $O(n^3)$ total, robusto |
| Autovalor dominante de una matriz enorme y rala | [[Metodo de las potencias\|potencias]] (si $\vert \lambda_2/\lambda_1\vert $ es chico) | 1 mat-vec por paso |
| Un autovalor cerca de un valor conocido | [[Metodo de las potencias inversas\|potencias inversas]] con ese shift | tasa $\propto\vert \lambda-\mu\vert $, la elegimos |

Volver a [[Inicio]].
