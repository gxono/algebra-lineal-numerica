---
tags: [moc, parcial]
---

# Parcial Modelo 2026 - Primer Parcial

Examen modelo del primer parcial (Cálculo Numérico II / Análisis Numérico, 2026), con enunciado completo y resolución justificada, verificada por ejecución (no solo leída) - mismo estándar que el resto del vault.

Fuente del enunciado: `algebra-lineal-numerica/parciales/modelo-parcial-1.pdf`. La resolución de 1(a)/1(b) partió de un borrador propio (`resolucion-parcial1.typ`, `script-parcial1.jl`); el resto se completó acá.

## Problemas

- [[Parcial Modelo - Problema 1|Problema 1]] (15 pts) - normas de una matriz triangular, convergencia de $x_{k+1}=Mx_k$.
- [[Parcial Modelo - Problema 2|Problema 2]] (20 pts) - matriz de Läuchli: ecuaciones normales vs. QR (Householder y Gram-Schmidt clásico) en un caso extremo de mal condicionamiento.
- [[Parcial Modelo - Problema 3|Problema 3]] (30 pts) - estimaciones de iteraciones para Richardson y GC a partir del espectro, con y sin explotar autovalores aislados.
- [[Parcial Modelo - Problema 4|Problema 4]] (35 pts) - el problema modelo $-u''+u=f$: Gauss-Seidel vs. descenso vs. GC, experimentos numéricos.

## Qué conceptos puso en juego (y dónde ampliarlos)

| Problema | Conceptos centrales | Si hace falta repasar |
|---|---|---|
| 1 | Normas matriciales, radio espectral, convergencia de $M^k$ | [[Normas de vectores y matrices]], [[Espectro y radio espectral]] |
| 2 | Mínimos cuadrados, ecuaciones normales, QR (Householder y Gram-Schmidt), número de condición, estabilidad numérica | [[Minimos cuadrados y ecuaciones normales]], [[Matrices de Householder]], [[Ortogonalizacion de Gram-Schmidt]], [[Estabilidad y numero de condicion]] |
| 3 | Número de condición, cotas de Richardson y GC, autovalores aislados | [[Numero de condicion y error relativo]], [[Cotas de error del gradiente conjugado]], [[5.8 - Cota pesimista autovalor aislado\|Problema 5.8]] |
| 4 | Gauss-Seidel y e.d.d., número de condición vía autovalores explícitos, SD/GC, terminación finita | [[Metodo de Gauss-Seidel]], [[GC como metodo directo]], [[Multigrid en 1D - problema modelo y suavizado]] |

## Ver también
- [[Inicio]]
- [[Glosario - Indice|Glosario]]
- [[Demostraciones - Indice]]
