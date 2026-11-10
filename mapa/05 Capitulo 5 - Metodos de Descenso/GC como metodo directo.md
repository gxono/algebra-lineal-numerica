---
tags: [capitulo5, teorema]
---

# El [[Gradiente conjugado - definicion y propiedades|gradiente conjugado]] como método directo

> [!idea]
> Terminación finita: el GC resuelve $Ax=b$ **exactamente** en a lo sumo $n$ pasos (y a menudo mucho antes). En la práctica esto casi nunca se observa (el redondeo destruye la ortogonalidad de los residuos antes de llegar tan lejos), así que el GC se usa como método *iterativo*, con un criterio de parada por tolerancia - pero el resultado teórico es lo que garantiza que *algo* bueno tiene que pasar en a lo sumo $n$ pasos.

## Teorema 5.21

Si $A$ tiene $k$ autovalores distintos, el GC converge en a lo sumo $k$ pasos, partiendo de cualquier $x_0$. Demostración completa: [[Teorema 5.21]].

## Teorema 5.22

Si $b$ es combinación lineal de $k$ autovectores de $A$, el GC partiendo de $x_0=0$ converge en a lo sumo $k$ pasos. Demostración completa (mismo tipo de argumento que el Teorema 5.21, pero aprovechando que $x_0-x_*$ vive en un subespacio invariante chico - comparar con [[5.20 - Demostrar Corolario 5.41|el resultado análogo para GMRes, Corolario 5.41]]): [[Teorema 5.22]].

> [!idea]
> Más fuerte que el anterior: si la solución "vive" en un subespacio invariante chico (aunque $A$ tenga muchos autovalores distintos), alcanza con $k$ pasos.

## Por qué esto no se observa en la práctica

Los errores de redondeo destruyen gradualmente la ortogonalidad exacta de los residuos (ver también [[GMRes - implementacion y rotaciones de Givens|el mismo fenómeno en GMRes con Gram-Schmidt]]). Para $n$ grande ($10^5$, $10^6$), de todos modos no se pretende hacer $n$ iteraciones: se usa un criterio de parada por tolerancia (ver [[Implementacion del gradiente conjugado]]).

## Ver también
- [[Cotas de error del gradiente conjugado]]
- [[Implementacion del gradiente conjugado]]
