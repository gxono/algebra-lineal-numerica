---
tags: [moc, capitulo2]
---

# Capítulo 2 - Métodos directos para sistemas lineales

> [!idea]
> Idea central del capítulo: resolver $Ax=b$ reduciendo $A$ a forma triangular (Gauss/[[Factorizacion LU|LU]]) o, si $A$ es sdp, a $R^TR$ ([[Factorizacion de Cholesky|Cholesky]]), y después resolver dos sistemas triangulares (baratos, $O(n^2)$). Todo el capítulo es, en el fondo, una sola pregunta: *¿cuándo y cómo se puede triangular $A$ de forma barata y estable?*

## Teoría

- [[Matrices triangulares - sustitucion]] - definición, sustitución hacia adelante/atrás, conteo $O(n^2)$.
- [[Metodo de Gauss - eliminacion]] - eliminación por filas, conteo $O(2n^3/3)$.
- [[Factorizacion LU]] - matrices de Gauss, existencia y unicidad, costo de reusar la factorización para varios lados derechos.
- [[Pivoteo y matrices de permutacion]] - por qué hace falta, pivoteo parcial vs. total, teorema $PA=LU$.
- [[Factorizacion de Cholesky]] - matrices sdp, teorema de existencia (demostración constructiva por inducción), relación con LU.
- [[Costo computacional metodos directos]] - tabla comparativa de todos los costos del capítulo.

## Qué método usar, según $A$

```mermaid
flowchart TD
    A["Resolver Ax=b"] --> B{"¿A es sdp?"}
    B -- Sí --> C["Cholesky: A = RᵀR<br/>costo n³/3"]
    B -- No --> D{"¿A tridiagonal?"}
    D -- Sí --> E["tridilu: LU sin fill-in<br/>costo O(n)"]
    D -- No --> F["Gauss / LU con pivoteo<br/>costo 2n³/3"]
    C --> G["2 sustituciones triangulares<br/>O(n²) cada una"]
    E --> G
    F --> G
    G --> H["x"]
```

La pregunta que organiza todo el capítulo: cuanta más estructura tiene $A$ (simetría, banda), más barato sale triangularla. Ver [[Costo computacional metodos directos]] para los números exactos.

## Cómo se relacionan Gauss, LU y Cholesky

```mermaid
flowchart LR
    G["Eliminación de Gauss<br/>(triangulariza A, paso a paso)"] -- "guardar los multiplicadores m_kj" --> LU["Factorización LU<br/>A = LU"]
    LU -- "A sdp ⟹ L=Rᵀ, U=DR<br/>(explota simetría)" --> CH["Factorización de Cholesky<br/>A = RᵀR, mitad de costo"]
    LU -- "reusar L,U<br/>para cada b nuevo" --> MULTI["k sistemas:<br/>2n³/3 + k·n² en vez de k·2n³/3"]
```

Gauss **es** el algoritmo; LU es Gauss **recordando** lo que hizo (los multiplicadores) para no repetir trabajo con otro $b$; Cholesky es LU **aprovechando** que, si $A$ es simétrica, $L$ y $U$ no son independientes.

## Problemas (2.1 a 2.14)

Ver [[Capitulo 2 - Problemas - Indice]].

## Teoremas y demostraciones

Ver [[Demostraciones - Capitulo 2]] (un archivo por teorema/propiedad, en `Teoremas/`).

## Conexiones con otros capítulos

- `\` de Julia usa estas ideas por detrás - ver [[Resolver sistemas lineales en Julia]] (Capítulo 1).
- Cholesky reaparece como **precondicionador de factorización incompleta** en el Capítulo 5 - ver [[Precondicionadores genericos]].
- El número de condición, apenas mencionado acá, se estudia en profundidad en el [[Capitulo 4 - MOC|Capítulo 4]].
- Las ecuaciones normales del Capítulo 3 ($A^TAx=A^Tb$) se resuelven, precisamente, con Cholesky.
