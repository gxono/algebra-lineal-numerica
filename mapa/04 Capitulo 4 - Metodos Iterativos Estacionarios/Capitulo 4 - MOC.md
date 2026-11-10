---
tags: [moc, capitulo4]
---

# Capítulo 4 - Métodos iterativos estacionarios

> [!idea]
> Para sistemas grandes y sparse, factorizar ($LU$/[[Factorizacion de Cholesky|Cholesky]], Capítulo 2) puede ser demasiado caro (el *fill-in* destruye la estructura sparse). La alternativa: iterar $x_{k+1}=Mx_k+c$ desde algún $x_0$, con $M$ construida a partir de $A$ de forma barata. Converge sii $ρ(M)<1$ ([[Espectro y radio espectral|radio espectral]]), y la velocidad de convergencia la determina exactamente ese radio espectral.

## Teoría

- [[Normas de vectores y matrices]] - $‖·‖_1,2,∞$, norma inducida, submultiplicatividad.
- [[Numero de condicion y error relativo]] - $κ(A)=‖A‖‖A⁻¹‖$; cuán sensible es la solución a perturbaciones en los datos.
- [[Espectro y radio espectral]] - $ρ(A)≤‖A‖$; $ρ(A)=\lim\|A^k\|^{1/k}$.
- [[Metodos iterativos estacionarios - marco general]] - $x_{k+1}=Mx_k+c$; convergencia $\iff\rho(M)<1$.
- [[Metodo de Jacobi]] - $M_J=-D^{-1}(L+U)$.
- [[Metodo de Gauss-Seidel]] - $M_{GS}=-(L+D)^{-1}U$.
- [[Metodo SOR]] - generaliza [[Metodo de Gauss-Seidel|Gauss-Seidel]] con un parámetro de relajación $ω$.
- [[Metodo de Richardson]] - $M_\alpha=I-\alpha A$; el más simple, y el que conecta directo con el [[Capitulo 5 - MOC|Capítulo 5]] (steepest descent = Richardson con $α$ óptimo en cada paso).

## Problemas (4.1 a 4.16)

Ver [[Capitulo 4 - Problemas - Indice]].

## Teoremas y demostraciones

Ver [[Demostraciones - Capitulo 4]] (un archivo por teorema/propiedad, en `Teoremas/`).

## La tabla que resume todo el capítulo

Ver [[Comparacion metodos iterativos estacionarios]].

## Conexiones con otros capítulos

- $κ₂(A)$ generaliza el $κ₂$ rectangular del [[Capitulo 3 - MOC|Capítulo 3]].
- e.d.d. (Capítulo 2) es condición **suficiente** (no necesaria) para convergencia de [[Metodo de Jacobi|Jacobi]]/Gauss-Seidel - ver [[4.11 - Convergencia Jacobi y Gauss-Seidel caso particular|Problema 4.11]].
- Todo el [[Capitulo 5 - MOC|Capítulo 5]] ([[Gradiente conjugado - definicion y propiedades|gradiente conjugado]], [[GMRes - definicion y propiedad de minimizacion|GMRes]]) nace de la pregunta "¿se puede hacer mejor que un método estacionario con $M$ fija?" - la respuesta es sí, cambiando $α$/la dirección en cada paso.
