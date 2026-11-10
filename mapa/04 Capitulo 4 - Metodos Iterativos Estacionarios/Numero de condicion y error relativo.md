---
tags: [capitulo4, definicion, teorema]
---

# Número de condición y error relativo

> [!idea]
> $κ(A)=‖A‖‖A⁻¹‖≥1$ mide cuánto puede amplificarse, en el peor caso, un error relativo en los datos ($A$ o $b$) al pasar a la solución $x$. Un sistema con $κ(A)$ grande es *intrínsecamente* sensible: ningún algoritmo, por más preciso que sea, puede devolver una solución mucho mejor que lo que permite $κ(A)$.

## El resultado central (Lema 4.10, ya usado en capítulos 3 y 5)

Si $\tilde x$ resuelve $A\tilde x=\tilde b$ (aproximadamente, $\tilde b\approx b$) y $x$ resuelve $Ax=b$:
$$
\frac{1}{\kappa(A)}\frac{\|\tilde b-b\|}{\|b\|} \leq \frac{\|\tilde x-x\|}{\|x\|} \leq \kappa(A)\frac{\|\tilde b-b\|}{\|b\|}.
$$

Demostración completa: [[Lema 4.10]].

**Interpretación:** si el residuo relativo $\|r_k\|/\|b\|$ es chico, el error relativo $\|x_k-x_*\|/\|x_*\|$ puede seguir siendo grande si $κ(A)$ es grande - *residuo chico no implica error chico* si la matriz está mal condicionada. Este es el motivo por el que todas las cotas de convergencia del Capítulo 5 (GC, [[GMRes - definicion y propiedad de minimizacion|GMRes]]) terminan con un factor $κ(A)$ o $κ_2(A)$ multiplicando.

## $κ₂(A)$ para $A$ simétrica

$\kappa_2(A)=\lambda_{\max}/\lambda_{\min}$ (en valor absoluto). Siempre $κ(A)≥1$, con igualdad para $A$ ortogonal (ver [[Matrices ortogonales]]): por eso los métodos basados en matrices ortogonales son estables. Demostración completa: [[Lema 4.14]].

## $κ$ para matrices rectangulares

Generaliza el $κ₂(A)=σ₁/σₘ$ del Capítulo 3 (ver [[Estabilidad y numero de condicion]]): mide la sensibilidad del problema de [[Minimos cuadrados y ecuaciones normales|mínimos cuadrados]] $\min\|Ax-b\|_2$. Además $\kappa_2(A^TA)=\kappa_2(A)^2$ para $A$ de rango completo por columnas - demostración completa: [[Lema 4.17]].

## Ver también
- [[Espectro y radio espectral]]
- [[Estabilidad y numero de condicion]] (Capítulo 3)
- [[Capitulo 5 - MOC|Capítulo 5]]: todas las cotas de convergencia de GC/GMRes dependen de $κ₂(A)$ - el motivo central del precondicionamiento.
