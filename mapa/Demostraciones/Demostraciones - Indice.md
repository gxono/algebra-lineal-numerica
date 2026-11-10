---
tags: [moc, demostracion]
---

# Demostraciones - Índice

Todos los teoremas, lemas, proposiciones y corolarios del curso con su demostración completa, un archivo por capítulo - para repasar "¿puedo reproducir esta prueba?" de punta a punta, sin tener que saltar entre la nota de teoría, el algoritmo y los diagramas. El Capítulo 1 no tiene (es solo sintaxis de Julia).

- [[Demostraciones - Capitulo 2]] - 21 demostraciones (triangulares, LU, pivoteo, Cholesky, e.d.d., matrices sdp).
- [[Demostraciones - Capitulo 3]] - 3 demostraciones (Gram-Schmidt, Householder).
- [[Demostraciones - Capitulo 4]] - 36 demostraciones (Cauchy-Schwarz, normas, Schur, radio espectral (incl. fórmula de Gelfand), número de condición, convergencia de Jacobi/Gauss-Seidel/SOR/Richardson, cota cuantitativa de Gauss-Seidel).
- [[Demostraciones - Capitulo 5]] - 33 demostraciones (producto escalar $A$, funcional de energía, GC, implementación de GC (recurrencia de dos términos), cotas de Chebyshev, precondicionador Gauss-Seidel simétrico, multigrid, operador de Galerkin, GMRes, breakdown de Arnoldi).
- [[Demostraciones - Capitulo 6]] - 17 demostraciones (Ostrowski, contracciones, TFC vectorial y desigualdad del valor medio, lemas de perturbación, Newton local/cuadrático/simplificado, $G'(x_*)=0$, Newton-Baluev con unicidad) + 7 problemas que son demostraciones.
- [[Demostraciones - Capitulo 7]] - 16 resultados (matriz compañera, hermíticas, potencias, cociente de Rayleigh, potencias inversas, $A^{k+1}=P_kU_k$, esquinas de QR, shifts, autopares de $A_N$, unicidad de QR, Hessenberg).

> [!pendiente]
> Los resultados del curso que **no** se demuestran son dos. El Teorema 5.34 (convergencia del ciclo en V de multigrid): el propio apunte dice explícitamente que su demostración excede los alcances del curso. Lo que sí hay que saber reproducir es la cadena de deducciones que va *desde* ese teorema (tomado como dato) *hasta* la cota de $\kappa_2$ del sistema precondicionado, y esa parte está completa en [[Demostraciones - Capitulo 5]]. Y el [[Teorema 7.11]] (convergencia del algoritmo QR), que remite a Stoer–Bulirsch; lo que sí se demuestra es **por qué** funciona (esquinas = potencias / potencias inversas, [[Observacion 7.12]]).

## Ver también
- [[Inicio]]
- [[Formulario - Indice]]
- [[Glosario - Indice|Glosario]]
- [[Pendientes - repasar antes del parcial]]
