---
tags: [capitulo5, teorema]
---

# Lema 5.12 (espacios de Krylov)

Si $r_k\neq0$, $\text{span}\{r_0,\dots,r_k\}=\text{span}\{r_0,Ar_0,\dots,A^kr_0\}=:\mathcal{K}_{k+1}$ (el **espacio de Krylov**).

> [!demostracion]- Por inducción en $k$
> **Caso base $k=0$:** $\text{span}\{r_0\}=\mathcal{K}_1$, trivial.
>
> **Paso inductivo:** supongamos $\text{span}\{r_0,\dots,r_k\}=\mathcal{K}_{k+1}$ (H.I.), y veamos que $\text{span}\{r_0,\dots,r_{k+1}\}=\mathcal{K}_{k+2}$ cuando $r_{k+1}\neq0$.
>
> *(⊂)* Por H.I., $r_i\in\mathcal{K}_{k+1}\subset\mathcal{K}_{k+2}$ para $i=0,\dots,k$. Falta ver $r_{k+1}\in\mathcal{K}_{k+2}$: como $x_{k+1}=x_0+\xi$ con $\xi\in\text{span}\{r_0,\dots,r_k\}\overset{\text{H.I.}}{=}\mathcal{K}_{k+1}$,
> $$r_{k+1}=b-Ax_{k+1}=b-Ax_0-A\xi=\underbrace{r_0}_{\in\mathcal K_{k+2}}-\underbrace{A\xi}_{\in\mathcal K_{k+2},\ \text{pues}\ \xi\in\mathcal K_{k+1}}\in\mathcal{K}_{k+2}.$$
> Así, $\text{span}\{r_0,\dots,r_{k+1}\}\subset\mathcal{K}_{k+2}$.
>
> *(=)* Como $r_{k+1}\neq0$, por la [[Gradiente conjugado - definicion y propiedades#Observación 5.11|Observación 5.11]] $\{r_0,\dots,r_{k+1}\}$ es linealmente independiente, así que $\text{span}\{r_0,\dots,r_{k+1}\}$ tiene dimensión $k+2$. Pero $\mathcal{K}_{k+2}$ tiene dimensión a lo sumo $k+2$ (generado por $k+2$ vectores) y acabamos de ver que lo contiene - luego son iguales.

> [!idea]
> Esto es sorprendente: el espacio que generan los residuos (que a priori parecen "cualquier cosa") coincide **exactamente** con el espacio de Krylov $\{r_0,Ar_0,A^2r_0,\dots\}$, que solo depende de $A$ y $r_0$ - no hace falta calcular los residuos intermedios para saber qué espacio están generando. Esto es lo que permite, más adelante, expresar el error como un polinomio evaluado en $A$ (ver [[Cotas de error del gradiente conjugado]]).

## Ver también
- [[Gradiente conjugado - definicion y propiedades]] (nota de tema completa)
- [[Demostraciones - Capitulo 5|Demostraciones]]
