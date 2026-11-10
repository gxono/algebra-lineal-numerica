---
tags: [capitulo5, teorema, definicion]
---

# Gradiente conjugado: definición y propiedades básicas

> [!idea]
> En vez de minimizar $\phi$ solo sobre la recta $x_k+\alpha r_k$ (steepest descent), minimizar sobre **todo** el afín $x_0+\mathcal{K}_{k+1}$ generado por todos los residuos vistos hasta ahora. Esto da una sucesión que converge mucho más rápido, sin "olvidarse" de las direcciones anteriores.

## La definición

$$
x_{k+1} := \operatorname*{argmin}_{x\in x_0+\mathcal{K}_{k+1}}\phi(x), \qquad \mathcal{K}_{k+1}:=\text{span}\{r_0,r_1,\dots,r_k\}, \quad r_i=b-Ax_i.
$$

## Lema 5.10 (residuos ortogonales)

$r_k\perp\{r_0,\dots,r_{k-1}\}$ para todo $k\geq1$. Demostración completa: [[Lema 5.10]].

### Observación 5.11

Si $r_k\neq0$, $\{r_0,\dots,r_k\}$ es linealmente independiente (vectores no nulos, ortogonales dos a dos).

## Lema 5.12 (espacios de Krylov)

Si $r_k\neq0$, $\text{span}\{r_0,\dots,r_k\}=\text{span}\{r_0,Ar_0,\dots,A^kr_0\}=:\mathcal{K}_{k+1}$ (el **espacio de Krylov**). Demostración completa (inducción en $k$): [[Lema 5.12]].

> [!idea]
> Esto es sorprendente: el espacio que generan los residuos (que a priori parecen "cualquier cosa") coincide **exactamente** con el espacio de Krylov $\{r_0,Ar_0,A^2r_0,\dots\}$, que solo depende de $A$ y $r_0$ - no hace falta calcular los residuos intermedios para saber qué espacio están generando. Esto es lo que permite, más adelante, expresar el error como un polinomio evaluado en $A$ (ver [[Cotas de error del gradiente conjugado]]).

## Ver también
- [[Funcional de energia y caracterizacion]]
- [[Cotas de error del gradiente conjugado]]
- [[Implementacion del gradiente conjugado]]
