---
tags: [glosario, capitulo6, capitulo7]
aliases: [convergencia lineal, convergencia cuadrática, convergencia superlineal]
---

# Orden de convergencia

## Definición

Una sucesión $x_k\to x_*$ converge con **orden $p\geq1$** si existe $C>0$ con
$$
\|x_{k+1}-x_*\|\leq C\|x_k-x_*\|^p\quad\text{(para }k\text{ grande)},
$$
con $C<1$ si $p=1$. Casos: **lineal** ($p=1$, el error se multiplica por $\approx C$ en cada paso), **superlineal** ($\|x_{k+1}-x_*\|/\|x_k-x_*\|\to0$), **cuadrática** ($p=2$, se duplican los dígitos correctos), **cúbica** ($p=3$, se triplican). **Sublineal**: el cociente tiende a $1$ (ej. $x_k\sim1/\sqrt{k}$).

## Cómo estimarlo con datos

Con tres errores consecutivos $e_k=\|x_k-x_*\|$:
$$
p\approx\frac{\log(e_{k+1}/e_k)}{\log(e_k/e_{k-1})}.
$$
Para órdenes altos conviene `BigFloat`: en `Float64` el error llega a $10^{-16}$ en 2–3 pasos y quedan muy pocos cocientes.

## Dónde aparece

| Método | Orden |
|---|---|
| Punto fijo con contracción / métodos estacionarios / potencias / QR sin shift | 1 (lineal) |
| Newton con $F'(x_*)$ invertible | superlineal; 2 si $F'$ es Lipschitz |
| Newton con raíz múltiple ($F'(x_*)$ singular) | 1 |
| Iteración del cociente de Rayleigh / QR con shift de Rayleigh, $A$ simétrica | 3 |

## Ver también
- [[Metodo de Newton - convergencia cuadratica]]
- [[Metodo de las potencias inversas]]
- [[Glosario - Indice|Glosario]]
