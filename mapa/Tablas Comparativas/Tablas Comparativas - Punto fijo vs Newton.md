---
tags: [comparacion, capitulo6]
---

# Punto fijo vs. Newton (vs. Newton simplificado vs. minimizar sin derivadas)

> [!comparacion]
> | | Punto fijo $x=G(x)$ | Newton | Newton simplificado | Nelder–Mead sobre $\Vert F\Vert ^2$ |
> |---|---|---|---|---|
> | Qué necesita | una reformulación $x=G(x)$ (no única) | $F$ y $F'$ | $F$ y $F'(x_0)$ una vez | solo $\Vert F\Vert ^2$ |
> | Costo por paso | evaluar $G$ | evaluar $F$, $F'$ + resolver $n\times n$ ($O(n^3)$ denso) | evaluar $F$ + 2 triangulares ($O(n^2)$) | varias evaluaciones de $h$ |
> | Orden | lineal, factor $\gamma$ o $\rho(G'(x_*))$ | **cuadrático** (si $F'$ Lipschitz, $F'(x_*)$ invertible) | lineal, factor chico si $x_0\approx x_*$ | lineal o peor |
> | Garantía de existencia/unicidad | **sí**, si es contracción en $D$ (Teo. 6.12) | no | no | no |
> | Convergencia | global en $D$ (contracción) o local (Ostrowski) | **local** (global solo con convexidad: Newton-Baluev) | local | depende del símplex inicial |
> | Precisión alcanzable | máquina | máquina | máquina | $\sim\sqrt{\varepsilon_{\text{maq}}}$ ($h$ es cuadrática en el error) |
> | Ejemplo del vault | [[6.4 - Region para el teorema de contraccion\|6.4]]: tasa $\rho(G'(x_*))\approx0.32$, 24 pasos a $10^{-12}$ | [[6.10 - Orden experimental de Newton\|6.10]]: 4 pasos a $10^{-16}$ | [[6.16 - Newton simplificado\|6.16]]: tasa $0.02$, 8 pasos | [[6.11 - Cuatro ceros y fminsearch\|6.11]]: ~120 evaluaciones para $10^{-9}$ |

## Qué falla y cómo se nota

| Síntoma | Causa probable | Ejemplo |
|---|---|---|
| Newton lineal en vez de cuadrático | $F'(x_*)$ singular (raíz múltiple) | [[6.14 - Newton para x al cuadrado\|6.14]] |
| Newton superlineal pero no cuadrático | $F'$ no Lipschitz en $x_*$ | [[6.15 - Newton para x mas x a la 1 mas alfa\|6.15]] |
| Newton más que cuadrático | $F''(x_*)=0$ | [[6.17 - Newton en problema discretizado\|6.17]] (orden 3) |
| Newton oscila / diverge | $x_0$ lejos y sin convexidad | [[6.18 - Newton-Baluev sin convexidad\|6.18]] |
| Newton cae en un ciclo | dinámica de Newton cerca de puntos con $F'\approx0$ | [[6.13 - Minimo con Newton\|6.13]] (período 5) |
| Punto fijo no converge con $\rho(G'(x_*))=1$ | Ostrowski no decide | [[6.1 - Rho igual a 1 no decide\|6.1]] |
| Iteración se va a infinito pese a "acercar" puntos | contracción no uniforme, sin punto fijo | [[6.3 - Contraccion en acotados\|6.3]] |

## Estrategia práctica

1. Graficar curvas de nivel cero ($n=2$) para ubicar raíces y elegir $x_0$.
2. Si hay que **garantizar** existencia y unicidad en una región: punto fijo + Corolario 6.19 + Teorema 6.12.
3. Para **calcular** con precisión: Newton desde un $x_0$ cercano (o desde la salida de unas iteraciones de punto fijo).
4. Si $F'$ es cara y $n$ grande: Newton simplificado (factorizar una vez) o reactualizar $F'$ cada algunos pasos.

## Ver también
- [[Capitulo 6 - MOC]]
- [[Tablas Comparativas - Indice]]
