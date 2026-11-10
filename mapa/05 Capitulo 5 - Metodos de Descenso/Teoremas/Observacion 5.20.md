---
tags: [capitulo5, teorema]
---

# Observación 5.20 (cuántas iteraciones hacen falta)

Para $\|x_k-x_*\|_A\leq\varepsilon\|x_0-x_*\|_A$ alcanza con
$$
k \geq \frac{\sqrt{\kappa_2(A)}+1}{2}\ln\frac2\varepsilon = O\big(\sqrt{\kappa_2(A)}\big).
$$

> [!demostracion]-
> Sea $q:=\dfrac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}\in[0,1)$, de modo que la cota del [[Teorema 5.19]] es $\|x_k-x_*\|_A\leq2q^k\|x_0-x_*\|_A$. Alcanza con pedir $2q^k\leq\varepsilon$, es decir (tomando $\ln$, con $q<1$ así que $\ln q<0$):
> $$k\ln\frac1q \geq \ln\frac2\varepsilon \iff k\geq\frac{\ln(2/\varepsilon)}{\ln(1/q)}.$$
> Usamos la desigualdad elemental $\ln(y)\geq1-\dfrac1y$ para $y>0$ (la recta tangente a $\ln$ en $y=1$ queda por debajo de $\ln$, por concavidad), con $y=1/q$:
> $$\ln\frac1q \geq 1-q = 1-\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1} = \frac{2}{\sqrt{\kappa_2(A)}+1}.$$
> Como esta es una cota **inferior** de $\ln(1/q)$, reemplazarla en el denominador da una cota **superior** (más conservadora, pero válida) para el $k$ necesario:
> $$k \geq \frac{\ln(2/\varepsilon)}{2/(\sqrt{\kappa_2(A)}+1)} = \frac{\sqrt{\kappa_2(A)}+1}{2}\ln\frac2\varepsilon.$$

La dependencia en la tolerancia es solo logarítmica: pedir una cifra más de precisión cuesta un número *fijo* de iteraciones; duplicar $\sqrt\kappa$ duplica el total.

## Ver también
- [[Cotas de error del gradiente conjugado]] (nota de tema completa)
- [[Teorema 5.19]]
- [[Demostraciones - Capitulo 5|Demostraciones]]
