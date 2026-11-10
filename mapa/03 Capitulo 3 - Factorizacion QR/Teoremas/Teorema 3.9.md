---
tags: [capitulo3, teorema]
---

# Teorema 3.9: construcción de $u$ para anular un vector

Si $\|x\|_2=1$, tomando $\sigma=\text{sgn}(x_1)$ y
$$
u = \frac{1}{\sqrt{1+\sigma x_1}}(x+\sigma e_1), \qquad H=I-uu^T,
$$
resulta $Hx=-\sigma e_1$. *(El signo se elige para que $1+\sigma x_1$ no sea chico - ver Observación 3.11.)* Si $\|x\|_2\neq1$ se normaliza antes: $Hx=-\sigma\|x\|_2e_1$.

> [!demostracion]- Esquema de la demostración
> $\|u\|_2^2=\frac{1}{1+\sigma x_1}(x^Tx+2\sigma x_1+1)=\frac{2(1+\sigma x_1)}{1+\sigma x_1}=2$. Y $u^Tx=\frac{1+\sigma x_1}{\sqrt{1+\sigma x_1}}=\sqrt{1+\sigma x_1}$, de donde $Hx=x-u(u^Tx)=x-\sqrt{1+\sigma x_1}\,u=x-(x+\sigma e_1)=-\sigma e_1$.

## Ver también
- [[Matrices de Householder]] (nota de tema completa)
- [[Demostraciones - Capitulo 3|Demostraciones]]
