---
tags: [capitulo4, teorema]
---

# Teorema: convergencia de un método iterativo estacionario $\iff \rho(M)<1$

> [!demostracion]- La iteración converge para todo $x_0$ $\iff\rho(M)<1$
> Si $x_*$ es la solución (punto fijo, $x_*=Mx_*+c$), el error $e_k:=x_k-x_*$ satisface $e_{k+1}=Mx_k+c-(Mx_*+c)=M(x_k-x_*)=Me_k$, así que $e_k=M^ke_0$. La convergencia **para todo** $x_0$ (equivalentemente, para todo $e_0$) es exactamente la condición (4) del [[Teorema 4.23]] aplicado a $M$: $\lim_k M^ke_0=0\ \forall e_0$. Por ese teorema, esto es equivalente a $\rho(M)<1$ (condición (1)). Si $\rho(M)\geq1$, existe algún $x_0$ para el cual la iteración no converge (ver [[4.9 - No convergencia si rho mayor igual 1|Problema 4.9]]).

**Condición suficiente más simple (pero no necesaria):** si existe alguna norma inducida con $‖M‖<1$, converge (porque $ρ(M)≤‖M‖<1$). Esto es lo que se usa para probar convergencia de Jacobi/Gauss-Seidel en [[2.11 - Matrices e.d.d. son no singulares|matrices e.d.d.]]

## Ver también
- [[Metodos iterativos estacionarios - marco general]] (nota de tema completa)
- [[Demostraciones - Capitulo 4|Demostraciones]]
