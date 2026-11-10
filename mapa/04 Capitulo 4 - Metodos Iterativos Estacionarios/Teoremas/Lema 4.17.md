---
tags: [capitulo4, teorema]
---

# Lema 4.17

Sea $A\in\mathbb R^{n\times m}$ con $n\geq m$ y rango completo por columnas. Entonces
$$
\kappa_2(A^TA) = \kappa_2(A)^2.
$$

> [!demostracion]-
> Sea $A=U\Sigma V^T$ la SVD de $A$ (valores singulares $\sigma_1\geq\cdots\geq\sigma_m>0$, positivos porque $A$ tiene rango completo por columnas). Entonces
> $$A^TA = V\Sigma^TU^TU\Sigma V^T = V(\Sigma^T\Sigma)V^T, \qquad \Sigma^T\Sigma=\text{diag}(\sigma_1^2,\dots,\sigma_m^2).$$
> Esto es una diagonalización ortogonal de $A^TA$ ($V$ ortogonal), así que los autovalores de $A^TA$ son exactamente $\sigma_1^2,\dots,\sigma_m^2$. Como $A^TA$ es simétrica definida positiva, por el [[Lema 4.14]]:
> $$\kappa_2(A^TA) = \frac{\lambda_{\max}(A^TA)}{\lambda_{\min}(A^TA)} = \frac{\sigma_1^2}{\sigma_m^2} = \Big(\frac{\sigma_1}{\sigma_m}\Big)^2 = \kappa_2(A)^2,$$
> usando en el último paso que $\kappa_2(A)=\sigma_1/\sigma_m$ para matrices rectangulares de rango completo (Capítulo 3, [[Estabilidad y numero de condicion]]).

> [!idea]
> Esta es la razón matemática exacta por la que **nunca** conviene resolver mínimos cuadrados formando las ecuaciones normales $A^TAx=A^Tb$ cuando $A$ está mal condicionada: el condicionamiento del sistema que efectivamente se resuelve es el **cuadrado** del condicionamiento del problema original - y eso es antes de contar el error de redondeo de formar $A^TA$ explícitamente. Ver la comparación completa en [[Estabilidad y numero de condicion]] (Capítulo 3), incluido el ejemplo de la matriz de Läuchli donde $A^TA$ calculada en punto flotante resulta directamente singular.

## Ver también
- [[Lema 4.14]] (usado acá)
- [[Estabilidad y numero de condicion]] (Capítulo 3, nota de tema completa)
- [[Numero de condicion y error relativo]] (nota de tema completa, Capítulo 4)
- [[Demostraciones - Capitulo 4|Demostraciones]]
