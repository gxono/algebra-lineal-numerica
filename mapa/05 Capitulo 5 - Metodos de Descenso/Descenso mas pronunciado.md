---
tags: [capitulo5, algoritmo, teorema]
---

# Método del descenso más pronunciado (Steepest Descent)

> [!idea]
> Moverse en la dirección $-\nabla\phi(x_k)=r_k=b-Ax_k$ (el residuo) con el $α$ óptimo en cada paso. Es [[Metodo de Richardson|Richardson]] (Capítulo 4) con $α$ recalculado cada vez. Converge, pero puede *zigzaguear* mucho si $A$ está mal condicionada - eso es exactamente lo que el [[Gradiente conjugado - definicion y propiedades|gradiente conjugado]] arregla.

## La iteración

$$
x_{k+1}=x_k+\alpha r_k, \qquad \alpha=\frac{r_k^Tr_k}{r_k^TAr_k} \quad(\text{minimiza }\phi(x_k+\alpha r_k)\text{ exactamente}).
$$

```julia
function descmaspron(A, b, x0; tol=1e-8, max_iter=1000)
    x = copy(x0)
    r = b - A*x
    while norm(r) > tol && max_iter > 0
        Ar = A*r
        α = (r'r) / (r'Ar)
        x += α*r
        r -= α*Ar          # = b - A*x, sin volver a multiplicar por A
        max_iter -= 1
    end
    return x
end
```

## Por qué zigzaguea (Ejemplo 5.8)

Las curvas de nivel de $\phi$ son **elipses** centradas en $x_*$, con semieje mayor/menor $1/\sqrt{\lambda_{\min}}$, $1/\sqrt{\lambda_{\max}}$. El residuo $r_k$ es ortogonal (euclídeo) a la curva de nivel en $x_k$, pero **no apunta hacia $x_*$** salvo que la elipse sea una esfera ($A=cI$) o $x_k-x_*$ sea autovector de $A$. Cuanto más alargada la elipse (mayor $κ_2(A)$), peor el desvío, y el método avanza en zigzag, deshaciendo parte del paso anterior en cada iteración.

**Dos residuos consecutivos son siempre ortogonales:** $r_{k+1}^Tr_k=0$ (eq. 5.2) - consecuencia directa de que $\alpha$ se eligió para minimizar exactamente en esa dirección.

**El peor caso se alcanza exactamente**, con $A=\text{diag}(\lambda_1,\lambda_2)$ y $x_0-x_*=(1/\lambda_1,1/\lambda_2)$.

## Teorema 5.9 (cota de convergencia)

$$
\|x_{k+1}-x_*\|_A \leq \frac{\kappa_2(A)-1}{\kappa_2(A)+1}\|x_k-x_*\|_A, \qquad k\geq0.
$$

Misma cota que Richardson con $α$ óptimo global - la "ventaja" de recalcular $α$ en cada paso no mejora la cota en el peor caso, pero sí el comportamiento típico (SD es, en cada paso, **al menos tan bueno** como Richardson con el mejor $\alpha$ fijo posible, porque SD elige el óptimo para la dirección actual). Demostración completa: [[Teorema 5.9]].

## Ver también
- [[Metodo de Richardson]] (Capítulo 4)
- [[Gradiente conjugado - definicion y propiedades]] (la mejora real: no solo cambiar $α$, cambiar la **dirección**)
- [[Tablas Comparativas - GC vs SD|GC vs. SD]]
