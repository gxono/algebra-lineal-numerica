---
tags: [capitulo4, algoritmo, teorema]
---

# Método de Richardson

> [!idea]
> El método más simple de todos: $x_{k+1}=x_k-\alpha(Ax_k-b)=x_k+\alpha r_k$, con $α$ **fijo** elegido de antemano (no cambia entre iteraciones - a diferencia del [[Capitulo 5 - MOC|descenso más pronunciado]] del Capítulo 5, que es "Richardson con $α$ recalculado óptimamente en cada paso"). $M_\alpha=I-\alpha A$.

## Construcción y convergencia (A sdp)

Si $A$ es sdp con autovalores $0<\lambda_{\min}\leq\cdots\leq\lambda_{\max}$, los autovalores de $M_\alpha=I-\alpha A$ son $1-\alpha\lambda_i$. Converge sii $|1-\alpha\lambda_i|<1$ para todo $i$, es decir, $0<\alpha<2/\lambda_{\max}$.

**$α$ óptimo:** minimizar $\rho(M_\alpha)=\max_i|1-\alpha\lambda_i|=\max\{|1-\alpha\lambda_{\min}|,|1-\alpha\lambda_{\max}|\}$ sobre $α$ da
$$
\alpha_{\text{opt}} = \frac{2}{\lambda_{\min}+\lambda_{\max}}, \qquad \rho(M_{\alpha_{\text{opt}}}) = \frac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\max}+\lambda_{\min}} = \frac{\kappa_2(A)-1}{\kappa_2(A)+1}.
$$
Demostración completa (condición de convergencia + derivación del $\alpha$ óptimo vía el argumento minimax): [[Lema 4.35]]. Con ese $\alpha_{\text{opt}}$, el error decrece con ese mismo factor tanto en norma $2$ ([[Teorema 4.36]]) como en norma $A$ ([[Teorema 4.37]]) - idéntico a la cota del descenso más pronunciado del Capítulo 5 (no es casualidad, ver abajo).

> [!idea]
> Esta última fórmula, $\frac{\kappa_2(A)-1}{\kappa_2(A)+1}$, es **exactamente** la tasa de convergencia del [[Descenso mas pronunciado|descenso más pronunciado]] del Capítulo 5 (ver el método del descenso más pronunciado) - porque steepest descent es Richardson recalculando este $α$ óptimo (para la dirección actual) en cada paso, en vez de fijarlo una sola vez.

## Cuándo no hay ningún $α$ que funcione

Si $A$ tiene algún autovalor $\leq0$ o **autovalores complejos** con parte real de ambos signos, no existe $α$ real que deje a todos los $1-\alpha\lambda_i$ dentro del disco unitario simultáneamente para *todo* $x_0$ - ver [[4.16 - Matriz sin alpha de Richardson que converja|Problema 4.16]].

## Algoritmo

```julia
"""
    richardson(A, b, α; tol=1e-6, maxiter=10_000)
"""
function richardson(A, b, α; tol=1e-6, maxiter=10_000)
    x = zeros(length(b))
    err = Float64[]
    for k in 1:maxiter
        r = b - A*x
        xn = x + α*r
        push!(err, norm(xn - x, 1))
        x = xn
        err[end] < tol && break
    end
    return x, err
end
```

## Ver también
- [[Metodos iterativos estacionarios - marco general]]
- En el Capítulo 5: el [[Capitulo 5 - MOC|método del descenso más pronunciado]] generaliza esto, y el [[Gradiente conjugado - definicion y propiedades|gradiente conjugado]] va más allá todavía.
