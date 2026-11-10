---
tags: [capitulo5, algoritmo, teorema, definicion]
---

# Precondicionamiento y [[Gradiente conjugado - definicion y propiedades|gradiente conjugado]] precondicionado (GCP)

> [!idea]
> Si $M\approx A^{-1}$, resolver $SASy=Sb$ (con $S=M^{1/2}$) en vez de $Ax=b$ mejora la tasa de convergencia a $\sqrt{c_2/c_1}$ en vez de $\sqrt{\kappa_2(A)}$. La parte elegante: **$S$ nunca se calcula**. Reescribiendo el algoritmo en términos de $x_k=Sy_k$, $S$ siempre aparece de a pares ($S\cdot S=M$), y lo único que hace falta es una rutina $ξ↦Mξ$.

## Buen precondicionador

$M$ sdp, con $c_1\leq\dfrac{x^TMx}{x^TA^{-1}x}\leq c_2$ para todo $x\neq0$. Bueno si $c_2/c_1\approx1$. Caso extremo $M=I$: acción gratis, pero $c_2/c_1=\kappa_2(A)$ (no ayuda). Caso extremo $M=A^{-1}$: $c_2/c_1=1$ (ideal), pero la acción es tan cara como resolver el problema original.

## De $SASy=Sb$ a un algoritmo sin $S$

Aplicando GC a $SASy=Sb$ y sustituyendo $x_k=Sy_k$, $r_k=b-Ax_k$, se obtiene (con $z_k:=Mr_k$, el "residuo precondicionado"):
$$
d_1=Mr_0, \qquad d_{k+1}=Mr_k+\tilde\beta_{k+1}d_k, \qquad \tilde\alpha_{k+1}=\frac{r_k^TMr_k}{d_{k+1}^TAd_{k+1}}, \qquad \tilde\beta_{k+1}=\frac{r_k^TMr_k}{r_{k-1}^TMr_{k-1}}.
$$
$S$ desapareció: solo queda la acción $ξ\mapsto M\xi$.

> [!demostracion]- Cómo desaparece $S$ (la sustitución completa)
> Sea $\{y_k,\tilde r_k,\tilde d_k\}$ la sucesión que genera el GC **sin modificar** ([[Implementacion del gradiente conjugado]]) aplicado al sistema $SASy=Sb$: $\tilde r_0=Sb-SASy_0$, $\tilde d_1=\tilde r_0$, $\tilde d_{k+1}=\tilde r_k+\tilde\beta_{k+1}\tilde d_k$, $\tilde\alpha_{k+1}=\dfrac{\tilde r_k^T\tilde r_k}{\tilde d_{k+1}^T(SAS)\tilde d_{k+1}}$, $\tilde\beta_{k+1}=\dfrac{\tilde r_k^T\tilde r_k}{\tilde r_{k-1}^T\tilde r_{k-1}}$.
>
> Definimos $x_k:=Sy_k$ y $r_k:=b-Ax_k$ (el residuo del problema **original**).
>
> **Los residuos se corresponden vía $S$.** $r_k=b-Ax_k=b-ASy_k=S^{-1}(Sb-SASy_k)=S^{-1}\tilde r_k$, es decir $\tilde r_k=Sr_k$.
>
> **Las direcciones se corresponden vía $S$.** De $y_{k+1}=y_k+\tilde\alpha_{k+1}\tilde d_{k+1}$, multiplicando por $S$: $x_{k+1}=x_k+\tilde\alpha_{k+1}(S\tilde d_{k+1})$, así que definimos $d_{k+1}:=S\tilde d_{k+1}$ (mismo paso $\tilde\alpha_{k+1}$ que en $y$, solo cambia la dirección).
>
> **La recurrencia de las direcciones, en términos de $M=S^2$.** Multiplicando $\tilde d_{k+1}=\tilde r_k+\tilde\beta_{k+1}\tilde d_k$ por $S$, y usando $\tilde r_k=Sr_k$:
> $$d_1 = S\tilde d_1 = S\tilde r_0 = S(Sr_0) = S^2r_0 = Mr_0,$$
> $$d_{k+1} = S\tilde d_{k+1} = S\tilde r_k+\tilde\beta_{k+1}S\tilde d_k = S(Sr_k)+\tilde\beta_{k+1}d_k = Mr_k+\tilde\beta_{k+1}d_k.$$
> Cada vez que aparece $S$, aparece **dos veces** ($S\cdot S=M$) - la dirección de búsqueda ya no es el residuo sino el residuo precondicionado $z_k:=Mr_k$. (Con $M=I$, $S=I$ y esto colapsa exactamente al GC sin precondicionar.)
>
> **Los coeficientes, en términos de $M$.** Usando $\tilde r_k=Sr_k$, $\tilde d_{k+1}=S^{-1}d_{k+1}$, y que $S$ es simétrica ($S$ es sdp):
> $$\tilde\alpha_{k+1} = \frac{\tilde r_k^T\tilde r_k}{\tilde d_{k+1}^TSAS\tilde d_{k+1}} = \frac{(Sr_k)^T(Sr_k)}{(S\tilde d_{k+1})^TA(S\tilde d_{k+1})} = \frac{r_k^TS^TSr_k}{d_{k+1}^TAd_{k+1}} = \frac{r_k^TMr_k}{d_{k+1}^TAd_{k+1}},$$
> $$\tilde\beta_{k+1} = \frac{\tilde r_k^T\tilde r_k}{\tilde r_{k-1}^T\tilde r_{k-1}} = \frac{(Sr_k)^T(Sr_k)}{(Sr_{k-1})^T(Sr_{k-1})} = \frac{r_k^TMr_k}{r_{k-1}^TMr_{k-1}}.$$
> En ambas fórmulas $S$ desapareció de a pares - solo queda la acción $\xi\mapsto M\xi$, nunca $S$ por sí sola.

> [!idea]
> El punto central de toda la deducción: $S=M^{1/2}$ es una herramienta puramente **teórica** para justificar que $SASy=Sb$ es un sistema sdp al que se le puede aplicar GC sin modificar - nunca se calcula en la práctica (calcular una raíz cuadrada de matriz sería mucho más caro que resolver el problema original). El álgebra garantiza que, al traducir todo a la variable $x_k=Sy_k$, $S$ nunca aparece sola, solo como $S\cdot S=M$.

## Lema 5.29 (equivalencia de normas) - necesario para el criterio de parada

$$
\sqrt{c_1}\|v\|_{M^{-1}} \leq \|v\|_A \leq \sqrt{c_2}\|v\|_{M^{-1}}.
$$

Demostración completa (vía Cauchy-Schwarz, $v^TB^{-1}v=\max_{x\neq0}(x^Tv)^2/x^TBx$): [[Lema 5.29]].

Esto permite trasladar la cota de convergencia (naturalmente en norma $M^{-1}$, porque el criterio de parada usa $\|r_k\|_M$) de vuelta a la norma $A$ (en la que está escrita toda la teoría): $\|x_k-x_*\|_A\leq2\big(\frac{\sqrt{c_2/c_1}-1}{\sqrt{c_2/c_1}+1}\big)^k\|x_0-x_*\|_A$.

## Algoritmo

```julia
"""
    gradconjp(A, b, precon, x0; tol=1e-8, max_iter=1000)
"""
function gradconjp(A, b, precon, x0; tol=1e-8, max_iter=1000)
    Av(v) = A isa Function ? A(v) : A*v
    Mv(v) = precon isa Function ? precon(v) : precon*v
    x = copy(x0)
    r = b - Av(x)
    z = Mv(r)
    ρ = r'z
    d = copy(z)
    ρold = 0.0
    iter = 0
    while norm(r) > tol && iter < max_iter
        if iter > 0
            β = ρ/ρold
            d = z + β*d
        end
        w = Av(d)
        α = ρ/(d'w)
        x += α*d
        r -= α*w
        z = Mv(r)
        ρold = ρ
        ρ = r'z
        iter += 1
    end
    return x, iter
end
```

**Costo extra sobre GC:** solo **una** acción $Mξ$ más por iteración (1 vector extra en memoria).

## Ver también
- [[Implementacion del gradiente conjugado]]
- [[Precondicionadores genericos]]
- [[Tablas Comparativas - GC vs GCP]]
