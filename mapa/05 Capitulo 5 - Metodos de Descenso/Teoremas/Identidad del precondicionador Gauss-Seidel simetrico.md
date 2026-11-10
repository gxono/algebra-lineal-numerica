---
tags: [capitulo5, teorema]
---

# El precondicionador de Gauss-Seidel simétrico coincide con un barrido de ida y vuelta

Sea $A$ simétrica, $A=L+D+L^T$ ($L$ triangular inferior con diagonal nula, $D$ diagonal). Sean
$$
B_{GS}:=-(L+D)^{-1}L^T \ \text{(Gauss-Seidel hacia adelante)}, \qquad B_{BGS}:=-(L^T+D)^{-1}L \ \text{(Gauss-Seidel hacia atrás)},
$$
$$
B_{SGS}:=B_{BGS}B_{GS}, \qquad M_{SGS}:=(L^T+D)^{-1}D(L+D)^{-1} \ \text{(el precondicionador usado en [[Precondicionadores genericos]])}.
$$
Entonces
$$
I-M_{SGS}A = B_{SGS}.
$$

> [!demostracion]-
> Usando $A=L+D+L^T$:
> $$I-M_{SGS}A = I-(L^T+D)^{-1}D(L+D)^{-1}(L+D+L^T) = I-(L^T+D)^{-1}D\underbrace{(L+D)^{-1}(L+D)}_{=I}-(L^T+D)^{-1}D(L+D)^{-1}L^T.$$
> El término del medio se simplifica a $(L^T+D)^{-1}D$, así que
> $$I-M_{SGS}A = \big[I-(L^T+D)^{-1}D\big] - (L^T+D)^{-1}D(L+D)^{-1}L^T. \tag{i}$$
> **Primer corchete.** Como $D=(L^T+D)-L^T$:
> $$(L^T+D)^{-1}D = (L^T+D)^{-1}\big[(L^T+D)-L^T\big] = I-(L^T+D)^{-1}L^T \quad\Longrightarrow\quad I-(L^T+D)^{-1}D = (L^T+D)^{-1}L^T.$$
> Sustituyendo en (i):
> $$I-M_{SGS}A = (L^T+D)^{-1}L^T - (L^T+D)^{-1}D(L+D)^{-1}L^T = (L^T+D)^{-1}\Big[I-D(L+D)^{-1}\Big]L^T. \tag{ii}$$
> **Corchete de (ii).** Como $L=(L+D)-D$:
> $$I-D(L+D)^{-1} = (L+D)(L+D)^{-1}-D(L+D)^{-1} = \big[(L+D)-D\big](L+D)^{-1} = L(L+D)^{-1}.$$
> Sustituyendo en (ii):
> $$I-M_{SGS}A = (L^T+D)^{-1}L(L+D)^{-1}L^T = B_{BGS}\,B_{GS} = B_{SGS}.$$

> [!idea]
> Esto confirma que aplicar $M_{SGS}$ como precondicionador es **exactamente** equivalente, en términos de la matriz de iteración resultante, a encadenar un barrido de Gauss-Seidel hacia adelante seguido de uno hacia atrás - no es solo una fórmula que "se parece" a Gauss-Seidel, es *literalmente* el método estacionario correspondiente a dos medios-pasos de GS. Esto es lo que garantiza que $M_{SGS}$ sea simétrica y, bajo las hipótesis usuales ($A$ sdp), definida positiva - condición necesaria para que tenga sentido usarla como precondicionador del GC (ver [[Precondicionamiento y GCP]]).

## Ver también
- [[Precondicionadores genericos]] (nota de tema completa, donde se define $M_{SGS}$ y se usa en el experimento numérico)
- [[Metodo de Gauss-Seidel]] (Capítulo 4)
- [[Demostraciones - Capitulo 5|Demostraciones]]
