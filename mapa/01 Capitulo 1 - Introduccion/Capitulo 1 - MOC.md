---
tags: [moc, capitulo1]
---

# Capítulo 1 - Introducción a Julia/Octave

> [!idea]
> Este capítulo es solo manejo de herramienta (sintaxis de Julia, matrices, gráficos, control de flujo). No suele pesar en el parcial como contenido teórico propio, pero todo lo que se usa acá (indexado, slicing, `\` para resolver sistemas) se usa todo el curso. Tratamiento liviano a propósito.

## Temas

- [[Sintaxis basica de Julia]] - expresiones, matrices, indexado, operaciones elemento a elemento.
- [[Graficos con CairoMakie]] - `lines`, `lines!`, `scatter!`, curvas paramétricas y 3D.
- [[Resolver sistemas lineales en Julia]] - el operador `\`, primer contacto con matrices tridiagonales.
- [[Sentencias de control en Julia]] - `for`, `while`, funciones, recursión.

## Por qué importa para el resto del curso

- El operador `\` (backslash) es el punto de referencia contra el que se comparan *todos* los métodos de los capítulos 2 a 5 (ver [[Resolver sistemas lineales en Julia]]).
- `tril`/`triu`/`diag`/`Diagonal`/`Tridiagonal`/`Bidiagonal` aparecen constantemente desde el Capítulo 2 en adelante.
- `@benchmark` (de `BenchmarkTools.jl`) es la herramienta que se usa en los capítulos 2 y 3 para comparar el costo real de las implementaciones.
