# Álgebra Lineal Numérica

Cuadernos de Pluto.jl y apuntes del curso de Álgebra Lineal Numérica (Análisis Numérico - LMA / Métodos Numéricos II - LCD), Facultad de Ingeniería Química, UNL.

## Estructura

- `capitulo1/` a `capitulo4/`: cuadernos de Pluto.jl con la teoría del apunte (transcripta y, donde corresponde, adaptada a Julia) más las resoluciones de los ejercicios.
  - Capítulo 1: introducción a MATLAB/Octave y Julia.
  - Capítulo 2: métodos directos para sistemas lineales.
  - Capítulo 3: factorización QR.
  - Capítulo 4: métodos iterativos estacionarios para sistemas lineales.
- `diapositivas/`: diapositivas de las clases, provistas por la cátedra.
- `libros/`: material bibliográfico.
  - `apunte-numerico.pdf`: apunte de la cátedra (Garau, Morin, Zocola), material de curso.

## Ver los cuadernos online (sin instalar Julia)

Versiones estáticas exportadas a [pluto.land](https://pluto.land/) (las celdas interactivas quedan fijas en su último valor, no son editables):

- [Capítulo 2: Métodos directos para sistemas lineales](https://pluto.land/n/w6upnp9h)
- [Capítulo 3: Factorización QR](https://pluto.land/n/8xinyzyx)

## Cómo abrir los cuadernos

Cada archivo `.jl` dentro de `capituloN/` es un cuaderno de [Pluto.jl](https://plutojl.org/). Para abrirlo:

```julia
using Pluto
Pluto.run()
```

y seleccionar el archivo desde la interfaz de Pluto.
