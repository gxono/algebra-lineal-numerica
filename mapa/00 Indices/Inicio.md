---
tags: [moc]
---

# Álgebra Lineal Numérica - Mapa de repaso para el parcial

Bóveda de repaso para el parcial de Análisis Numérico (LMA) / Métodos Numéricos II (LCD), basada en el apunte de la cátedra (`apunte-numerico.pdf`), las diapositivas de todas las clases, y mis propios cuadernos de Pluto (teoría + actividades resueltas), capítulos 1 a 7. El **parcial 2** abarca hasta el Capítulo 7 inclusive.

## Cómo está organizado esto

- **Una nota por concepto** (teorema, definición, algoritmo, observación importante), enlazada con `[[wikilinks]]`. El objetivo es poder navegar saltando de un concepto a otro, no leer todo de corrido.
- Cada nota de **Problema/Actividad** tiene la resolución completa. Las que ya había resuelto yo misma están transcriptas tal cual (a veces prolijizadas); las que estaban pendientes las resolví para esta bóveda.
- Todo el código está en **Julia** (no Octave), embebido y comentado.
- Las comparaciones entre métodos/algoritmos viven en dos lugares: dentro de la nota del tema correspondiente, y centralizadas en [[Tablas Comparativas - Indice|Tablas Comparativas]].

## Cómo navegar y buscar esto en Obsidian

No hace falta ningún plugin - todo lo de acá abajo es Obsidian nativo.

**"Sé el nombre, quiero ir directo"** → `Ctrl/Cmd+O` (quick switcher): escribí parte del nombre y Enter. Funciona para notas de teoría ("multigrid"), problemas ("5.12") o los índices ("formulario", "demostraciones").

**"No sé el nombre, busco por tema"** → `Ctrl/Cmd+Shift+F` (búsqueda global, busca en el *contenido* de todas las notas, no solo títulos). Operadores útiles:
- `tag:#teorema` → todas las notas con ese tag.
- `tag:#teorema cholesky` → notas con el tag Y la palabra.
- `path:"Capitulo 5"` → restringe a una carpeta.
- `"frase exacta entre comillas"` → búsqueda literal.

**"Quiero ver todo lo de un tipo (todos los teoremas, todos los problemas)"** → panel de **tags** (ícono de etiqueta en la barra lateral izquierda). Lista cada tag con su cantidad de notas; click para filtrar. Los tags están en [[#Convenciones|Convenciones]] más abajo.

**"Quiero una demostración puntual sin andar buscando en qué nota vive"** → mejor que el tag `#teorema` (que es por archivo, no por teorema individual): andá directo a [[Demostraciones - Indice|Demostraciones]] del capítulo que corresponda y buscá el número/nombre con `Ctrl/Cmd+F` dentro de esa nota.

**"¿Qué otras notas mencionan esto?"** → al abrir cualquier nota, el panel derecho (o abajo del todo en vista móvil) tiene **Linked mentions / Backlinks**: todo lo que apunta a esta nota. Arriba de eso, **Outgoing links** muestra a dónde apunta esta nota.

**"Quiero ver el mapa completo de un vistazo"** → ícono de **grafo** en la barra lateral izquierda (coloreado por capítulo y tipo), o [[Canvas/Mapa del curso.canvas|el Canvas armado a mano]] con las conexiones entre capítulos (más curado que el grafo automático, mejor para entender el panorama general).

**"No entiendo este término"** → [[Glosario - Indice|Glosario]]: definiciones cortas de la jerga que se usa sin reexplicar en cada nota (sparse, fill-in, precondicionador, etc.).

**Para tener siempre a mano los 5 índices** (este archivo, [[Demostraciones - Indice|Demostraciones]], [[Formulario - Indice|Formulario]], [[Glosario - Indice|Glosario]], [[Pendientes - repasar antes del parcial|Pendientes]]): clic derecho sobre cada uno → "Bookmark", quedan fijos en el panel de **Bookmarks** (ícono de marcador en la barra lateral) sin tener que buscarlos cada vez.

## Convenciones

- Tags por capítulo: `#capitulo1` … `#capitulo7`.
- Tags por tipo: `#teorema` `#lema` `#corolario` `#definicion` `#observacion` `#algoritmo` `#problema` `#comparacion`.
- Callouts: `> [!demostracion]` (demostración completa, plegable), `> [!idea]` (la idea clave en una línea), `> [!algoritmo]` (pseudocódigo o código), `> [!comparacion]` (tabla comparativa), `> [!pendiente]` (algo que vale la pena repasar con más cuidado antes del parcial).
- Cada nota de teoría tiene, arriba de todo, un bloque `> [!idea]` con el resumen de 2-3 líneas para repaso rápido el día del examen, y la demostración completa más abajo dentro de `> [!demostracion]-` (el `-` la deja plegada por default).

## Capítulos

- [[Capitulo 1 - MOC|Capítulo 1 - Introducción a Julia/Octave]] *(tratamiento liviano)*
- [[Capitulo 2 - MOC|Capítulo 2 - Métodos directos para sistemas lineales]]
- [[Capitulo 3 - MOC|Capítulo 3 - Factorización QR]]
- [[Capitulo 4 - MOC|Capítulo 4 - Métodos iterativos estacionarios]]
- [[Capitulo 5 - MOC|Capítulo 5 - Métodos de descenso (GC, GMRes, multigrid)]]
- [[Capitulo 6 - MOC|Capítulo 6 - Sistemas de ecuaciones no lineales (punto fijo, Newton)]]
- [[Capitulo 7 - MOC|Capítulo 7 - Problemas de autovalores (potencias, potencias inversas, QR)]]

## Para repasar de punta a punta

- [[Tablas Comparativas - Indice|Tablas comparativas entre todos los métodos]] - incluye, al final, una guía rápida de "qué método usar según el problema".
- [[Formulario - Indice|Formulario]] - solo fórmulas, sin demostraciones, para consulta rápida el día del examen.
- [[Demostraciones - Indice|Demostraciones]] - todos los teoremas/lemas/proposiciones/corolarios con su prueba completa, un archivo por capítulo, para practicar "¿puedo reproducir esto?".
- [[Glosario - Indice|Glosario]] - jerga recurrente (sparse, fill-in, precondicionador, número de condición...) explicada una sola vez, en un solo lugar.
- [[Pendientes - repasar antes del parcial|Pendientes]] - todos los errores típicos y sutilezas marcadas en el vault, juntas en un solo lugar.
- [[Parcial Modelo - Indice|Parcial Modelo]] - un examen modelo completo, enunciado y resolución, con links a cada concepto que hizo falta.
- [[Referencia de Julia]] - qué función llamar para cada cuenta (normas, autovalores, factorizaciones...), built-in vs. implementación propia, y lo mínimo de Symbolics.jl.

## Fuentes

- `apunte-numerico.pdf` (Garau, Morin, Zocola), capítulos 1 a 7.
- Diapositivas de las clases 1 a 14 (`diapositivas/`), que llegan hasta el Capítulo 5.
- Mis cuadernos de Pluto: `capitulo1/` a `capitulo6/`. El Capítulo 7 no tiene cuaderno: sus notas salen directamente del apunte.
- Los problemas de los Capítulos 6 y 7 están todos resueltos acá (del cuaderno del Cap. 6 solo venían resueltos el 6.9 y el 6.10), con cálculos y gráficos verificados ejecutando Julia.
