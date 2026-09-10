using LinearAlgebra


### EJERCICIO 1

#A)
(15/4 - (1 - 2/9) * sqrt(81/16))^3

#B)
(-3 + 2^(-5) * 4) / (cbrt(0.88 - 7/8))

#C)
log(sin(pi/4))^2 - 2*exp(1-3pi)



### EJERCICIO 2
a = 5 * 8
a = a + 5
a = 2 * a
a = a + 10




###EJERCICIO 3
#A) 
A = [1 2 3; 4 5 6; 7 8 9; 10 11 12]

#Esto devuelve lo mismo. Averiguá vos por qué (probá ejecutar desde afuera para adentro).


#B) 
u = -8:5 # Averigua que tipo tiene un vector cargado a mano y uno así con typeof(). Que pasa si usas collect(u)?
v = -8:1.5:5
w = range(-8, 5, 10)


#C)
length(v)
A
size(A)


#D)
w[3]
A[2,3]


#E)
A[3,:]
A[:,3]


#F)
A[[2,4], :]


#G)
A[2:4, :]


#H)
A[[1,4], [2,3]]

A

#I)
A[2,:] = -4 * A[1,:] + A[2,:]


#J)
A[:,[1,3]] = A[:, [3,1]]


#K)
A[2,:] = [] 
#ERROR: DimensionMismatch: tried to assign 0-element array to 1×3 destination


#L) 
u = transpose(u)
#tambien funciona transpose(u) para vectores numericos o permutedims(u) para casos generales.
# Revisa la documentacion de ?permutedims para ver las diferencias entre estos dos.




#M) 
b = v[[2,5,9]]
# Algo como v[2,5,9] no funciona porque implicaria que v es una matriz de R3 y que estas accediendo a la componente 9 de la tercera dimension (cosa que es posible).


#N)
C = [A b']
# DimensionMismatch: number of rows of each array must match (got (4, 1))

#Esto es concatenar horizontalmente vectores o matrices. Lo podes hacer tambien con hcat(X, Y).

C = [A; b']
#Esto es contatenar verticalmente. Tambien funciona vcat(A, b')
#Fijate por que [A; b] daría error.





#EJERCICIO 4
a = rand() #Numero aleatorio entre cero y uno
v = rand(5) #Vector de 5 componentes aleatorias entre cero y 1
A = rand(3,4) #Matriz de numeros aleatorios de 3×4 entre cero y uno.

A = rand(5:10, 2, 3) #matriz de numeros aleatorios enteros entre 5 y 10 de dos filas y tres columnas.
rand(["coso1", "coso2", "coso3"], 2, 5) #matrix de 2×5 de cosos.
#Podes ver mas en ?rand



A = zeros(5)
A = zeros(3,4)
A = zeros(Int, 3, 4)
#Podes ver mas en ?zeros


A = ones(5)
A = ones(3,4)
A = ones(Int, 3, 4)
A = ones(ComplexF64, 3, 4)
#Podes ver mas en ?ones


A = I(1)
A = I(2)
A = I(4)
# Devuelve un objeto del tipo Diagonal. Para obtenerla con ceros y demas usa collect.


# En julia diag hace algo diferente a diag en octave.
# En Octave, diag es de doble propósito: diag(matriz) extrae la diagonal
# y diag(vector) construye una matriz diagonal a partir del vector.
# En Julia, diag SOLO extrae la diagonal de una matriz (no acepta un vector).

M = [1 2 3; 4 5 6; 7 8 9]
diag(M) #Extrae la diagonal principal de M como vector: [1, 5, 9]
diag(M, 1) #Extrae la diagonal ubicada un lugar por encima de la principal: [2, 6]
diag(M, -1) #Extrae la diagonal ubicada un lugar por debajo de la principal: [4, 8]

v = [1, 2, 3]
#diag(v) #Esto tira error: no existe el método diag para un vector.

Diagonal(v) #Construye un objeto de tipo Diagonal a partir de v (eficiente, guarda solo la diagonal).
diagm(v) #Construye una matriz densa (Matrix) con v en la diagonal principal y ceros en el resto.
diagm(1 => v) #Ubica v en la diagonal justo encima de la principal, en una matriz 4x4.
#Y mira
diagm(-1 => [1,2,3], 1 => [4,5,6])




A = rand(4,4)
triu(A)
triu(A,1)
triu(A,2)
triu(A,-1) #:D

A = rand(4,4)
tril(A)
tril(A,1)
tril(A,2)
tril(A,-1)



#B)
u = transpose(rand(-8:8, 6))

v = fill(2, 6)
v = 2 * ones(6)

A = rand(0:10, 4, 6)

B = rand(-10:10, 4, 6)

C = triu(rand(6,6))

D = diagm(v)

E = rand(Int, 6, 6) #Aca depende del tipo de Int que se use. Int8 da enteros entre -255 y 255 por ejemplo. 

F = rand(-7.0:eps():5.0, 6, 8)

I = LinearAlgebra.I(6)


Z = zeros(4,6)




#EJERCICIO 5
A + F

A + B

B + A

A + Z

Z + A


A * Z 

-3A+3B

3(B-A)

(A+B)F

A * F + B * F

C * E

E * C

I * E

E * I

v * u

u * v

u * v'

u .* v'

u ./ v'

u .^(v')

3 .^ v

D^3

A^2

A .^ 2

F + 2

F / 2

F ./ 2

2 / F

2 ./ F

(C + E)^2

C^2 + 2C * E + E^2

log.(A)
exp.(C)
sqrt.(B)
cos.(Z)


inv(C)
inv(E)
det(C)
det(E)