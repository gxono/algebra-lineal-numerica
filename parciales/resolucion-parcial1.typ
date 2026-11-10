#set text(font: "Fira Sans")


+ Dada $M$
  + sabemos por teorema que $||M||_1$ es el máximo de la suma por columnas, $||M||_infinity$ es el máximo de la suma por filas y
    $ ||M||_2 = max{sqrt(lambda) : lambda in sigma(A^T A)} $
    Es fácil ver que $||M||_1 = 1.8$, $||M||_infinity=1.8$ mientras que para calcular $||M||_2$ utilizo el operador `opnorm` de Julia y obtengo `1.6180160461808182` aproximadamente.
  
  + La iteracion $k+1$ viene dada por
    $ x_(k+1) = M x_k = M^k x_0 $
    Es decir, 
    $ lim_(k arrow infinity) x_(k+1) = lim_(k arrow infinity) M^k x_0 $
    Como $0 <= rho(M)=0.8<1$, tenemos que $lim_(k arrow infinity) M^k = 0$ por lo que el límite vale cero, independientemente del valor de $x_0$.

+ 
