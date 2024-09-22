import numpy as np
import scipy.linalg as sc

def calcularLU(A):
    m=A.shape[0] #filas
    n=A.shape[1] #columnas
    Ac= A.copy() 
    Ad =  A.copy() 
    
    if m!=n:
        print('Matriz no cuadrada')
        return 

    ## Aca calculo los valores de los multiplicadores y lo actualizo en A
    for j in range (n): #columnas
        pivote = Ad[j,j]
        for i in range(1, m): #filas
            if j < i:
                k = calculo_k(Ad[i], pivote, j) #calculo k
                if k != 0:
                    Ad[i] = Ad[i] - k*Ad[j]
                    Ac[i][j] = k #agrego k 
                    cant_op += 1
                else:
                    continue
                
    ## Aca actualizo los valores arriba de la diagonal
    for j in range (n): #columnas
        for i in range(1, m): #filas
            if j >= i:
                Ac[i][j] = Ad[i][j]
     
    L = np.tril(Ac,-1) + np.eye(A.shape[0]) 
    U = np.triu(Ac)
    return L, U, cant_op
    
    ###########
    return L, U


def inversaLU(L, U):
    Inv = []
    L, U, cant_op = elim_gaussiana(A)
    y = sc.solve_triangular(L,e, lower=True)
    x = sc.solve_triangular(U,y)
    return Inv


###codigo nari 
def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0] #filas
    n=A.shape[1] #columnas
    Ac= A.copy() 
    Ad =  A.copy() 
    
    if m!=n:
        print('Matriz no cuadrada')
        return 

    ## Aca calculo los valores de los multiplicadores y lo actualizo en A
    for j in range (n): #columnas
        pivote = Ad[j,j]
        for i in range(1, m): #filas
            if j < i:
                k = calculo_k(Ad[i], pivote, j) #calculo k
                if k != 0:
                    Ad[i] = Ad[i] - k*Ad[j]
                    Ac[i][j] = k #agrego k 
                    cant_op += 1
                else:
                    continue
                
    ## Aca actualizo los valores arriba de la diagonal
    for j in range (n): #columnas
        for i in range(1, m): #filas
            if j >= i:
                Ac[i][j] = Ad[i][j]
     
    L = np.tril(Ac,-1) + np.eye(A.shape[0]) 
    U = np.triu(Ac)
    return L, U, cant_op
    

def calculo_k(fila_actual, divisor, iterador):
    multiplicador = 0
    if divisor != 0:
        multiplicador = fila_actual[iterador] / divisor
    return multiplicador


#### habia armado esta de permutar filas pero no llegue a hacer que ande bien

def permutarFilas(matriz, i):
    while matriz[i][i] == 0:
        if i>= matriz.shape[0]-1:
            return "no se puede hacer la LU"
        else:
            matriz[[i, i + 1]] = matriz[[i + 1, i]]
    return matriz
   


def coefTec(z,p):
    p = np.diag(p.values) #Diagonalizo el vector
    per,l,u = sc.lu(p) #Lu de p
    inv_p = sc.inv(per@l@u) #Inversa de p
    return z@inv_p

def Leont2Reg(A,P): #Funcion de Leontief para 2 regiones, usando la formula (I-A)P = D
    m=A.shape[0] #filas A
    Id = np.identity(m)
    _A = Id - A
    return _A @ P
