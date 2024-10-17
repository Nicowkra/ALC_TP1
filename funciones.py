"""
Materia: Algebra Lineal Computacional - FCEyN - UBA
Motivo  : 1er Trabajo Practico
Autor  : Nicolas, Valentin Carcamo, Nadina Soler
"""

# =============================================================================
# IMPORTS
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg as sc
# =============================================================================
# FUNCIONES PARA CALCULAR LU E INVERSA DE UNA MATRIZ
# =============================================================================
#Aca intercambio las filas, entra una matriz y las dos filas a intercambiar
def intercambiarfilas(A, fila1, fila2):
    A[[fila1, fila2]] = A[[fila2, fila1]]
    return A
#En calcularLU entra una matriz A y devuelve su factorizacion LU y un verctor de permutacion
def calcularLU(A):
    m, n = A.shape
    '''Si no es matriz cuadrada, no es invertible, 
    entonces no podemos calcular la factorización LU'''
    if m != n:
        print('Matriz no cuadrada')
        return
    '''
    Iniciamos el vector de permutaciones
    '''
    P = np.arange(n) 
    Ac = A.copy()
    '''
    Recorremos las filas de la matriz A y si el pivote es cero, intercambiamos
    la fila con la siguiente
    '''  
    for fila in range(m):
        if Ac[fila, fila] == 0:
            '''
            Nos aseguramos de no estar en la última fila
            '''
            if fila + 1 < m: 
                intercambiarfilas(Ac, fila, fila + 1)
                P[fila] = fila + 1
                P[fila + 1] = fila 
            else:
                print("La matriz no tiene factorización LU.")
        '''Recorremos la matriz Ac. En cada paso, se calcula un factor 
        y se utiliza para restar las filas y obtener la eliminación gaussiana'''
        for i in range(fila + 1, m):
            factor = Ac[i, fila] / Ac[fila, fila]
            Ac[i, fila] = factor  
            Ac[i, fila + 1:] -= factor * Ac[fila, fila + 1:]
        '''Calculamos las matrices L y U que componen la factorización LU de la matriz original.
        L toma la parte triangular inferior estricta de la matriz Ac y le añadimos una matriz identidad''' 
        L = np.tril(Ac, -1) + np.eye(m) 
        U = np.triu(Ac) 
    return L, U, P

#inversaLU recibe la factorizacion LU y vector de permutacion y devuelve la inversa de A
def inversaLU (L, U,  P):
    filas, columnas = L.shape
    Inv = np.zeros((filas, columnas))  # Inicializa una matriz de ceros
    id = np.eye(filas)  # Crea una matriz identidad

    for i in range(columnas):
        y = sc.solve_triangular(L, id[:, i], lower=True)  # Resuelve L * y = e_i
        x = sc.solve_triangular(U, y)  # Resuelve U * x = y
        Inv[:, i] = x  # Almacena la columna en Inv
    Inv = Inv[:, P]
    return Inv
