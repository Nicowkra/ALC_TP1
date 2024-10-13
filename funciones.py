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
def intercambiarfilas(A, fila1, fila2):
    A[[fila1, fila2]] = A[[fila2, fila1]]
    return A
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

# =============================================================================
# --
# =============================================================================
                
def crearMatrizA():
    matriz = pd.read_excel("matriz.xlsx", sheet_name ="LAC_IOT_2011",)
    Nic_col = []    
    Pry_col = []
    for i in range(1,41): #Crea la lista de columnas a filtrar
        Nic_col.append('NICs'+str(i))
        Pry_col.append('PRYs'+str(i))
        
    Pry = matriz[matriz["Country_iso3"] == "PRY"] # Crea la tabla con filas de PRY
    Nic = matriz[matriz["Country_iso3"] == "NIC"] # Crea la tabla con filas de NIC
    
    
    # Crea matrices intra-regionales
    Pry_int= Pry.loc[:,Pry_col] 
    Nic_int = Nic.loc[:,Nic_col] 
    
    #Crea matrices intre-regionales
    Nic_ext = Nic.loc[:,Pry_col] 
    Pry_ext = Pry.loc[:,Nic_col]
    
    #Columna con los nombres de los sectores para despues mantener los indices
    colnames = pd.DataFrame({'Sectores':Pry_col + Nic_col})
    colnamesPry = pd.DataFrame({'Sectores':Pry_col})
    colnamesNic = pd.DataFrame({'Sectores':Nic_col})
    # Se cambian los indices a los nombres del sector
    Pry_int.index = colnamesPry['Sectores']
    Nic_int.index = colnamesNic['Sectores']
    Nic_ext.index = colnamesNic['Sectores']
    Pry_ext.index = colnamesPry['Sectores']
    
    #Concateno las submatrices para crear mi A
    A_Pry = pd.concat([Pry_int,Nic_ext])
    A_Nic = pd.concat([Pry_ext,Nic_int])
    A = pd.concat([A_Pry,A_Nic], axis=1)
    return A
def coefTec(inter,out):
    matriz = pd.read_excel("matriz.xlsx", sheet_name ="LAC_IOT_2011",)
    Nic_col = []    
    Pry_col = []
    for i in range(1,41): #Crea la lista de columnas a filtrar
        Nic_col.append('NICs'+str(i))
        Pry_col.append('PRYs'+str(i))
    Pry = matriz[matriz["Country_iso3"] == "PRY"] # Crea la tabla con filas de PRY
    Nic = matriz[matriz["Country_iso3"] == "NIC"] # Crea la tabla con filas de NIC

    Pry_int= Pry.loc[:,Pry_col] 
    Nic_int = Nic.loc[:,Nic_col]

    # Se cambian los indices a los nombres del sector
    colnamesPry = pd.DataFrame({'Sectores':Pry_col})
    colnamesNic = pd.DataFrame({'Sectores':Nic_col})
    Pry_int.index = colnamesPry['Sectores']
    Nic_int.index = colnamesNic['Sectores']
    
    #Creo los vectores de produccion total para luego usar como P en la formula A = ZP^(-1)
    Pry_out = Pry["Output"]
    Pry_out = Pry_out.replace(0,1) #remplazo 0 por 1
    Nic_out = Nic["Output"]
    Nic_out = Nic_out.replace(0,1) #remplazo 0 por 1
    
    if inter == "Nic": 
        z = Nic_int
    else:
        z = Nic_int
        
    if out == "Nic": 
        p = Nic_out
    else:
        p = Pry_out
    A = calcCoefTec(z,p)


def calcCoefTec(z,p):
    p = np.diag(p.values) #Diagonalizo el vector
    per,l,u = sc.lu(p) #Lu de p
    inv_p = sc.inv(per@l@u) #Inversa de p
    return z@inv_p
def shock():
    matriz = pd.read_excel("matriz.xlsx", sheet_name ="LAC_IOT_2011",)
    A = crearMatrizA()
    Pry = matriz[matriz["Country_iso3"] == "PRY"] # Crea la tabla con filas de PRY
    Nic = matriz[matriz["Country_iso3"] == "NIC"] # Crea la tabla con filas de NIC
    Nic_col = []    
    Pry_col = []
    for i in range(1,41): #Crea la lista de columnas a filtrar
        Nic_col.append('NICs'+str(i))
        Pry_col.append('PRYs'+str(i))
    #Creo los vectores de produccion total para luego usar como P en la formula A = ZP^(-1)
    # Crea matrices intra-regionales
    Pry_int= Pry.loc[:,Pry_col] 
    Nic_int = Nic.loc[:,Nic_col] 
        
    #Crea matrices intre-regionales
    Nic_ext = Nic.loc[:,Pry_col] 
    Pry_ext = Pry.loc[:,Nic_col]

    Pry_out = Pry["Output"]
    Pry_out = Pry_out.replace(0,1) #remplazo 0 por 1
    Nic_out = Nic["Output"]
    Nic_out = Nic_out.replace(0,1) #remplazo 0 por 1

    colnames = pd.DataFrame({'Sectores':Pry_col + Nic_col})
    colnamesPry = pd.DataFrame({'Sectores':Pry_col})
    colnamesNic = pd.DataFrame({'Sectores':Nic_col})
    # Se cambian los indices a los nombres del sector
    Pry_int.index = colnamesPry['Sectores']
    Nic_int.index = colnamesNic['Sectores']
    Nic_ext.index = colnamesNic['Sectores']
    Pry_ext.index = colnamesPry['Sectores']

        
    P1 = pd.concat([Pry_out,Nic_out]) #Vector P
    P1.index = colnames['Sectores']
           
    D1 = Leont2Reg(A,P1) # Demanda para las dos regiones originales
    D2 = D1.copy()
    D2["PRYs5"] = D2["PRYs5"]*0.9
    D2["PRYs6"] = D2["PRYs6"]*1.033
    D2["PRYs7"] = D2["PRYs7"]*1.033
    D2["PRYs8"] = D2["PRYs8"]*1.033
    Delta_Demanda = D2 - D1 # Diferencia en la demanda
    Delta_Demanda = Delta_Demanda.loc[Pry_col] #Para usar Delta_Demanda en el calculo de Delta_P tiene que ser solo de PRY

    #Calculo Delta_P con Delta_Demanda con la ecuacion de variacion de produccion considerando las relaciones inter-regionales
    Id = np.identity(Pry_int.shape[0])
    Id_p = Id - Pry_int
    Id_n = Id - Nic_int
    _Id_n = pd.DataFrame(sc.inv(Id_n),columns=Nic_col,index=Nic_col) #Lo invierto y convierto en Dataframe para tener los mismos objetos
    res = Id_p - (Pry_ext @ _Id_n @ Nic_ext)
    _res = pd.DataFrame(sc.inv(res),columns=Pry_col,index=Pry_col) #Lo invierto y convierto en Dataframe para tener los mismos objetos
    Delta_Prod = _res @ Delta_Demanda

    #Delta_Prod = P2 - P1 # Diferencia en la producción
    
    Delta_Prod.plot(
        kind="bar", 
        rot=45, 
        title='Variación de producción', 
        color=np.where(Delta_Prod < 0, 'crimson', 'steelblue'), 
        figsize=(20, 5)
    )
    
    # Mejoro la letra del título
    plt.title('Variación de producción', fontsize=20, fontweight='bold')
    
    # Agrego los valores encima de cada barra
    for idx, value in enumerate(Delta_Prod):
        plt.text(idx, value + (0.01 if value >= 0 else -0.05), 
                 f'{value:.2f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=9)
    
    # Mejoro el estilo de los ejes
    plt.xlabel('Sectores', fontsize=15)
    plt.ylabel('Variación', fontsize=15)
    
    # Ajusto los márgenes para que no se corte el gráfico
    plt.tight_layout()

    # Muestro el gráfico
    plt.show()

def regionSimple():
    matriz = pd.read_excel("matriz.xlsx", sheet_name ="LAC_IOT_2011",)
    A = crearMatrizA()
    Pry = matriz[matriz["Country_iso3"] == "PRY"] # Crea la tabla con filas de PRY
    Nic = matriz[matriz["Country_iso3"] == "NIC"] # Crea la tabla con filas de NIC
    Nic_col = []    
    Pry_col = []
    for i in range(1,41): #Crea la lista de columnas a filtrar
        Nic_col.append('NICs'+str(i))
        Pry_col.append('PRYs'+str(i))
    #Creo los vectores de produccion total para luego usar como P en la formula A = ZP^(-1)
    # Crea matrices intra-regionales
    Pry_int= Pry.loc[:,Pry_col] 
    Nic_int = Nic.loc[:,Nic_col] 
        
    #Crea matrices intre-regionales
    Nic_ext = Nic.loc[:,Pry_col] 
    Pry_ext = Pry.loc[:,Nic_col]

    Pry_out = Pry["Output"]
    Pry_out = Pry_out.replace(0,1) #remplazo 0 por 1
    Nic_out = Nic["Output"]
    Nic_out = Nic_out.replace(0,1) #remplazo 0 por 1

    colnames = pd.DataFrame({'Sectores':Pry_col + Nic_col})
    colnamesPry = pd.DataFrame({'Sectores':Pry_col})
    colnamesNic = pd.DataFrame({'Sectores':Nic_col})
    # Se cambian los indices a los nombres del sector
    Pry_int.index = colnamesPry['Sectores']
    Nic_int.index = colnamesNic['Sectores']
    Nic_ext.index = colnamesNic['Sectores']
    Pry_ext.index = colnamesPry['Sectores']

        
    P1 = pd.concat([Pry_out,Nic_out]) #Vector P
    P1.index = colnames['Sectores']
           
    D1 = Leont2Reg(A,P1) # Demanda para las dos regiones originales
    D2 = D1.copy()
    D2["PRYs5"] = D2["PRYs5"]*0.9
    D2["PRYs6"] = D2["PRYs6"]*1.033
    D2["PRYs7"] = D2["PRYs7"]*1.033
    D2["PRYs8"] = D2["PRYs8"]*1.033
    Delta_Demanda = D2 - D1 # Diferencia en la demanda
    Delta_Demanda = Delta_Demanda.loc[Pry_col]
    
    #Calculo Delta_P_Simple con la ecuacion del modelo de región simple
    Id = np.identity(Pry_int.shape[0])
    A = Id - Pry_int
    _A = pd.DataFrame(sc.inv(A),columns=Pry_col,index=Pry_col)
    Delta_P_Simple = _A@ Delta_Demanda
    Delta_P_Simple.plot(
        kind="bar", 
        rot=45, 
        title='Variación de producción', 
        color=np.where(Delta_P_Simple < 0, 'crimson', 'steelblue'), 
        figsize=(20, 5)
    )
    
    # Mejoro la letra del título
    plt.title('Variación de producción', fontsize=20, fontweight='bold')
    
    # Agrego los valores encima de cada barra
    for idx, value in enumerate(Delta_P_Simple):
        plt.text(idx, value + (0.01 if value >= 0 else -0.05), 
                 f'{value:.2f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=9)
    
    # Mejoro el estilo de los ejes
    plt.xlabel('Sectores', fontsize=15)
    plt.ylabel('Variación', fontsize=15)
    
    # Ajusto los márgenes para que no se corte el gráfico
    plt.tight_layout()

    # Muestro el gráfico
    plt.show()
    
        
def Leont2Reg(A,P): #Funcion de Leontief para 2 regiones, usando la formula (I-A)P = D
    m=A.shape[0] #filas A
    Id = np.identity(m)
    res = Id - A
    return res @ P


# =============================================================================
# FUNCION PARA CONSIGNA 5
# =============================================================================

def construirZyP():
    # Construyo Z y P
    Z = pd.DataFrame({'S1':[350,50,200],'S2':[0,250,150],'S3':[0,150,550]})
    P = pd.Series({'S1':1000,'S2':500,'S3':1000}) 

    A = calcCoefTec(Z, P) # Funcion que calcula los coeficientes tecnicos 

    L = np.identity(A.shape[0]) - A # Busco la matriz de Leontief
    L = sc.inv(L)

    print("Imprimimos la matriz A:\n", A)
    print("Imprimimos la matriz L de Leontief:\n", L)
