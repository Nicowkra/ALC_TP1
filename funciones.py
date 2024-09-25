import numpy as np
import pandas as pd
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
    filas, columnas = L.shape
    Inv = np.zeros((filas, columnas))  # Inicializa una matriz de ceros
    id = np.eye(filas)  # Crea una matriz identidad

    for i in range(columnas):
        y = sc.solve_triangular(L, id[:, i], lower=True)  # Resuelve L * y = e_i
        x = sc.solve_triangular(U, y)  # Resuelve U * x = y
        Inv[:, i] = x  # Almacena la columna en Inv

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
   
def crearMatrizA(matriz):
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
    A = crearMatrizA(matriz)
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
    Delta_Prod.plot(kind="bar",rot = 90,title ='Variación de producción',
                           color=np.where(Delta_Prod<0,'crimson','steelblue'),figsize=(20, 5))


    
def Leont2Reg(A,P): #Funcion de Leontief para 2 regiones, usando la formula (I-A)P = D
    m=A.shape[0] #filas A
    Id = np.identity(m)
    res = Id - A
    return res @ P
