import pandas as pd 
import numpy as np
import scipy.linalg as sc

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

#Crea vectores de produccion total
Pry_out = Pry["Output"]
Pry_out = Pry_out.replace(0,1) #remplazo 0 por 1
Nic_out = Nic["Output"]
Nic_out = Nic_out.replace(0,1) #remplazo 0 por 1

#----Coeficientes Tecnicos----#
def coefTec(z,p):
    p = np.diag(p.values) #Diagonalizo el vector
    per,l,u = sc.lu(p) #Lu de p
    inv_p = sc.inv(per@l@u) #Inversa de p
    return z@inv_p

#Coef intra-regionales
cT_NxN = coefTec (Nic_int,Nic_out)
cT_PxP = coefTec (Pry_int,Pry_out)
#Coef intre-regionales
cT_NxP = coefTec (Nic_int,Pry_out)
cT_PxN = coefTec (Pry_int,Nic_out)



count_row = Pry_int.shape[0]  # cant de filas
count_col = Pry_int.shape[1]  # cant de columnas

