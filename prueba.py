import pandas as pd 
import numpy as np
import funciones as f

matriz = pd.read_excel("matriz.xlsx", sheet_name ="LAC_IOT_2011",)

Nic_col = []    
Pry_col = []
for i in range(1,41): #Crea la lista de columnas a filtrar
    Nic_col.append('NICs'+str(i))
    Pry_col.append('PRYs'+str(i))
    
Pry = matriz[matriz["Country_iso3"] == "PRY"] # Crea la tabla con filas de PRY
Nic = matriz[matriz["Country_iso3"] == "NIC"] # Crea la tabla con filas de NIC

#Columna con los nombres de los sectores para despues mantener los indices
colnames = pd.DataFrame({'Sectores':Pry_col + Nic_col})
colnamesPry = pd.DataFrame({'Sectores':Pry_col})
colnamesNic = pd.DataFrame({'Sectores':Nic_col})

# Crea matrices intra-regionales
Pry_int= Pry.loc[:,Pry_col] 
Nic_int = Nic.loc[:,Nic_col] 

#Crea matrices intre-regionales
Nic_ext = Nic.loc[:,Pry_col] 
Pry_ext = Pry.loc[:,Nic_col]

# Se cambian los indices a los nombres del sector
Pry_int.index = colnamesPry['Sectores']
Nic_int.index = colnamesNic['Sectores']
Nic_ext.index = colnamesNic['Sectores']
Pry_ext.index = colnamesPry['Sectores']

A1 = pd.concat([Pry_int,Nic_ext])
A2 = pd.concat([Pry_ext,Nic_int])
A = pd.concat([A1,A2], axis=1)

#Crea vectores de produccion total
Pry_out = Pry["Output"]
Pry_out = Pry_out.replace(0,1) #remplazo 0 por 1
Nic_out = Nic["Output"]
Nic_out = Nic_out.replace(0,1) #remplazo 0 por 1



#----Coeficientes Tecnicos----#

#Coef intra-regionales
cT_NxN = f.coefTec (Nic_int,Nic_out)
#cT_PxP = f.coefTec (Pry_int,Pry_out)
#Coef intre-regionales
#cT_NxP = f.coefTec (Nic_int,Pry_out)
#cT_PxN = f.coefTec (Pry_int,Nic_out)

P1 = pd.concat([Pry_out,Nic_out]) #Vector P
P1.index = colnames['Sectores']

P2 = P1.copy() #Vector P2 con las variaciones del shock aplicadas
P2['PRYs5'] = P2['PRYs5']*0.9
P2['PRYs6'] = P2['PRYs6']*1.033
P2['PRYs7'] = P2['PRYs7']*1.033
P2['PRYs8'] = P2['PRYs8']*1.033

D1 = f.Leont2Reg(A,P1) # Demanda para las dos regiones originales
D2 = f.Leont2Reg(A,P2) # Demanda para las dos regiones con los shocks aplicados

Delta_Demanda = D2 - D1

df = pd.DataFrame({'Original':D1,'Variación':D2},index = colnames['Sectores'])
df_Delta = pd.Series(Delta_Demanda,index = colnames['Sectores'])
df_Delta.plot.bar(rot = 90,title ='Variación de demanda en unidades',
                       color=np.where(df_Delta<0,'crimson','steelblue'),figsize=(20, 5))
df.plot.bar(rot = 90,title ='Comparación de demanda',figsize=(20, 5))

