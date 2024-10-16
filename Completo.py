import matplotlib.pyplot as plt
import scipy.linalg as sc
import numpy as np
import pandas as pd
import funciones as f

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


#Convierto las submatrices a arrays de numpy para trabajar mas comodo
A = A.to_numpy()
Pry_int = Pry_int.to_numpy()
Pry_ext = Pry_ext.to_numpy()
Nic_int = Nic_int.to_numpy()
Nic_ext = Nic_ext.to_numpy()
#Creo los vectores de produccion total para luego usar como P en la formula A = ZP^(-1)
Pry_out = Pry["Output"].to_numpy()
Pry_out_copy = Pry_out.copy()
Pry_out_copy[Pry_out_copy == 0] = 1 #remplazo 0 por 1
Nic_out = Nic["Output"].to_numpy()
Nic_out_copy = Nic_out.copy()
Nic_out_copy[Nic_out_copy == 0] = 1 #remplazo 0 por 1

#Coef intra-regionales PRY
z = Pry_int
p = Pry_out_copy
p = np.diag(p) #Diagonalizo el vector
l,u,per = f.calcularLU(p) #Lu de p
inv_p = f.inversaLU(l,u,per) #Inversa de p
res_intra_PRY = z@inv_p
#Coef intra-regionales NIC
z = Nic_int
p = Nic_out_copy
p = np.diag(p) #Diagonalizo el vector
l,u,per = f.calcularLU(p) #Lu de p
inv_p = f.inversaLU(l,u,per) #Inversa de p
res_intra_NIC = z@inv_p
#Coef inter-regionales PRY
z = Pry_int
p = Nic_out_copy
p = np.diag(p) #Diagonalizo el vector
l,u,per = f.calcularLU(p) #Lu de p
inv_p = f.inversaLU(l,u,per) #Inversa de p
res_inter_PRY = z@inv_p
#Coef inter-regionales NIC
z = Nic_int
p = Pry_out_copy
p = np.diag(p) #Diagonalizo el vector
l,u,per = f.calcularLU(p) #Lu de p
inv_p = f.inversaLU(l,u,per) #Inversa de p
res_inter_NIC = z@inv_p


#Shock
D_completa = ((np.identity(Pry_int.shape[0]) - res_intra_PRY) @ Pry_out_copy) - res_inter_PRY @Nic_out_copy #Calculo Demanda de Pry
D2_completa = D_completa.copy()
D2_completa[4] = D2_completa[4]*0.9
D2_completa[5] = D2_completa[5]*1.033
D2_completa[6] = D2_completa[6]*1.033
D2_completa[7] = D2_completa[7]*1.033
Delta_Demanda_Completa = D2_completa - D_completa # Diferencia en la demanda


#Calculo Delta_P con Delta_Demanda_Completa con la ecuacion de variacion de produccion considerando las relaciones inter-regionales
Id = np.identity(Pry_int.shape[0])
Id_p = Id - res_intra_PRY
Id_n = Id - res_intra_NIC
l,u,p = f.calcularLU(Id_n)
Id_n_inv = f.inversaLU(l,u,p)
multi = res_inter_NIC @ Id_n_inv @ res_inter_PRY
res = Id_p - multi
l,u,p = f.calcularLU(res)
res_inv = f.inversaLU(l,u,p)
Delta_Prod = res_inv @ Delta_Demanda_Completa

#Calculo Demanda simple
D_simple = ((np.identity(Pry_int.shape[0]) - res_intra_PRY) @ Pry_out_copy)
D2_simple = D_simple.copy()
D2_simple[4] = D2_simple[4]*0.9
D2_simple[5] = D2_simple[5]*1.033
D2_simple[6] = D2_simple[6]*1.033
D2_simple[7] = D2_simple[7]*1.033
Delta_Demanda_Simple = D2_simple - D_simple # Diferencia en la demanda

#Calculo Delta_P con la ecuacion del modelo simple
Id = np.identity(Pry_int.shape[0])
Id_p = Id - res_intra_PRY
l,u,p = f.calcularLU(Id_p)
Id_p_inv = f.inversaLU(l,u,p)
Delta_Prod_Simple = Id_p_inv @ Delta_Demanda_Simple



res_Delta = Delta_Prod - Delta_Prod_Simple

plt.figure(figsize=(20,5))

plt.bar(range(len(Delta_Prod)),Delta_Prod, 
color=np.where(Delta_Prod < 0, 'crimson', 'steelblue') #Color dependiendo de si es positivo o negativo
)

plt.xticks(range(len(Pry_col)),Pry_col,rotation=45)

plt.title('Variación de producción', fontsize=20, fontweight='bold')

# Agrego los valores encima de cada barra
for idx, value in enumerate(Delta_Prod):
    plt.text(idx, value + (0.01 if value >= 0 else -0.05), 
             f'{value:.2f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=9)

# Mejoro el estilo de los ejes
plt.xlabel('Sectores', fontsize=15)
plt.ylabel('Variación', fontsize=15)

plt.show()

plt.figure(figsize=(20,5))

plt.bar(range(len(res_Delta)),res_Delta, 
color=np.where(Delta_Prod < 0, 'crimson', 'steelblue') #Color dependiendo de si es positivo o negativo
)

plt.xticks(range(len(Pry_col)),Pry_col,rotation=45)

plt.title('Variación de producción', fontsize=20, fontweight='bold')

# Agrego los valores encima de cada barra
for idx, value in enumerate(res_Delta):
    plt.text(idx, value + (0.01 if value >= 0 else -0.05), 
             f'{value:.2f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=9)

# Mejoro el estilo de los ejes
plt.xlabel('Sectores', fontsize=15)
plt.ylabel('Variación', fontsize=15)

plt.show()