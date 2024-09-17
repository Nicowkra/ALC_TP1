import pandas as pd 
matriz = pd.read_excel("matriz.xlsx", sheet_name ="LAC_IOT_2011",)

Nic_col = []    
Pry_col = []
for i in range(1,41): #Crea la lista de columnas a filtrar
    Nic_col.append('NICs'+str(i))
    Pry_col.append('PRYs'+str(i))
    
Pry_temp = matriz[matriz["Country_iso3"] == "PRY"] # Crea la tabla con filas de PRY
Nic_temp = matriz[matriz["Country_iso3"] == "NIC"] # Crea la tabla con filas de NIC

Pry_int= Pry_temp.loc[:,Pry_col] # Crea la matriz de demanda interna
Pry_ext = Pry_temp.loc[:,Nic_col] # Crea la matriz de demanda externa

Nic_int = Nic_temp.loc[:,Nic_col]  # Crea la matriz de demanda interna
Nic_ext = Nic_temp.loc[:,Pry_col]  # Crea la matriz de demanda externa

count_row = Pry_int.shape[0]  # cant de filas
count_col = Pry_int.shape[1]  # cant de columnas

Pry = pd.DataFrame(index = range(1,41), columns = Pry_col)
#df_Nic = 
for i in range(0,count_row):
    for j in range(0,count_col):
        Pry.iat[i,j] = Pry_int.iat[i,j] + Nic_ext.iat[i,j]

print(Pry)