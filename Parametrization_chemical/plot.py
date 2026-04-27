import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd

plt.rcParams.update({'font.size': 11})

if len(sys.argv) > 1:
    path = sys.argv[1]
else:
    # Fallback para pruebas si no se pasa argumento
    print("Error: No file specified. Usage: python plot.py <path_to_csv>")
    exit()

# Nuevo argumento: si es 1, combinar las dos primeras curvas en una gráfica
combine = int(sys.argv[2]) if len(sys.argv) > 2 else 0

file_name = os.path.basename(path)

print("Plotting file from", file_name)
data = pd.read_csv(path, delimiter=" ", low_memory=False)

headers = data.keys()
# rows es el total de columnas, pero el numero de plots es rows - 1 (quitando Time)
rows = data.shape[1]
num_plots = rows - 1 if combine == 0 else rows - 2  # Reducir si combinamos

colors = ['teal', 'brown', 'blue', 'green', 'maroon', 'teal', 'brown', 'blue', 'green', 'maroon']

# Cambio 1: constrained_layout maneja mejor los espacios y figsize más razonable
plt.figure(figsize=(10, 15), constrained_layout=True)

for i in range(1, rows):
    if combine == 1 and i == 2:
        continue  # Saltar la segunda curva, ya se grafica con la primera
    
    # Calcular índice de subplot
    subplot_idx = i if combine == 0 else (i - 1 if i > 2 else 1)
    
    if subplot_idx == 1:
        ax1 = plt.subplot(num_plots, 1, subplot_idx)
    else:
        plt.subplot(num_plots, 1, subplot_idx, sharex=ax1)

    col_name = headers[i]
    
    # Cambio 3: Detección de unidades (Voltaje vs Corriente)
    if 'V' in col_name:
        y_label = "Adim voltage"
    elif 'i' in col_name:
        y_label = "Adim current"
    else:
        y_label = "Value"

    if subplot_idx == num_plots:  # Último plot
        plt.xlabel("Time (ms)")

    plt.ylabel(y_label, multialignment='center')
    
    # Cambio 4: Plotear Tiempo vs Dato (no Indice vs Dato)
    plt.plot(data[headers[0]], data[col_name], color=colors[i-1])
    
    if combine == 1 and i == 1:
        # Graficar también la segunda curva en el mismo subplot
        col_name2 = headers[2]
        plt.plot(data[headers[0]], data[col_name2], color=colors[1])
        plt.legend([col_name, col_name2])
        plt.title(f"{col_name} and {col_name2}")
    else:
        plt.title(col_name)
    
    plt.grid(True, alpha=0.3)  # Opcional: Ayuda a cuadrar visualmente

# No es necesario tight_layout si usamos constrained_layout=True arriba
output_file = "./images/" + os.path.splitext(file_name)[0] + ".eps"
# Crear directorio si no existe para evitar error
os.makedirs("./images/", exist_ok=True) 

plt.savefig(output_file, format='eps')
plt.show()