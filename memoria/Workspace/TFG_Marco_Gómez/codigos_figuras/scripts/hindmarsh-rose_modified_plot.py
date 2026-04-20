import pandas as pd
import matplotlib.pyplot as plt

# Carga de datos generados por C++
df = pd.read_csv('data/hindmarsh-rose_modified.csv')

# Configuración y renderizado de la gráfica
plt.figure(figsize=(6.5, 3.5), dpi=600)
plt.plot(df['t'], df['x'], color='dodgerblue', linewidth=0.8)

plt.title('Hindmarsh-Rose Modificada', fontsize=8)
plt.xlabel('Tiempo/t (unidades adimensionales)', fontsize=7)
plt.ylabel('Voltaje/x(t) (unidades adimensionales)', fontsize=7)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.xlim(df['t'].min(), df['t'].max())
plt.ylim(df['x'].min(), df['x'].max())
plt.tight_layout()
plt.savefig('figures/hindmarsh-rose_modified.png', dpi=600, bbox_inches='tight')
