import pandas as pd
import matplotlib.pyplot as plt

# Carga de datos generados por C++
df = pd.read_csv('data/hindmarsh_rose_improved.csv')

# Configuración y renderizado de la gráfica
plt.figure(figsize=(12, 5), dpi=300)
plt.plot(df['t'], df['x'], color='dodgerblue', linewidth=0.8)

plt.title('Hindmarsh Rose Mejorada', fontsize=10)
plt.xlabel('Tiempo (unidades adimensionales)', fontsize=8)
plt.ylabel('Voltaje/x (unidades adimensionales)', fontsize=8)

plt.xlim(df['t'].min(), df['t'].max())
plt.ylim(df['x'].min(), df['x'].max())
plt.tight_layout()
plt.savefig('figures/hindmarsh_rose_improved.png', dpi=300, bbox_inches='tight')
