import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

# Carga de datos generados por C++
df = pd.read_csv('data/hindmarsh-rose_modified.csv')

# Configuración y renderizado de la gráfica
plt.figure(figsize=(5.88, 2.94), dpi=300)
plt.plot(df['t'], df['x'], color='dodgerblue', linewidth=0.6)

plt.title('Hindmarsh-Rose Modificada', fontsize=8)
plt.xlabel('Tiempo/t (u. a.)', fontsize=7)
plt.ylabel('Voltaje/x(t) (u. a.)', fontsize=7)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

ax = plt.gca()
ax.margins(x=0.0, y=0.02)
plt.tight_layout()
plt.savefig('figures/hindmarsh-rose_modified.png',
            dpi=300, bbox_inches='tight')
