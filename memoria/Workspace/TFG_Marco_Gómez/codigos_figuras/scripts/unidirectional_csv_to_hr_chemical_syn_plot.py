import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

parser = argparse.ArgumentParser(
    description="Graficar interacciones sinapticas.")
parser.add_argument("input_csv", help="Ruta al archivo CSV de entrada")
parser.add_argument("output_png", help="Ruta al archivo PNG de salida")
parser.add_argument("plot_mode", choices=[
                    '0', '1', '2'], help="Corriente a mostrar: '0' (i_fast), '1' (i_slow), o '2' (i, i_fast e i_slow)")
parser.add_argument("label_pre", help="Etiqueta para la neurona presináptica")
args = parser.parse_args()

data_path = args.input_csv
out_png = args.output_png
plot_mode = args.plot_mode.lower()
label_pre = args.label_pre

if not os.path.exists(data_path):
    print("Data file not found!")
    exit(1)

df = pd.read_csv(data_path)

if plot_mode == '2':
    fig, axs = plt.subplots(4, 1, figsize=(6.5, 7.5), dpi=600, sharex=True)
elif plot_mode in ['0', '1']:
    fig, axs = plt.subplots(2, 1, figsize=(6.5, 4.5), dpi=600, sharex=True)
else:
    print(f"Unknown mode: {plot_mode}")
    exit(1)

fig.suptitle(
    'Potenciales de membrana y corriente sináptica de sinapsis química unidireccional', fontsize=8)

lw = 0.6

# 1. v_pre and v_post superimposed
axs[0].plot(df['t'], df['v_post'], label='Neurona postsináptica (modelo)',
            color='darkorange', linewidth=lw)
axs[0].plot(df['t'], df['v_pre'], label=label_pre,
            color='dodgerblue', linewidth=lw, alpha=0.6)
axs[0].set_ylabel('Voltaje (V adim.)', fontsize=7)
axs[0].legend(fontsize=7)
axs[0].tick_params(axis='both', labelsize=7)
axs[0].set_ylim(min(df['v_pre'].min(), df['v_post'].min()),
                max(df['v_pre'].max(), df['v_post'].max()))
plt.xlim(df['t'].min(), df['t'].max())

if plot_mode == '2':
    # 2. total i
    axs[1].plot(df['t'], df['i'], label='I', color='green', linewidth=lw)
    axs[1].set_ylabel('Corriente (uA adim.)', fontsize=7)
    axs[1].legend(fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    axs[1].set_ylim(df['i'].min(), df['i'].max())

    # 3. i_fast
    axs[2].plot(df['t'], df['i_fast'], label='I_fast',
                color='darkred', linewidth=lw)
    axs[2].set_ylabel('Corriente (uA adim.)', fontsize=7)
    axs[2].legend(fontsize=7)
    axs[2].tick_params(axis='both', labelsize=7)
    axs[2].set_ylim(df['i_fast'].min(), df['i_fast'].max())

    # 4. i_slow
    axs[3].plot(df['t'], df['i_slow'], label='I_slow',
                color='purple', linewidth=lw)
    axs[3].set_ylabel('Corriente (uA adim.)', fontsize=7)
    axs[3].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[3].legend(fontsize=7)
    axs[3].tick_params(axis='both', labelsize=7)
    axs[3].set_ylim(df['i_slow'].min(), df['i_slow'].max())

elif plot_mode == '0':
    axs[1].plot(df['t'], df['i_fast'], label='I_fast',
                color='orange', linewidth=lw)
    axs[1].set_ylabel('Corriente (uA adim.)', fontsize=7)
    axs[1].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[1].legend(fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    axs[1].set_ylim(df['i_fast'].min(), df['i_fast'].max())

elif plot_mode == '1':
    axs[1].plot(df['t'], df['i_slow'], label='I_slow',
                color='purple', linewidth=lw)
    axs[1].set_ylabel('Corriente (uA adim.)', fontsize=7)
    axs[1].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[1].legend(fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    axs[1].set_ylim(df['i_slow'].min(), df['i_slow'].max())

plt.tight_layout()
plt.savefig(out_png, dpi=600, bbox_inches='tight')
plt.close()
