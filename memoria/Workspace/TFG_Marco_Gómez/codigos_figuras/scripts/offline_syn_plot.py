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
parser.add_argument("plot_mode", type=int, choices=[-1, 0, 1, 2], help="Corriente a mostrar: -1 (ninguna), 0 (i_fast), 1 (i_slow), o 2 (i, i_fast e i_slow)")
parser.add_argument(
    "pre_is_live",
    type=int,
    choices=[0, 1],
    help="Flag de neurona presináptica: 0 (modelo) o 1 (viva)"
)
args = parser.parse_args()

data_path = args.input_csv
out_png = args.output_png
plot_mode = args.plot_mode
pre_is_live = args.pre_is_live
label_pre = 'Neur. presin. (viva)' if pre_is_live == 1 else 'Neur. presin. (modelo)'

if not os.path.exists(data_path):
    print("Archivo de datos no encontrado.")
    exit(1)

df = pd.read_csv(data_path)

n_plots = 1 if plot_mode == -1 else (2 if plot_mode in [0, 1] else 4)
heights = {1: 3.5, 2: 5.0, 4: 8.0}
fig, axs = plt.subplots(n_plots, 1, figsize=(
    6.5, heights[n_plots]), dpi=600, sharex=True)

# Asegurar que axs sea siempre un iterable
if n_plots == 1:
    axs = [axs]

title = 'Potenciales de membrana de sinapsis química unidireccional' if plot_mode == -1 else 'Potenciales de membrana y corriente sináptica de sinapsis química unidireccional'
fig.suptitle(title, fontsize=8)

lw = 0.6

# 1. v_pre y v_post
axs[0].plot(df['t'], df['v_post'], label='Neur. postsin. (modelo)',
            color='darkorange', linewidth=lw)
axs[0].plot(df['t'], df['v_pre'], label=label_pre,
            color='dodgerblue', linewidth=lw, alpha=0.6)
axs[0].set_ylabel('Voltaje (mV adim.)', fontsize=7)
axs[0].legend(fontsize=5)
axs[0].tick_params(axis='both', labelsize=7)

axs[0].autoscale(enable=True, axis='both', tight=True)

if plot_mode == 2:
    # 2. total i
    axs[1].plot(df['t'], df['i'], color='green', linewidth=lw)
    axs[1].set_ylabel('I (nA adim.)', fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    axs[1].autoscale(enable=True, axis='both', tight=True)

    # 3. i_fast
    axs[2].plot(df['t'], df['i_fast'], color='darkred', linewidth=lw)
    axs[2].set_ylabel('I_fast (nA adim.)', fontsize=7)
    axs[2].tick_params(axis='both', labelsize=7)
    axs[2].autoscale(enable=True, axis='both', tight=True)

    # 4. i_slow
    axs[3].plot(df['t'], df['i_slow'], color='purple', linewidth=lw)
    axs[3].set_ylabel('I_slow (nA adim.)', fontsize=7)
    axs[3].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[3].tick_params(axis='both', labelsize=7)
    axs[3].autoscale(enable=True, axis='both', tight=True)

elif plot_mode == 0:
    axs[1].plot(df['t'], df['i_fast'], color='darkred', linewidth=lw)
    axs[1].set_ylabel('I_fast (nA adim.)', fontsize=7)
    axs[1].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    axs[1].autoscale(enable=True, axis='both', tight=True)

elif plot_mode == 1:
    axs[1].plot(df['t'], df['i_slow'], color='purple', linewidth=lw)
    axs[1].set_ylabel('I_slow (nA adim.)', fontsize=7)
    axs[1].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    axs[1].autoscale(enable=True, axis='both', tight=True)

plt.tight_layout()
plt.savefig(out_png, dpi=600, bbox_inches='tight')
plt.close()
