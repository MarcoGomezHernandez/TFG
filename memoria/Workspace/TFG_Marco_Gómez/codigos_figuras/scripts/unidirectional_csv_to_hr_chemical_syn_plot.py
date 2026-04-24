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
parser.add_argument(
    "pre_is_live",
    choices=['0', '1'],
    help="Flag de neurona presináptica: '0' (modelo) o '1' (viva)"
)
args = parser.parse_args()

data_path = args.input_csv
out_png = args.output_png
plot_mode = args.plot_mode.lower()
pre_is_live = args.pre_is_live
label_pre = 'Voltaje neurona presináptica (viva)' if pre_is_live == '1' else 'Voltaje neurona presináptica (modelo)'

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
v_min = min(df['v_pre'].min(), df['v_post'].min())
v_max = max(df['v_pre'].max(), df['v_post'].max())

axs[0].plot(df['t'], df['v_post'], label='Voltaje neurona postsináptica (modelo)',
            color='darkorange', linewidth=lw)
axs[0].plot(df['t'], df['v_pre'], label=label_pre,
            color='dodgerblue', linewidth=lw, alpha=0.6)
axs[0].set_ylabel('Voltaje (mV adim.)', fontsize=7)
axs[0].legend(fontsize=7)
axs[0].tick_params(axis='both', labelsize=7)

if v_min != v_max:
    axs[0].set_ylim(v_min, v_max)

plt.xlim(df['t'].min(), df['t'].max())

if plot_mode == '2':
    # 2. total i
    i_min = df['i'].min()
    i_max = df['i'].max()
    axs[1].plot(df['t'], df['i'], label='I', color='green', linewidth=lw)
    axs[1].set_ylabel('Corriente (nA adim.)', fontsize=7)
    axs[1].legend(fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    if i_min != i_max:
        axs[1].set_ylim(i_min, i_max)

    # 3. i_fast
    i_fast_min = df['i_fast'].min()
    i_fast_max = df['i_fast'].max()
    axs[2].plot(df['t'], df['i_fast'], label='I_fast',
                color='darkred', linewidth=lw)
    axs[2].set_ylabel('Corriente (nA adim.)', fontsize=7)
    axs[2].legend(fontsize=7)
    axs[2].tick_params(axis='both', labelsize=7)
    if i_fast_min != i_fast_max:
        axs[2].set_ylim(i_fast_min, i_fast_max)

    # 4. i_slow
    i_slow_min = df['i_slow'].min()
    i_slow_max = df['i_slow'].max()
    axs[3].plot(df['t'], df['i_slow'], label='I_slow',
                color='purple', linewidth=lw)
    axs[3].set_ylabel('Corriente (nA adim.)', fontsize=7)
    axs[3].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[3].legend(fontsize=7)
    axs[3].tick_params(axis='both', labelsize=7)
    if i_slow_min != i_slow_max:
        axs[3].set_ylim(i_slow_min, i_slow_max)

elif plot_mode == '0':
    i_fast_min = df['i_fast'].min()
    i_fast_max = df['i_fast'].max()
    axs[1].plot(df['t'], df['i_fast'], label='I_fast',
                color='darkred', linewidth=lw)
    axs[1].set_ylabel('Corriente (nA adim.)', fontsize=7)
    axs[1].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[1].legend(fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    if i_fast_min != i_fast_max:
        axs[1].set_ylim(i_fast_min, i_fast_max)

elif plot_mode == '1':
    i_slow_min = df['i_slow'].min()
    i_slow_max = df['i_slow'].max()
    axs[1].plot(df['t'], df['i_slow'], label='I_slow',
                color='purple', linewidth=lw)
    axs[1].set_ylabel('Corriente (nA adim.)', fontsize=7)
    axs[1].set_xlabel('Tiempo (ms)', fontsize=7)
    axs[1].legend(fontsize=7)
    axs[1].tick_params(axis='both', labelsize=7)
    if i_slow_min != i_slow_max:
        axs[1].set_ylim(i_slow_min, i_slow_max)

plt.tight_layout()
plt.savefig(out_png, dpi=600, bbox_inches='tight')
plt.close()
