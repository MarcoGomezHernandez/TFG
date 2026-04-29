import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import re

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

# 1. Parser de argumentos
parser = argparse.ArgumentParser(
    description="Graficar interacciones sinapticas bidireccionales.")
parser.add_argument("input_h5", help="Ruta al archivo HDF5 de entrada")
parser.add_argument("output_png", help="Ruta al archivo PNG de salida")
parser.add_argument("plot_mode_12", type=int, choices=[-1, 0, 1, 2],
                    help="Corriente a mostrar en la dirección 1->2: -1 (ninguna), 0 (i_fast), 1 (i_slow), o 2 (i, i_fast e i_slow)")
parser.add_argument("plot_mode_21", type=int, choices=[-1, 0, 1, 2],
                    help="Corriente a mostrar en la dirección 2->1: -1 (ninguna), 0 (i_fast), 1 (i_slow), o 2 (i, i_fast e i_slow)")
parser.add_argument("n1_is_live", type=int, choices=[0, 1],
                    help="Flag de la neurona 1: 0 (modelo) o 1 (viva)")
parser.add_argument("n2_is_live", type=int, choices=[0, 1],
                    help="Flag de la neurona 2: 0 (modelo) o 1 (viva)")
parser.add_argument("single_axis", type=int, choices=[0, 1],
                    help="1 para un único eje Y, 0 para ejes separados")
args = parser.parse_args()

data_path = args.input_h5
out_png = args.output_png
single_axis = args.single_axis == 1

# Etiquetas y sufijos de unidad según si la neurona es viva o modelo
label_n1 = "Neur. viva" if args.n1_is_live == 1 else "Neur. modelo"
label_n2 = "Neur. viva" if args.n2_is_live == 1 else "Neur. modelo"
unit_suffix_n1 = "" if args.n1_is_live == 1 else " adim."
unit_suffix_n2 = "" if args.n2_is_live == 1 else " adim."

# 2. Validar que el archivo existe
if not os.path.exists(data_path):
    print("Archivo de datos no encontrado.")
    exit(1)

# 3. Leer datos del archivo HDF5
cols = {}
with h5py.File(data_path, 'r') as f:
    trial_key = list(f.keys())[0] if 'Trial1' not in f else 'Trial1'
    sync = f[f'{trial_key}/Synchronous Data']
    period_ns = f[f'{trial_key}/Period (ns)'][()
                                              ] if 'Period (ns)' in f[trial_key] else 1000000
    dt_ms = period_ns * 1e-6
    data = sync['Channel Data'][:]

    for k in sync.keys():
        if k == 'Channel Data':
            continue
        m_idx = re.match(r'^(\d+)\s+', k)
        if not m_idx:
            continue
        idx = int(m_idx.group(1)) - 1

        m_type = re.search(r'(Vpre|Vpost|Current|Ifast|Islow)\s+(12|21)', k)
        if m_type:
            cols[f"{m_type.group(1).lower()}_{m_type.group(2)}"] = idx

t = np.arange(data.shape[0]) * dt_ms


def get_d(key):
    return data[:, cols[key]] if key in cols else np.zeros_like(t)


# 4. Elegir qué gráficas incluir según el modo de corriente
plots = ['v']
modes = [args.plot_mode_12, args.plot_mode_21]
if 2 in modes:
    plots.append('i')
if 0 in modes or 2 in modes:
    plots.append('ifast')
if 1 in modes or 2 in modes:
    plots.append('islow')

# 5. Configurar figura y ejes dinámicamente
n_plots = len(plots)
heights = {1: 3.5, 2: 5.0, 3: 6.5, 4: 8.0}
fig, axs = plt.subplots(n_plots, 1, figsize=(
    6.5, heights[n_plots]), dpi=600, sharex=True)

# Asegurar que axs sea siempre un iterable
if n_plots == 1:
    axs = [axs]

title = 'Potenciales de membrana y corriente sináptica de sinapsis química biidireccional en RTXI' if any(
    mode != -1 for mode in modes) else 'Potenciales de membrana de sinapsis química bidireccional en RTXI'
fig.suptitle(title, fontsize=8)
lw = 0.6
ax_idx = 0

# 6. Graficar voltajes
vpre12 = get_d('vpre_12')
vpost12 = get_d('vpost_12')

axs[ax_idx].plot(
    t, vpost12, label=f"{label_n2} (2)", color='darkorange', linewidth=lw)
axs[ax_idx].plot(t, vpre12, label=f"{label_n1} (1)",
                 color='dodgerblue', linewidth=lw, alpha=0.6)
axs[ax_idx].set_ylabel(f'Voltaje (mV{unit_suffix_n2})', fontsize=7)
axs[ax_idx].legend(fontsize=5)
axs[ax_idx].tick_params(axis='both', labelsize=7)

if not single_axis:
    vpre21 = get_d('vpre_21')
    vpost21 = get_d('vpost_21')

    ax_v_twin = axs[ax_idx].twinx()
    ax_v_twin.set_ylabel(f'Voltaje (mV{unit_suffix_n1})', fontsize=7)
    ax_v_twin.tick_params(axis='both', labelsize=7)

    axs[ax_idx].set_zorder(ax_v_twin.get_zorder() + 1)
    axs[ax_idx].patch.set_visible(False)

    ax_v_twin.autoscale(enable=True, axis='both', tight=True)

axs[ax_idx].autoscale(enable=True, axis='both', tight=True)

ax_idx += 1

# 7. Función auxiliar para graficar corrientes por dirección


def plot_current(ax, name, mode_target, c_dark, c_light, lbl):
    show_12 = args.plot_mode_12 in [mode_target, 2]
    show_21 = args.plot_mode_21 in [mode_target, 2]

    if not show_12 and not show_21:
        return

    lns = []

    # Condición de doble eje: solicitado por el usuario Y ambas señales presentes
    use_twin = (not single_axis) and show_12 and show_21
    use_legend = show_12 and show_21 and not use_twin

    # --- Eje Izquierdo ---
    if show_12:
        d12 = get_d(f'{name}_12')

        lbl_12 = "1->2" if use_legend else "_nolegend_"
        l1 = ax.plot(t, d12, label=lbl_12, color=c_dark,
                     linewidth=lw, alpha=0.6 if show_21 else 1.0)
        if use_legend:
            lns += l1

    if show_21 and not use_twin:
        d21 = get_d(f'{name}_21')

        lbl_21 = "2->1" if use_legend else "_nolegend_"
        l2 = ax.plot(t, d21, label=lbl_21, color=c_light, linewidth=lw)
        if use_legend:
            lns += l2

    # Ajuste de límites del eje izquierdo
    if show_12 and (show_21 and not use_twin):
        ax.set_ylabel(f"{lbl} (nA{unit_suffix_n2})", fontsize=7)
    elif show_12:
        ax.set_ylabel(f"{lbl} 1->2 (nA{unit_suffix_n2})", fontsize=7)

        # --- Eje Derecho (Twin) ---
        if show_21 and use_twin:
            d21 = get_d(f'{name}_21')

            ax_twin = ax.twinx()
            l2 = ax_twin.plot(
                t, d21, label="_nolegend_", color=c_light, linewidth=lw)

            ax_twin.set_ylabel(f"{lbl} 2->1 (nA{unit_suffix_n1})", fontsize=7)
            ax_twin.tick_params(axis='both', labelsize=7)

            ax.set_zorder(ax_twin.get_zorder() + 1)
            ax.patch.set_visible(False)

            ax_twin.autoscale(enable=True, axis='both', tight=True)
    elif show_21 and not use_twin:
        ax.set_ylabel(f"{lbl} 2->1 (nA{unit_suffix_n1})", fontsize=7)

    ax.autoscale(enable=True, axis='both', tight=True)
    ax.tick_params(axis='both', labelsize=7)

    if lns and use_legend:
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, fontsize=5)


if 'i' in plots:
    plot_current(axs[ax_idx], 'current', 2, 'darkgreen', 'limegreen', 'I')
    ax_idx += 1
if 'ifast' in plots:
    plot_current(axs[ax_idx], 'ifast', 0, 'darkred', 'red', 'I_fast')
    ax_idx += 1
if 'islow' in plots:
    plot_current(axs[ax_idx], 'islow', 1, 'indigo', 'mediumorchid', 'I_slow')
    ax_idx += 1

axs[-1].set_xlabel('Tiempo (ms)', fontsize=7)

plt.tight_layout()
plt.savefig(out_png, dpi=600, bbox_inches='tight')
plt.close()
