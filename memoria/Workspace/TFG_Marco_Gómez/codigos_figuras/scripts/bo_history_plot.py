import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

# Analiza un historial de optimización bayesiana (BO) en JSONL y genera gráficos PNG.
parser = argparse.ArgumentParser(
    description="Genera 4 gráficas a partir del historial BO en formato JSONL.")
parser.add_argument(
    "input_jsonl", help="Ruta al archivo JSONL de enrada")
parser.add_argument(
    "--output_max_score", help="Ruta al archio PNG de salida para la gráfica de máximos encontrados")
parser.add_argument(
    "--output_avg_scores", help="Ruta al archivo PNG de salida para la gráfica de promedios de puntuaciones")
parser.add_argument(
    "--output_scores", help="Ruta al archivo PNG de salida para la gráfica de puntuaciones de una única ejecución")
parser.add_argument(
    "--output_ard_lss", help="Ruta al archivo PNG de salida para la gráfica del inverso de lengthscales")
parser.add_argument(
    "--output_avg_execution_time", help="Ruta al archivo TXT de salida para guardar el tiempo de optimización medio y su desviación")
parser.add_argument(
    "--output_avg_best_params", help="Ruta al archivo CSV de salida para guardar la media y desviación de los mejores parámetros")
parser.add_argument("--row_index", type=int, default=0,
                    help="Índice de la fila/ejecución para la gráfica de las puntuaciones (por defecto 0)")
parser.add_argument("--initial_samples", type=int,
                    help="Número de muestras iniciales para dibujar una línea divisoria")
parser.add_argument("--xtick_step", type=int,
                    help="Paso de los ticks mayores en el eje x para todas las gráficas de puntuaciones")
parser.add_argument("--xminor_step", type=int,
                    help="Paso de los ticks menores en el eje x para la gráfica de puntuaciones individuales")

args = parser.parse_args()

need_only_single_run = bool(
    args.output_scores and not args.output_max_score
    and not args.output_avg_scores and not args.output_ard_lss
    and not args.output_avg_execution_time and not args.output_avg_best_params)

# Determina qué datos deben cargarse según los gráficos que se pedirán.
need_scores = args.output_max_score or args.output_avg_scores or args.output_scores
need_score_components = args.output_avg_scores or args.output_scores
need_ard_lss = bool(args.output_ard_lss)
need_avg_execution_time = bool(args.output_avg_execution_time)
need_avg_best_params = bool(args.output_avg_best_params)

if need_scores:
    scores_list = []
if need_score_components:
    range_scores_list = []
    shape_scores_list = []
if need_ard_lss:
    ard_lss_list = []
if need_avg_execution_time:
    opt_times_list = []
if need_avg_best_params:
    best_params_list = []

with open(args.input_jsonl, 'r') as f:
    # Recorre cada línea JSONL y procesa los resultados de cada ejecución.
    for idx, line in enumerate(f):
        line = line.strip()
        if not line:
            continue

        data = json.loads(line)

        if need_only_single_run:
            if idx == args.row_index:
                score_history = data["score_history"]
                scores_list.append(score_history["scores"])
                range_scores_list.append(score_history["range_scores"])
                shape_scores_list.append(score_history["shape_scores"])
                break
        else:
            if need_scores:
                scores_list.append(data["score_history"]["scores"])
            if need_score_components:
                score_history = data["score_history"]
                range_scores_list.append(score_history["range_scores"])
                shape_scores_list.append(score_history["shape_scores"])
            if need_ard_lss:
                ard_lss_list.append(data["ARD_lss"])
            if need_avg_execution_time:
                opt_times_list.append(data["optimization_time"])
            if need_avg_best_params:
                best_params_list.append(data["best_params"])

if need_scores:
    # Convierte las listas recopiladas a arreglos NumPy para facilitar los cálculos.
    scores = np.array(scores_list)
    if need_score_components:
        range_scores = np.array(range_scores_list)
        shape_scores = np.array(shape_scores_list)
    n_execs, n_evals = scores.shape
    eval_indices = range(1, n_evals + 1)
else:
    if need_ard_lss:
        n_execs = len(ard_lss_list)
    elif need_avg_execution_time:
        n_execs = len(opt_times_list)
    elif need_avg_best_params:
        n_execs = len(best_params_list)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

lw = 0.6
markersize = 1.2
margin = 0.01


def setup_plot(title, xlabel, ylabel):
    # Configura el estilo básico del gráfico antes de trazar los datos.
    plt.figure(figsize=(6.5, 3.5), dpi=600)
    plt.title(title, fontsize=8)
    plt.xlabel(xlabel, fontsize=7)
    plt.ylabel(ylabel, fontsize=7)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)


def save_plot(path):
    # Ajusta el diseño y guarda la figura en archivo.
    plt.tight_layout()
    plt.savefig(path, dpi=600, bbox_inches='tight')
    plt.close()


def set_extras(xtick_step, legend_fontsize=7):
    if xtick_step:
        plt.xticks(range(1, n_evals + 2, xtick_step))
    else:
        plt.xticks(range(1, n_evals + 2, max(1, n_evals // 10)))
    plt.margins(margin)
    plt.grid(True, linewidth=0.2)
    plt.gca().set_axisbelow(True)
    if args.initial_samples:
        plt.axvline(x=args.initial_samples + 1, color='black', linestyle='--',
                    linewidth=lw, label='Fin de muestras iniciales', zorder=0.51)
    plt.legend(fontsize=legend_fontsize)


def plot_with_shade(mean, std, label, color, alpha):
    # Traza la media de una serie con su desviación estándar como región sombreada.
    plt.plot(eval_indices, mean, label=label,
             color=color, linewidth=lw, alpha=alpha)
    mean_minus_std = mean - std
    mean_plus_std = mean + std
    plt.fill_between(eval_indices, mean_minus_std, mean_plus_std, color=color,
                     alpha=alpha/6, edgecolor='none')
    plt.plot(eval_indices, mean_minus_std,
             color=color, alpha=alpha/3, linewidth=lw/2)
    plt.plot(eval_indices, mean_plus_std,
             color=color, alpha=alpha/3, linewidth=lw/2)


if n_execs > 0:
    # Solo se crean figuras cuando hay ejecuciones y evaluaciones disponibles.
    if n_evals > 0:
        if args.output_max_score:
            # Calcula la puntuación máxima acumulada por evaluación y su media.
            max_scores = np.maximum.accumulate(scores, axis=1)

            setup_plot('Puntuación máxima acumulada por evaluación (convergencia)' if n_execs == 1 else 'Promedio de la puntuación máxima acumulada por evaluación (convergencia)',
                       'Evaluación', 'Puntuación máxima acumulada')
            if n_execs > 1:
                mean_max_scores = np.mean(max_scores, axis=0)
                std_max_scores = np.std(max_scores, axis=0)
                plot_with_shade(mean_max_scores,
                                std_max_scores, '_nolegend_', 'dodgerblue', 0.6)
            else:
                plt.plot(eval_indices, max_scores[0], label='_nolegend_',
                         color='dodgerblue', linewidth=lw, alpha=0.6)
            set_extras(args.xtick_step)
            save_plot(args.output_max_score)

        if args.output_avg_scores:
            mean_scores = np.mean(scores, axis=0)
            std_scores = np.std(scores, axis=0)

            mean_range = np.mean(range_scores, axis=0)
            std_range = np.std(range_scores, axis=0)

            mean_shape = np.mean(shape_scores, axis=0)
            std_shape = np.std(shape_scores, axis=0)

            setup_plot('Promedio de la puntuación y sus componentes por evaluación',
                       'Evaluación', 'Puntuación')
            plot_with_shade(mean_scores,
                            std_scores, 'Puntuación total', 'C0', 0.6)
            plot_with_shade(mean_range,
                            std_range, 'Puntuación del rango de la corriente', 'C1', 0.6)
            plot_with_shade(mean_shape,
                            std_shape, 'Puntuación de la forma', 'C2', 0.6)
            set_extras(args.xtick_step, legend_fontsize=5)
            save_plot(args.output_avg_scores)

        if args.output_scores:
            row_index = 0 if need_only_single_run else args.row_index
            if row_index < n_execs:
                setup_plot('Puntuación y sus componentes por evaluación',
                           'Evaluación', 'Puntuación')
                plt.plot(eval_indices, scores[row_index], 'o', label='Puntuación total',
                         color='C0', markersize=markersize, linestyle='None', alpha=0.6)
                plt.plot(eval_indices, range_scores[row_index], 'o', label='Puntuación del rango',
                         color='C1', markersize=markersize, linestyle='None', alpha=0.6)
                plt.plot(eval_indices, shape_scores[row_index], 'o', label='Puntuación de la forma',
                         color='C2', markersize=markersize, linestyle='None', alpha=0.6)
                set_extras(args.xtick_step, legend_fontsize=5)
                if args.xtick_step:
                    plt.gca().set_xticks(range(1, n_evals + 2, args.xminor_step), minor=True)
                    plt.grid(True, which='minor', linewidth=0.2)
                save_plot(args.output_scores)
            else:
                print(
                    f"Advertencia: el índice de fila {row_index} es mayor o igual que el número de ejecuciones disponibles ({n_execs}). No se generará la gráfica de puntuaciones individuales.")
    else:
        print("Advertencia: no se encontraron evaluaciones en los datos para generar las gráficas de puntuaciones.")

    if args.output_ard_lss:
        # Cada entrada ARD_lss contiene las lengthscales de los parámetros.
        if len(ard_lss_list) > 0:
            # Aplanar el diccionario para soportar tanto la versión offline (plano) como online (anidado por dirección)
            flat_ard_lss_list = []
            for lss in ard_lss_list:
                flat_ard_lss = {}
                for k, v in lss.items():
                    if isinstance(v, dict):
                        for sub_k, sub_v in v.items():
                            flat_ard_lss[f"{k} {sub_k}"] = sub_v
                    else:
                        flat_ard_lss[k] = v
                flat_ard_lss_list.append(flat_ard_lss)

            keys = sorted(list(flat_ard_lss_list[0].keys()))

            if len(keys) > 0:
                inv_ls_data = {k: [] for k in keys}
                for lss in flat_ard_lss_list:
                    for k in keys:
                        val = lss[k]
                        if val == 0:
                            inv_val = 1e10
                        else:
                            inv_val = 1.0 / val
                        inv_ls_data[k].append(inv_val)

                if n_execs > 1:
                    means = []
                    stds = []
                    for k in keys:
                        means.append(np.mean(inv_ls_data[k]))
                        stds.append(np.std(inv_ls_data[k]))
                else:
                    vals = [inv_ls_data[k][0] for k in keys]

                display_labels = []
                for k in keys:
                    label = k

                    if 'R' in k:
                        label += " (k2/k1, escala log.)"
                    elif any(p in k for p in ['k1', 'g_slow', 'g_fast']):
                        label += " (escala log.)"

                    display_labels.append(label)

                setup_plot('Importancia en la puntuación de los parámetros' if n_execs == 1 else 'Promedio de la importancia en la puntuación de los parámetros',
                           'Parámetro', 'Importancia (1 / ARD kernel lengthscale)')
                # Crea un gráfico de barras con la importancia media de cada parámetro.
                if n_execs > 1:
                    plt.bar(display_labels, means, yerr=stds, color='dodgerblue',
                            linewidth=4/len(keys), capsize=30/len(keys))
                else:
                    plt.bar(display_labels, vals, color='dodgerblue',
                            linewidth=4/len(keys))
                plt.grid(True, axis='y', linewidth=0.2)
                plt.gca().set_axisbelow(True)
                plt.xticks(rotation=45, ha='right')
                save_plot(args.output_ard_lss)
            else:
                print(
                    "Advertencia: No se encontraron datos de lengthscales (ARD_lss) para generar la gráfica 4.")
        else:
            print(
                "Advertencia: No se encontraron datos de lengthscales (ARD_lss) para generar la gráfica 4.")
else:
    print("Advertencia: no se encontraron ejecuciones en los datos para generar las gráficas de puntuaciones.")

if args.output_avg_execution_time:
    # Calcula e imprime en un archivo de texto la media y desviación del tiempo de optimización
    if len(opt_times_list) > 0:
        mean_time = np.mean(opt_times_list)
        std_time = np.std(opt_times_list)
        with open(args.output_avg_execution_time, 'w') as f:
            f.write(f"Mean: {mean_time}\n")
            f.write(f"Standard deviation: {std_time}\n")
    else:
        print("Advertencia: No se encontraron datos de tiempo de optimización (optimization_time).")

if args.output_avg_best_params:
    # Calcula y guarda en CSV la media y desviación de los mejores parámetros (best_params)
    if len(best_params_list) > 0:
        flat_params_list = []
        # Aplanar el diccionario para soportar tanto la versión offline (plano) como online (anidado por dirección)
        for params in best_params_list:
            flat_params = {}
            for k, v in params.items():
                if isinstance(v, dict):
                    for sub_k, sub_v in v.items():
                        flat_params[f"{k} {sub_k}"] = sub_v
                else:
                    flat_params[k] = v
            flat_params_list.append(flat_params)

        keys = sorted(list(flat_params_list[0].keys()))

        if len(keys) > 1:
            param_values = {k: [] for k in keys}
            for params in flat_params_list:
                for k in keys:
                    param_values[k].append(params[k])

            with open(args.output_avg_best_params, 'w') as f:
                f.write("param,mean,std\n")
                for k in keys:
                    vals = param_values[k]
                    if vals:
                        mean_val = np.mean(vals)
                        std_val = np.std(vals)
                        f.write(f"{k},{mean_val},{std_val}\n")
        else:
            print(
                "Advertencia: No se encontraron datos de los mejores parámetros (best_params).")
    else:
        print(
            "Advertencia: No se encontraron datos de los mejores parámetros (best_params).")
