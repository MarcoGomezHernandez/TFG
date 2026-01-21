import numpy as np
import matplotlib.pyplot as plt

# Parámetros globales del sistema (macros)
E = 3.0
MU = 0.0021
S = 4.0

# Parámetros de discretización
DELTA_T = 0.001

# Condiciones iniciales
X0 = -1.6
Y0 = -10.0
Z0 = 0.0

# Factor de diezmado (downsampling)
K = 10


def simulate_hindmarsh_rose(filename, t_total, t_offset, x0=X0, y0=Y0, z0=Z0):
    """
    Simula el sistema de Hindmarsh-Rose usando el método de Euler hacia adelante.
    Guarda los datos directamente en un archivo.

    Args:
        filename: Nombre del archivo donde guardar los datos
        t_total: Duración del intervalo temporal a guardar (después del offset)
        x0, y0, z0: Condiciones iniciales del sistema
        t_offset: Tiempo inicial desde el cual comenzar a guardar datos
    """
    # Tiempo total de simulación (offset + tiempo a graficar)
    t_sim = t_offset + t_total

    # Número total de iteraciones
    n_iterations = int(t_sim / DELTA_T)

    # Iteración a partir de la cual empezar a guardar
    n_offset = int(t_offset / DELTA_T)

    # Condiciones iniciales
    x_n, y_n, z_n = x0, y0, z0

    # Abrir archivo para escritura
    with open(filename, 'w') as f:
        # Escribir encabezado
        f.write("t,f\n")

        # Bucle principal de integración
        for n in range(1, n_iterations + 1):
            # Ecuaciones de recurrencia (Euler hacia adelante)
            x_n1 = x_n + DELTA_T * (y_n + 3*(x_n**2) - (x_n**3) - z_n + E)
            y_n1 = y_n + DELTA_T * (1 - 5*(x_n**2) - y_n)
            z_n1 = z_n + DELTA_T * MU * (-z_n + S * (x_n + 1.6))

            # Actualizar estado
            x_n, y_n, z_n = x_n1, y_n1, z_n1

            # Guardar cada K iteraciones, solo después del offset
            if n >= n_offset and n % K == 0:
                t_current = n * DELTA_T
                f.write(f"{t_current},{x_n}\n")


def plot_from_file(filename, title, xlabel='Tiempo (t)', ylabel='f(t)'):
    """
    Representa gráficamente datos de un archivo CSV con columnas t,f.

    Args:
        filename: Nombre del archivo a leer
        title: Título de la gráfica
        xlabel: Etiqueta del eje X
        ylabel: Etiqueta del eje Y
    """
    print(f"Leyendo datos desde {filename}...")

    # Leer datos del archivo
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    t_array = data[:, 0]
    x_array = data[:, 1]

    print(f"Datos cargados. Puntos: {len(t_array)}")

    # Crear figura
    plt.figure(figsize=(12, 6))
    plt.plot(t_array, x_array, linewidth=0.5, color='orange')
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def main():
    """
    Función principal para ejecutar la simulación y visualización.
    """
    filename = "hindmarsh_rose_data.csv"

    print(f"Simulando sistema de Hindmarsh-Rose...")
    print(f"Parámetros: e={E}, μ={MU}, S={S}, Δt={DELTA_T}")
    print(f"Condiciones iniciales: x₀={X0}, y₀={Y0}, z₀={Z0}")
    print(f"Factor de diezmado: k={K}")
    print(f"Guardando datos en: {filename}")

    # Simular y guardar en archivo
    simulate_hindmarsh_rose(filename, t_total=4000.0, t_offset=1000.0)

    print("Simulación completada.")

    # Graficar desde archivo
    plot_from_file(
        filename,
        title='Dinámica del modelo de Hindmarsh-Rose',
        ylabel='x(t) - Potencial de membrana'
    )


if __name__ == "__main__":
    main()
