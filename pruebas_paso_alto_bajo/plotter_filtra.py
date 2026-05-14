import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(csv_file):
    try:
        # Cargar el CSV generado por C++
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error al leer el archivo: {e}")
        return

    # Extraer columnas
    t = df['t']
    original = df['original']
    onda_lenta = df['onda_lenta']
    onda_rapida = df['onda_rapida']

    # Crear la figura (replicando tu estilo)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    
    # Gráfica 1: Original
    ax1.plot(t, original, color='black', linewidth=1)
    ax1.set_title("Señal Original (Leída por C++)")
    ax1.set_ylabel("Amplitud")
    ax1.grid(True)

    # Gráfica 2: Onda Lenta (LPF Zero-Phase)
    ax2.plot(t, onda_lenta, color='red', linewidth=1.5)
    ax2.set_title("Onda Lenta (C++ IIR1: Forward-Backward)")
    ax2.set_ylabel("Amplitud")
    ax2.grid(True)

    # Gráfica 3: Onda Rápida (Complementaria)
    ax3.plot(t, onda_rapida, color='blue', linewidth=1.5)
    ax3.set_title("Onda Rápida (Original - Onda Lenta)")
    ax3.set_xlabel("Tiempo (ms)")
    ax3.set_ylabel("Amplitud")
    ax3.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Si no pasas nombre de archivo, usa 'results.csv' por defecto
    file_to_plot = sys.argv[1] if len(sys.argv) > 1 else 'results.csv'
    main(file_to_plot)