# Para HR-HR, con extremos:

    # gfast_values = [0.001, 0.1, 0.3, 0.5, 10.0]
    # gslow_values = [0.001, 0.1, 0.3, 0.5, 10.0]
    # Esyn_values = [-5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0]
    # Vfast_values = [-2.0, -1.0, 0.0, 1.0, 2.0]
    # Vslow_values = [-2.0, -1.0, 0.0, 1.0, 2.0]
    # sfast_values = [0.5, 5.0, 7.5, 10.0, 100.0]
    # sslow_values = [0.5, 5.0, 7.5, 10.0, 100.0]
    # k1_values = [0.01, 1.0, 1.5, 2.0, 50.0]
    # k2_values = [0.001, 0.01, 0.03, 0.05, 2.0]

    # gfast_values = [0.3]
    # gslow_values = [0.3]
    # Esyn_values = [0.0]
    # Vfast_values = [0.0]
    # Vslow_values = [0.0]
    # sfast_values = [7.5]
    # sslow_values = [7.5]
    # k1_values = [1.5]
    # k2_values = [0.03]



    # Para HR-PD, con extremos:

    # gfast_values = [1e-12, 1e-9, 5e-8, 5e-7, 1e-3]
    # gslow_values = [1e-12, 1e-9, 5e-8, 5e-7, 1e-3]
    # Esyn_values = [-0.15, -0.08, -0.05, 0.0, 0.15]
    # Vfast_values = [-0.060, -0.050, -0.042, -0.035, -0.020]
    # Vslow_values = [-0.060, -0.050, -0.042, -0.035, -0.020]
    # sfast_values = [1.0, 50.0, 250.0, 1000.0, 2000.0]
    # sslow_values = [1.0, 50.0, 250.0, 1000.0, 2000.0]
    # k1_values = [1.0, 50.0, 500.0, 2000.0, 100000.0]
    # k2_values = [0.1, 10.0, 100.0, 200.0, 250.0]

    # gfast_values = [5e-8]
    # gslow_values = [5e-8]
    # Esyn_values = [-0.05]
    # Vfast_values = [-0.042]
    # Vslow_values = [-0.042]
    # sfast_values = [250.0]
    # sslow_values = [250.0]
    # k1_values = [500.0]
    # k2_values = [100.0]

    # best islow sin ir uno a uno, a lo bruto:
    # gfast_values = [5e-8] # da igual
    # gslow_values = [1e-9]
    # Esyn_values = [-0.15]
    # Vfast_values = [-0.042] # da igual
    # Vslow_values = [-0.035]
    # sfast_values = [250.0] # da igual
    # sslow_values = [1000.0]
    # k1_values = [1.0]
    # k2_values = [200.0]

    # bestislow llendo uno a uno (mejor puntuacion de todas):
    # gfast_values = [5e-8]
    # gslow_values = [1e-9]
    # Esyn_values = [-0.15]
    # Vfast_values = [-0.042]
    # Vslow_values = [-0.042]
    # sfast_values = [250.0]
    # sslow_values = [250.0]
    # k1_values = [100000.0]
    # k2_values = [10.0]

    # python plot.py yaml_configs/config_gf5e-08_gs5e-08_sf250.0_Vf-0.042_Es-0.05_ss250.0_Vs-0.042_k1500.0_k2100.0.asc