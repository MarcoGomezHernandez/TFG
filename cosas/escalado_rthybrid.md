Aquí tienes todo el código solicitado, incluyendo la función `signal_period` y especificando claramente el archivo de origen para cada bloque. Está organizado siguiendo el flujo lógico de ejecución: desde la configuración inicial, pasando por el análisis de la señal y el cálculo de factores, hasta su aplicación en el bucle de tiempo real.

### 1. Inicialización de Parámetros

**Archivo:** `clamp/src/clamp.c`
Antes de iniciar el hilo de tiempo real, se preparan los parámetros base (`freq`, `period`, etc.) a partir de la configuración del usuario.

```c
/* ... (Dentro de la función clamp) ... */

    /* Configuración de argumentos para el hilo de tiempo real (r_args) */
    
    r_args.time_var = args->time_var;           // Duración de la interacción
    r_args.before = args->before;               // Tiempo de control antes
    r_args.after = args->after;                 // Tiempo de control después
    r_args.observation = args->observation;     // Tiempo de observación para calibración inicial
    
    // Cálculo del periodo en nanosegundos a partir de la frecuencia en Hz (definida por usuario)
    r_args.period =  (1 / args->freq) * NSEC_PER_SEC; 
    r_args.freq = args->freq;                   // Frecuencia de muestreo (Hz)
    
    r_args.filename = args->filename;           // Nombre del archivo de salida
    r_args.input_factor = args->input_factor;   // Factor de conversión del hardware de entrada
    r_args.output_factor = args->output_factor; // Factor de conversión del hardware de salida
    
    r_args.sec_per_burst = args->sec_per_burst; // Duración forzada de ráfaga (opcional)
    r_args.check_drift = args->check_drift;     // Flag para activar corrección de deriva

```

### 2. Estructuras de Datos

**Archivo:** `clamp/includes/types_clamp.h`
Estructura utilizada para pasar punteros a las variables "vivas" del bucle de control hacia la función de corrección de deriva.

```c
typedef struct {
    double * scale_virtual_to_real;
    double * scale_real_to_virtual;
    double * offset_virtual_to_real;
    double * offset_real_to_virtual;
    double * max_window;
    double * min_window;
    double * max_rel_real;
    double * min_rel_real;
    double max_abs_model;
    double min_abs_model;
    synapse_model * sm_model_to_live;
    synapse_model * sm_live_to_model;
    synapse_model * sm_live_to_model_scaled;
} fix_drift_args;

```

### 3. Análisis Inicial de la Señal Viva

**Archivo:** `clamp/src/calibrate_functions_phase1.c`
Funciones para "escuchar" la neurona biológica, determinar sus rangos de voltaje y calcular su periodo de disparo (`signal_period`).

```c
int ini_recibido (double *min_rel_real, double *min_abs_real, double *max_abs_real, double *max_rel_real, double *period_signal, Daq_session * session, int chan, int period, int freq, char* filename, double input_factor, int observation_time){

    /*VARIABLES TO DETERMINATE RANGES*/
    int i=0;
    double retval=0.0;
    struct timespec ts_target, ts_start;
    double max_abs = -DBL_MAX;
    double min_abs =  DBL_MAX;
    double percentage_min = 0.10;
    double percentage_max = 0.90;
    double range;

    /*DAQ Config*/
    int n_channels = 1;
    int in_channels [1];
    double ret_values [1];
    in_channels[0] = chan;

    /*RT config*/
    clock_gettime(CLOCK_MONOTONIC, &ts_target);
    ts_assign (&ts_start,  ts_target);
    ts_add_time(&ts_target, 0, period);

    /*DECLARACIONES DE ARRAYS Y SUS TAMAÑOS*/
    int size_signal = freq*observation_time;
    double * signal = (double*) malloc(sizeof(double) * size_signal);

    for (i=0; i<size_signal; i++){

        /*SLEEP & READ DATA*/
        clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ts_target, NULL);
        if (daq_read(session, n_channels, in_channels, ret_values) != 0) {
            return -1;
        }
        retval = (ret_values[0] * 1000.0) / input_factor;
        signal[i] = retval;
        
        /*DETECT MAX/MIN*/
        if(retval>max_abs){
            max_abs=retval;
        }else if(retval<min_abs){
            min_abs=retval;
        }
        
        /*NEXT PERIOD*/
        ts_add_time(&ts_target, 0, period); 
    }
    
    /*RETURN*/
    range = max_abs - min_abs;
    *min_abs_real = min_abs;
    *max_abs_real = max_abs;
    *min_rel_real = percentage_min * range + min_abs;
    *max_rel_real = percentage_max * range + min_abs;

    /*SIGNAL PERIOD*/
    *period_signal = signal_period (observation_time, signal, size_signal, *max_rel_real, *min_rel_real);

    return OK;
}

/***
Get period using two thresholds
***/
double signal_period(int seg_observacion, double * signal, int size, double th_up, double th_on){
    int up=FALSE;
    if (signal[0]>th_up)
        up=TRUE;

    double changes=0;
    int i=0;
    for (i=0; i<size; i++){
        if(up==FALSE && signal[i]>th_up){
            //Cambio de tendencia
            changes++;
            up=TRUE;
        }else if(up==TRUE && signal[i]<th_on){
            up=FALSE;
        }
    }
    double period = 1.0 / (changes/seg_observacion);
    return period;
}

```

### 4. Matemáticas de Escalado y Corrección de Deriva

**Archivo:** `clamp/src/calibrate_functions_phase2_a.c`
Incluye `calcula_escala` (factores lineales) y `fix_drift` (recalculo en tiempo real).

```c
/* Ordenar los ficheros de calibraciones */
void calcula_escala (double min_virtual, double max_virtual, double min_viva, double max_viva, double *scale_virtual_to_real, double *scale_real_to_virtual, double *offset_virtual_to_real, double *offset_real_to_virtual){
    
    double rg_virtual, rg_viva;
    
    rg_virtual = max_virtual-min_virtual;
    rg_viva = max_viva-min_viva;
    
    *scale_virtual_to_real = rg_viva / rg_virtual;
    *scale_real_to_virtual = rg_virtual / rg_viva;
    
    *offset_virtual_to_real = min_viva - (min_virtual*(*scale_virtual_to_real));
    *offset_real_to_virtual = min_virtual - (min_viva*(*scale_real_to_virtual));

    return;
}


void fix_drift (fix_drift_args args) {
    double per_min = 0.1, per_max = 0.1;

    calcula_escala (args.min_abs_model, args.max_abs_model, *(args.min_window), *(args.max_window), args.scale_virtual_to_real, args.scale_real_to_virtual, args.offset_virtual_to_real, args.offset_real_to_virtual);

    args.sm_live_to_model->offset = *(args.offset_real_to_virtual);
    args.sm_live_to_model->scale = *(args.scale_real_to_virtual);
    args.sm_live_to_model->min = *(args.min_window);
    args.sm_live_to_model->max = *(args.max_window);

    args.sm_live_to_model_scaled->offset = *(args.offset_real_to_virtual);
    args.sm_live_to_model_scaled->scale = *(args.scale_real_to_virtual);
    args.sm_live_to_model_scaled->min = *(args.min_window);
    args.sm_live_to_model_scaled->max = *(args.max_window);

    args.sm_model_to_live->offset = *(args.offset_virtual_to_real);
    args.sm_model_to_live->scale = *(args.scale_virtual_to_real);

    if(*(args.min_window) > 0){
        *(args.min_rel_real) = *(args.min_window) + (*(args.min_window) * per_min);
    }else{
        *(args.min_rel_real) = *(args.min_window) - (*(args.min_window) * per_min);
    }

    if(*(args.max_window) > 0){
        *(args.max_rel_real) = *(args.max_window) - (*(args.max_window) * per_max);
    }else{
        *(args.max_rel_real) = *(args.max_window) + (*(args.max_window) * per_max);
    }

    return;
}

```

### 5. Configuración del Modelo Neuronal (Hindmarsh-Rose)

**Archivo:** `model_library/neuron/Hindmarsh_Rose_1986/nm_hindmarsh_rose_1986.c`
Definición de límites intrínsecos y tablas precomputadas de escalado temporal.

```c
void nm_hindmarsh_rose_1986_init (neuron_model * nm, double * vars, double * params) {
	nm->dim = 3;
	nm->vars = (double *) malloc (sizeof(double) * nm->dim);
	copy_1d_array(vars, nm->vars, nm->dim);

	nm->n_params = 10;
	nm->params = (double *) malloc (sizeof(double) * nm->n_params);
	copy_1d_array(params, nm->params, nm->n_params);

	nm->max = 1.797032; // Rango intrínseco Máximo
	nm->min = -1.608734; // Rango intrínseco Mínimo
	nm->pts_burst = -1.0;

	nm->func = &nm_hindmarsh_rose_1986;
	nm->set_pts_burst = &nm_hindmarsh_rose_1986_set_pts_burst;
	nm->method = integration_method_selector(params[NM_HINDMARSH_ROSE_1986_DT]);

	return;
}

double nm_hindmarsh_rose_1986_set_pts_burst (double pts_live, neuron_model * nm) {
	int length = 0;
	int method = nm->params[NM_HINDMARSH_ROSE_1986_DT];

	switch(method) {
		case RK4:
		{
			length = 144.000000;
			double dts[] = {0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500, 0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500, 0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300, 0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100, 0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600, 0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600, 0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900, 0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500, 0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000, 0.029400, 0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200, 0.032600, 0.033000, 0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900, 0.036400, 0.036900, 0.037400, 0.038000, 0.038600, 0.039200, 0.039800, 0.040400, 0.041000, 0.041700, 0.042400, 0.043100, 0.043800, 0.044500, 0.045300, 0.046100, 0.046900, 0.047700, 0.048600, 0.049500, 0.050400, 0.051400, 0.052400, 0.053400, 0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600, 0.062000, 0.063400, 0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600, 0.077700, 0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000};
			double pts[] = {577638.000000, 481366.000000, 412599.000000, 357615.500000, 317880.000000, 286092.500000, 259143.333333, 237548.000000, 218869.500000, 203236.000000, 189687.000000, 177634.000000, 157897.000000, 142001.833333, 129024.142857, 113496.125000, 101304.555556, 97811.222222, 94527.400000, 91478.200000, 88619.400000, 85916.636364, 83389.636364, 81007.090909, 78743.583333, 76615.416667, 74599.250000, 72676.000000, 70859.076923, 69130.846154, 67476.642857, 65907.357143, 64402.666667, 62971.466667, 61602.533333, 60286.187500, 59030.250000, 57825.562500, 56664.411765, 55553.294118, 54485.000000, 52463.222222, 50586.263158, 48841.842105, 47211.050000, 45685.666667, 44255.818182, 42914.772727, 41650.739130, 40459.083333, 39335.208333, 38270.680000, 36778.346154, 35398.000000, 34117.571429, 32926.517241, 31815.833333, 30777.612903, 29493.939394, 28313.588235, 27223.638889, 25974.405405, 24834.410256, 23790.268293, 22647.767442, 21609.977778, 20513.166667, 19388.627451, 18381.132075, 17365.719298, 16361.600000, 15299.937500, 14223.202899, 13164.400000, 12147.123457, 11098.876404, 10071.693878, 9965.282828, 9861.100000, 9759.059406, 9626.242718, 9497.009615, 9371.179245, 9248.672897, 9129.293578, 9012.981818, 8899.594595, 8789.000000, 8681.149123, 8575.896552, 8473.170940, 8348.176471, 8226.809917, 8108.934426, 7994.379032, 7883.015873, 7774.710938, 7669.348837, 7566.801527, 7447.308271, 7331.525926, 7219.291971, 7110.435714, 7004.816901, 6902.298611, 6786.417808, 6674.355705, 6565.940397, 6460.987013, 6359.339744, 6247.012579, 6138.592593, 6033.866667, 5932.660714, 5822.777778, 5716.896552, 5614.796610, 5505.541436, 5400.467391, 5299.324468, 5192.348958, 5089.615385, 4982.075000, 4878.985294, 4772.014354, 4669.633803, 4564.178899, 4463.381166, 4360.214912, 4255.294872, 4149.216667, 4048.296748, 3946.654762, 3844.764479, 3743.041353, 3641.872263, 3541.583630, 3438.300000, 3336.926421, 3233.951299, 3133.666667, 3032.899696, 2932.320588, 2829.684659};

			select_dt_neuron_model(dts, pts, length, pts_live, &(nm->params[NM_HINDMARSH_ROSE_1986_DT]), &(nm->pts_burst));
			break;
		}
	}

	return nm->params[NM_HINDMARSH_ROSE_1986_DT];
}

```

### 6. Lógica de Selección de DT (Drift Check Temporal)

**Archivo:** `model_library/neuron/neuron_models_aux_functions.c`
Función auxiliar que busca un múltiplo entero de la duración de la ráfaga biológica.

```c
void select_dt_neuron_model (double * dts, double * pts, unsigned int length, double pts_live, double * dt, double * pts_burst) {
    double aux = pts_live;
    double factor = 1;
    double intpart, fractpart;
    int flag = 0;
    int i;

    *dt = -1;
    *pts_burst = -1;

    while (aux < pts[0]) {
        aux = pts_live * factor;
        factor += 1;

        for (i = length - 1; i >= 0; i--) {
        	if (pts[i] > aux) {
        		*dt = dts[i];
	            *pts_burst = pts[i];

	            fractpart = modf(*pts_burst / pts_live, &intpart);

	            if (fractpart <= 0.1*intpart) flag = 1;

	            break;
        	}
            
        }

        if (flag == 1) break;
    }

    if (flag == 0) {
        for (i = length - 1; i >= 0; i--) {
        	if (pts[i] > aux) {
        		*dt = dts[i];
	            *pts_burst = pts[i];

	            break;	
        	}
        }
    }

    return;
}

```

### 7. Integración en el Bucle de Tiempo Real (RT Thread)

**Archivo:** `clamp/src/rt_thread_functions.c`
Bucle principal que integra todas las piezas anteriores, incluyendo la adquisición de datos, el cálculo de `s_points` con su protección y la ejecución periódica del chequeo de deriva.

```c
void * rt_thread(void * arg) {

    /* ... Declaración de variables ... */
    
    /* Inicialización de rangos del modelo */
    min_abs_model = args->nm.min;
    max_abs_model = args->nm.max;

    if (args->n_in_chan > 0) {

        /* Análisis inicial de la señal viva para obtener rangos */
        if ( ini_recibido (&min_rel_real, &min_abs_real, &max_abs_real, &max_rel_real, &external_firing_rate, session, calib_chan, args->period, args->freq, args->filename, args->input_factor, args->observation) == -1 ) {
            // Manejo de error y limpieza...
            /* ... */
            pthread_exit(NULL);
        }

        /* Si el usuario fijó la duración de ráfaga, se usa ese valor */
        if (args->sec_per_burst != -1) external_firing_rate = args->sec_per_burst;

        /* Cálculo inicial de factores de escala */
        calcula_escala (min_abs_model, max_abs_model, min_abs_real, max_abs_real, &scale_virtual_to_real, &scale_real_to_virtual, &offset_virtual_to_real, &offset_real_to_virtual);
    } 

    /* Cálculo de resolución temporal del modelo */
    external_pts_per_burst = args->freq * external_firing_rate;
    int_params.dt = args->nm.set_pts_burst(external_pts_per_burst, &(args->nm));
    
    /* Configuración del método de integración */
    if (int_params.dt != -1) {
        int_params.method = args->nm.method;
    } else {
        int_params.method = NULL;
    }

    /* Cálculo de relación de puntos modelo/real */
    s_points = args->nm.pts_burst / external_pts_per_burst;

    /* PROTECCIÓN: Asegura que la relación de puntos sea al menos 1:1 */
    if (s_points == 0) s_points = 1;

    /* ... Inicialización de estructuras de calibración y sinapsis ... */

    /* Bucle de Experimento */
    for (i = 0; i < n_loops; i++) {
        experiment_loop(&(lp[i]), s_points);
    }

    /* ... Limpieza y cierre ... */
}


void experiment_loop (struct Loop_params * lp, int s_points) {
    /* ... Variables locales ... */
    
	for (i = 0; i < loop_points || lp->infinite == TRUE; i++) {

        if (i % s_points == 0) {
            
            /* ... Sincronización de tiempo (clock_nanosleep) ... */

            /* Envío de voltaje escalado del modelo al DAQ */
            v_model_scaled = args->nm.vars[0] * scale_virtual_to_real + offset_virtual_to_real;
            /* ... Escritura en DAQ ... */

            /* Lectura del DAQ */
            if (daq_read(session, args->n_in_chan, args->in_channels, input_values) != 0) {
                /* ... Manejo de error ... */
            }

            /* Aplicación del factor de entrada (mV) */
            if (args->n_in_chan > 0) input_values[0] = (input_values[0] * 1000.0) / args->input_factor;


            /* Lógica de corrección de deriva (Drift Check & Fix) */
            if (args->check_drift == TRUE && args->sm_live_to_model.type != SM_EMPTY) {
                // Actualización de ventana deslizante (punto a punto)
                drift_aux_range = args->sm_live_to_model.max - args->sm_live_to_model.min;
                if ((min_window > input_values[0]) && (input_values[0] > (args->sm_live_to_model.min - drift_aux_range))) min_window = input_values[0];
                if ((max_window < input_values[0]) && (input_values[0] < (args->sm_live_to_model.max + drift_aux_range))) max_window = input_values[0];

                // Corrección periódica (cada 'drift_n_burst' ráfagas; definido constante a 2)
                if (drift_counter >= (drift_n_burst * external_pts_per_burst) && max_window != -999999 && min_window != 999999) {
                    drift_counter = 0;

                    // Preparación de argumentos para fix_drift (punteros a variables vivas)
                    fx_args.scale_virtual_to_real = &scale_virtual_to_real;
                    fx_args.scale_real_to_virtual = &scale_real_to_virtual;
                    fx_args.offset_virtual_to_real = &offset_virtual_to_real;
                    fx_args.offset_real_to_virtual = &offset_real_to_virtual;
                    fx_args.max_window = &max_window;
                    fx_args.min_window = &min_window;
                    fx_args.max_rel_real = &max_rel_real;
                    fx_args.min_rel_real = &min_rel_real;
                    fx_args.max_abs_model = max_abs_model;
                    fx_args.min_abs_model = min_abs_model;
                    fx_args.sm_model_to_live = &(args->sm_model_to_live);
                    fx_args.sm_live_to_model = &(args->sm_live_to_model);
                    fx_args.sm_live_to_model_scaled = &(args->sm_live_to_model_scaled);

                    // Ejecución de la corrección
                    fix_drift(fx_args);

                    // Reseteo de ventana
                    max_window = -999999;
                    min_window = 999999;
                }

                drift_counter++;
            }
        }

        if (lp->interaction == TRUE) {
            /* Cálculo de sinapsis usando los factores (scale/offset configurados en el struct) */
            
            // Entrada de viva a modelo (escalada)
            args->sm_live_to_model_scaled.calibrate = SYN_CALIB_POST;
            args->sm_live_to_model_scaled.func(args->nm.vars[0], input_values[0], &(args->sm_live_to_model_scaled), &c_external_scaled);
            
            // Entrada de viva a modelo (para inyección de corriente)
            args->sm_live_to_model.calibrate = SYN_CALIB_PRE;
            args->sm_live_to_model.func(args->nm.vars[0], input_values[0], &(args->sm_live_to_model), &c_external);
        }

        /* Ejecución del modelo neuronal */
        args->nm.func(args->nm, c_external);

        if (lp->interaction == TRUE) {
            /* Sinapsis de modelo a viva */
            args->sm_model_to_live.calibrate = SYN_CALIB_PRE;
            args->sm_model_to_live.func(input_values[0], args->nm.vars[0], &(args->sm_model_to_live), &c_model);
        }
    }
}

```