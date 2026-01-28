# Modelo de Sinapsis Química (ChemicalSynapsis)

Basado en los archivos proporcionados (`ChemicalSynapsis.h` y `ChemicalSynapsisModel.h`), este modelo simula una sinapsis química compuesta por dos componentes: una **transmisión rápida** (instantánea) y una **transmisión lenta** (que depende del tiempo y la acumulación de un neurotransmisor o modulador).

El modelo calcula la corriente total $I_{syn}$ que entra a la neurona postsináptica.

---

## 1. Fórmulas Completas del Modelo

La corriente sináptica total $I_{syn}$ es la suma de dos corrientes:

$$I_{syn} = I_{fast} + I_{slow}$$

### A. Componente Rápida ($I_{fast}$)
Es una función sigmoidea instantánea del voltaje de la neurona presináptica ($V_{pre}$). No requiere ecuaciones diferenciales.

$$I_{fast} = \frac{g_{fast} \cdot (V_{post} - E_{syn})}{1 + \exp(s_{fast} \cdot (V_{fast} - V_{pre}))}$$

### B. Componente Lenta ($I_{slow}$)
Depende de una variable de estado interna $m_{slow}$ (fracción de canales abiertos o concentración activa), que evoluciona en el tiempo mediante una ecuación diferencial.

$$I_{slow} = g_{slow} \cdot m_{slow} \cdot (V_{post} - E_{syn})$$

La dinámica de la variable $m_{slow}$ está definida por la siguiente **ecuación diferencial ordinaria (EDO)**:

$$\frac{dm_{slow}}{dt} = \left[ \frac{k_1 \cdot (1 - m_{slow})}{1 + \exp(s_{slow} \cdot (V_{slow} - V_{pre}))} \right] - k_2 \cdot m_{slow}$$

**Donde:**
* El primer término representa la **activación (producción)** dependiente del voltaje presináptico.
* El segundo término ($-k_2 \cdot m_{slow}$) representa la **relajación o degradación** del transmisor.

---

## 2. Discretización por el Método de Euler

Para implementar esto en una simulación paso a paso (computacional), convertimos la ecuación diferencial en una ecuación algebraica usando el paso de tiempo $h$ ($\Delta t$).

Dada la ecuación: $\frac{dm_{slow}}{dt} = f(V_{pre}, m_{slow})$

La actualización para el siguiente instante de tiempo ($t + h$) es:

$$m_{slow}^{(t+h)} = m_{slow}^{(t)} + h \cdot \left( \left[ \frac{k_1 \cdot (1 - m_{slow}^{(t)})}{1 + \exp(s_{slow} \cdot (V_{slow} - V_{pre}))} \right] - k_2 \cdot m_{slow}^{(t)} \right)$$

---

## 3. Explicación Técnica de los Parámetros

A continuación se detalla cada parámetro y el efecto físico de su variación en la transmisión sináptica.

### Parámetros Generales
* **Esyn (Potencial de Inversión):** Define el voltaje al que la corriente neta es cero. Determina la naturaleza de la sinapsis. El $E_{syn}$ actúa como un "atractor". Físicamente, la corriente sináptica siempre intenta llevar el potencial de la membrana ($V_{post}$) hacia el valor de $E_{syn}$
    * **Efecto:** Si $E_{syn}$ es mayor que el potencial de reposo (ej. $0\text{ mV}$), la sinapsis es **excitatoria**. Si es menor (ej. $-80\text{ mV}$), es **inhibitoria**.

### Parámetros de la Componente Rápida
* **gfast (Conductancia Máxima Rápida):** Es la "fuerza" máxima de la parte rápida de la sinapsis.
    * **Si crece:** Aumenta la magnitud de la corriente $I_{fast}$ inyectada instantáneamente tras el disparo.
    * **Si decrece:** La respuesta rápida se debilita.
* **Vfast (Umbral de Activación Rápida):** El voltaje presináptico ($V_{pre}$) en el cual la conductancia alcanza el 50% de su máximo.
    * **Si crece:** La sinapsis se vuelve menos sensible; requiere que la neurona presináptica esté más despolarizada para activarse.
    * **Si decrece:** Se activa con voltajes más bajos.
* **sfast (Pendiente Rápida):** Controla la abrupticidad de la activación.
    * **Si crece:** La transición es más abrupta (tipo interruptor).
    * **Si decrece:** La activación es más gradual y suave.

### Parámetros de la Componente Lenta (Cinética)
* **gslow (Conductancia Máxima Lenta):** La "fuerza" máxima de la parte lenta.
    * **Si crece:** Aumenta la influencia de la corriente sostenida $I_{slow}$.
* **Vslow (Umbral de Activación Lenta):** Voltaje presináptico necesario para activar la producción de $m_{slow}$.
* **sslow (Pendiente Lenta):** Sensibilidad de la activación de la variable lenta respecto al voltaje.
* **k1 (Tasa de Activación):** Velocidad a la que crece $m_{slow}$ cuando hay estímulo.
    * **Si crece:** La corriente lenta alcanza su máximo más rápido (menor latencia).
* **k2 (Tasa de Degradación):** Velocidad a la que desaparece $m_{slow}$ en ausencia de estímulo.
    * **Si crece:** La sinapsis lenta se "apaga" muy rápido.
    * **Si decrece:** El efecto sináptico perdura mucho más tiempo (memoria o integración temporal).

---










# Diferencia de algunos params confusos

| Parámetro | Nombre Técnico | ¿Qué controla? (Función) | Dominio | Efecto Principal al Aumentar |
| :--- | :--- | :--- | :--- | :--- |
| **g** | **Conductancia** | La **Magnitud** (Volumen). Escala la fuerza de la corriente sináptica. | Amplitud | La corriente es más fuerte (mayor amplitud) para el mismo estímulo. |
| **Esyn** | **Potencial de Inversión** | La **Dirección** y Límite. Define si excita o inhibe y hasta qué voltaje puede empujar. | Voltaje (Límite) | Cambia el punto de equilibrio. Si es muy alto, inhibe; si es muy bajo, excita (dependiendo del $V$ actual). |
| **s** | **Pendiente (Slope)** | La **Selectividad**. Define qué tan abrupta es la activación ante cambios de voltaje. | Voltaje (Forma) | Se vuelve "todo o nada". Ignora voltajes bajos y se activa de golpe al cruzar el umbral. |
| **k1** | **Tasa de Activación** | La **Velocidad**. Define la rapidez con la que se genera la corriente una vez activada. | Tiempo (Dinámica) | Reacción instantánea. Elimina el retraso (latencia) entre el disparo y la respuesta. |









# Pruebas con Neuronas Hindmarsh-Rose (HR)

Para realizar pruebas con neuronas **Hindmarsh-Rose (HR)**, es crucial entender que este modelo es **adimensional**. A diferencia de modelos biológicos reales, la variable de voltaje en HR ($x$) oscila típicamente entre **-1.6 (reposo)** y **+2.0 (pico del disparo)**.

### 1. Parámetros de Intensidad (gfast, gslow)
Controlan la fuerza con la que la sinapsis afecta a la neurona destino.
* **gfast (Conductancia Rápida)**
    * **Demasiado Bajo (0.001):** La sinapsis es imperceptible.
    * **Adecuado (0.1 - 0.5):** Suficiente para provocar potenciales postsinápticos visibles o inducir sincronización.
    * **Demasiado Alto (10.0):** Efecto de "bloqueo" (*clamping*). La corriente domina la ecuación, forzando a la neurona al valor de $E_{syn}$.

### 2. Parámetros de Voltaje (Esyn, Vfast, Vslow)
* **Esyn (Potencial de Inversión)**
    * **Nota Crítica:** En esta librería, con conductancia positiva, un $E_{syn}$ menor que el voltaje actual excita, y uno mayor inhibe.
    * **Adecuado (-2.0 a 2.0):**
        * Use **-2.0** para una sinapsis **Excitatoria**.
        * Use **+2.0** para una sinapsis **Inhibitoria**.
* **Vfast y Vslow (Umbral de Activación)**
    * **Demasiado Bajo (-2.0):** La sinapsis está siempre activa (ruido constante).
    * **Adecuado (0.0):** Justo en la mitad del disparo de la neurona HR.
    * **Demasiado Alto (+2.5):** La sinapsis nunca se activa.

### 3. Parámetros de Forma (sfast, sslow)
* **sfast y sslow (Pendiente de la Sigmoide)**
    * **Adecuado (5.0 - 10.0):** Comportamiento de "interruptor" biológico.
    * **Demasiado Alto (100.0):** Funciona como un escalón binario perfecto; puede causar inestabilidad numérica en el integrador.

### 4. Parámetros Temporales (k1, k2)
Controlan la "memoria" o duración de la componente lenta.
* **k1 (Velocidad de Activación)**
    * **Adecuado (1.0 - 2.0):** La variable interna sube rápidamente durante el disparo.
* **k2 (Velocidad de Degradación / Olvido)**
    * **Adecuado (0.01 - 0.05):** El efecto dura unos cuantos milisegundos (unidades HR) después del disparo.
    * **Demasiado Alto (2.0):** La componente lenta se vuelve indistinguible de la rápida.