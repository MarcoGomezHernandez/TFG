# Estabilidad de Redes a partir de la Regulación Dependiente de la Actividad de las Conductancias Neuronales

## Explicación Completa del Artículo de Golowasch et al. (1999)

---

## 1. Introducción: El Problema Fundamental

El sistema nervioso enfrenta una paradoja fundamental:

- **Debe ser plástico**: cambiar con la experiencia, el aprendizaje y el desarrollo.
- **Debe ser estable**: mantener propiedades funcionales consistentes a pesar de perturbaciones ambientales e internas.

A nivel de neurona individual, esto es particularmente evidente: los canales iónicos se sintetizan, degradan e insertan constantemente en la membrana, y aun así la neurona mantiene durante años propiedades eléctricas relativamente estables.

Esto sugiere la existencia de **mecanismos de retroalimentación** que:
1. Miden el nivel de actividad neuronal (por ejemplo, a través del flujo de iones de calcio).
2. Ajustan las conductancias de canales iónicos para mantener una actividad "objetivo".

Trabajos anteriores demostraron que neuronas individuales pueden **autorregular** sus corrientes iónicas. El artículo de Golowasch y colegas plantea una pregunta novedosa: **¿Pueden estos mecanismos individuales también estabilizar la actividad de una red completa de neuronas?**

---

## 2. Anatomía Básica de la Neurona

### 2.1. Estructura Fundamental

Toda neurona tiene tres componentes principales:

1. **Soma** (cuerpo celular)
   - Contiene el núcleo y la mayor parte del citoplasma.
   - Es donde se integran las señales de entrada.

2. **Dendritas**
   - Ramificaciones que salen del soma.
   - Reciben información de otros axones.

3. **Axón**
   - "Cable" largo que sale del soma.
   - Conduce la señal **lejos** del cuerpo celular.
   - Termina en pequeñas ramas (terminales sinápticas) que contactan con otras neuronas.

### 2.2. Las Sinapsis

Una **sinapsis** es la conexión funcional entre dos neuronas:

- La **neurona presináptica** es aquella cuyo axón envía la señal.
- La **neurona postsináptica** es la que recibe la información.
- El contacto ocurre entre el axón terminal de la presináptica y el soma (o dendritas) de la postsináptica.

En este modelo, **todas las sinapsis son inhibidoras**, lo que significa que la neurona presináptica intenta hacer que la postsináptica sea menos activa.

### 2.3. Corrientes Iónicas

Dentro y fuera de la neurona hay iones disueltos (Na⁺, K⁺, Ca²⁺, Cl⁻). La membrana neuronal actúa como una barrera selectiva con **puertas (canales iónicos)** que se abren y cierran.

Cuando se abre un canal:
- Los iones fluyen a través de la membrana (esto es una "corriente iónica").
- Cambian el **potencial de membrana** (voltaje eléctrico).

Tipos de corrientes principales en este modelo:

- **Corriente de K⁺ (potasio)**: cuando fluye HACIA AFUERA, hiperpolariza (voltaje más negativo).
- **Corriente de Na⁺ (sodio)**: cuando fluye HACIA ADENTRO, despolariza (voltaje más positivo).
- **Corriente de Ca²⁺ (calcio)**: cuando fluye HACIA ADENTRO, despolariza.
- **Corriente de fuga**: pequeña corriente constante que mantiene el equilibrio.

La **suma total de todas las corrientes** determina si el voltaje de la membrana sube, baja o se mantiene estable.

---

## 3. El Ganglio Estomatogástrico: Preparación Biológica

### 3.1. ¿Qué es?

El **ganglio estomatogástrico (STG)** es un pequeño cúmulo de neuronas en crustáceos (cangrejos, langostas) que controla parte del aparato digestivo. Es famoso en neurociencia porque:

- Tiene **muy pocas neuronas** (algunas decenas).
- Es relativamente estable y fácil de estudiar en preparación aislada.
- Genera **ritmos motores** muy claros y reproducibles.

### 3.2. El Ritmo Pilórico

Uno de los ritmos más estudiados es el **ritmo pilórico**, que coordina músculos para mover el contenido del estómago. Es un patrón cíclico que se repite regularmente.

### 3.3. Neuronas Principales

El artículo se enfoca en tres tipos de neuronas:

| Tipo | Descripción | Rol en el Ritmo |
|------|-------------|-----------------|
| **AB/PD** | Anterior Burster + Pyloric Dilator (eléctricamente acopladas) | Inicia cada ciclo |
| **LP** | Lateral Pyloric | Segunda fase del ritmo |
| **PY** | Pyloric neurons (grupo de neuronas acopladas) | Tercera fase del ritmo |

El ritmo resultante es **tripásico**: AB/PD → LP → PY → AB/PD (y se repite).

### 3.4. El Fenómeno Clave: Neuromodulación

En condiciones normales (in vivo), el STG recibe sustancias neuromoduladoras a través de un nervio llamado **stn** (stomatogastric nerve). Estas sustancias son neuropéptidos como la **proctolina**.

**Observación experimental crucial**: 
- Con moduladores intactos: ritmo pilórico robusto a ~1 Hz.
- Si se bloquea el stn: ritmo cae casi a cero.
- Si se mantiene bloqueado ~24 h: **el ritmo reaparece lentamente**, pero más lento (~0.3-0.4 Hz), ahora **sin necesidad de moduladores**.

Esto sugiere que la red se **reorganiza** durante el bloqueo, cambiando sus propiedades intrínsecas para poder oscilar de forma autónoma.

---

## 4. Modelo Matemático: Estructura General

### 4.1. Dos Compartimentos por Neurona

El modelo representa cada neurona como si tuviera dos "subpartes eléctricas" conectadas:

1. **Compartimento somático** ($V_s$)
   - Representa soma + neurita principal.
   - Genera oscilaciones lentas, "plateaus", y la parte lenta del bursting.
   - Contiene: canales de fuga, Ca²⁺, K⁺, K⁺ tipo A, y proctolina (en AB/PD y LP).

2. **Compartimento axonal** ($V_a$)
   - Representa la zona de inicio del axón.
   - Genera los potenciales de acción rápidos.
   - Contiene: canales de fuga, Na⁺ rápido, y K⁺ retardada.

Ambos compartimentos están **acoplados eléctricamente** mediante una conductancia de conexión.

### 4.2. Conectividad Sináptica

Las tres neuronas están interconectadas en la siguiente configuración:

- **AB/PD** inhibe LP y PY (ambas con componentes rápida y lenta).
- **LP** inhibe AB/PD y PY (rápidas).
- **PY** inhibe LP y AB/PD (rápidas).

**Todas las sinapsis son inhibidoras**, y la dirección del flujo es siempre: **axón presináptico → soma postsináptico**.

---

## 5. Ecuaciones del Modelo

### 5.1. Compartimento Somático

La ecuación que governa la dinámica del voltaje somático es:

$$C_s \frac{dV_s}{dt} = -I_{syn} - \bar{g}_{Ls}(V_s - E_L) - \bar{g}_{Ca} m_{Ca}^3 h_{Ca} (V_s - E_{Ca}) - \bar{g}_K m_K^4 (V_s - E_K)$$
$$- \bar{g}_A m_A^3 h_A (V_s - E_K) - \bar{g}_{Proc} m_{proc} (V_s - E_{proc}) - \bar{g}_E (V_s - V_a)$$

**Interpretación de cada término**:

- $C_s \frac{dV_s}{dt}$: Capacitancia del soma (como un condensador eléctrico).

- $-I_{syn}$: Corriente sináptica total desde otras neuronas.

- $-\bar{g}_{Ls}(V_s - E_L)$: **Corriente de fuga**: mantiene el equilibrio de reposo.

- $-\bar{g}_{Ca} m_{Ca}^3 h_{Ca} (V_s - E_{Ca})$: **Corriente de Ca²⁺**:
  - $\bar{g}_{Ca}$: conductancia máxima (sujeta a regulación dependiente de la actividad).
  - $m_{Ca}, h_{Ca}$: variables de compuerta (activación e inactivación).
  - $E_{Ca} = 120$ mV (potencial de equilibrio del calcio, positivo).
  - Esta corriente suele ser **entrante** (despolariza la membrana).

- $-\bar{g}_K m_K^4 (V_s - E_K)$: **Corriente de K⁺ lenta**:
  - $E_K = -80$ mV (negativo).
  - Cuando la membrana está despolarizada, esta corriente es **saliente** (hiperpolariza).

- $-\bar{g}_A m_A^3 h_A (V_s - E_K)$: **Corriente de K⁺ tipo A** (transitoria): contribuye a controlar el inicio del bursting.

- $-\bar{g}_{Proc} m_{proc} (V_s - E_{proc})$: **Corriente moduladora de proctolina**:
  - Solo presente en AB/PD y LP.
  - $E_{proc} = -10$ mV (cerca del reposo).
  - Esta es la corriente que simula la neuromodulación externa.

- $-\bar{g}_E (V_s - V_a)$: **Corriente de acoplamiento soma-axón**.

### 5.2. Compartimento Axonal

La ecuación que governa el voltaje axonal es:

$$C_a \frac{dV_a}{dt} = -\bar{g}_{La} (V_a - E_L) - \bar{g}_{Na} m_{Na}^3 h_{Na} (V_a - E_{Na}) - \bar{g}_{Kd} m_{Kd}^4 (V_a - E_K) - \bar{g}_E (V_a - V_s)$$

**Interpretación**:

- $-\bar{g}_{Na} m_{Na}^3 h_{Na} (V_a - E_{Na})$: **Corriente de Na⁺ rápida**:
  - $E_{Na} = 20$ mV (positivo).
  - Cuando se abre, entra Na⁺ rápidamente y despolariza: genera el **potencial de acción**.

- $-\bar{g}_{Kd} m_{Kd}^4 (V_a - E_K)$: **Corriente de K⁺ retardada**:
  - Se activa con retraso respecto a Na⁺.
  - Repolariza la membrana y finaliza el potencial de acción.

Todos los parámetros de conductancia en el axón son **fijos** (no cambian con la actividad).

### 5.3. Variables de Compuerta (Hodgkin-Huxley)

Cada corriente iónica tiene variables de compuerta que describen la fracción de canales abiertos:

- $m(V,t)$: variable de **activación** (abre rápidamente).
- $h(V,t)$: variable de **inactivación** (se cierra lentamente).

Evolucionan según ecuaciones diferenciales de primer orden:

$$\tau_m(V)\frac{dm}{dt} = m_\infty(V) - m$$

$$\tau_h(V)\frac{dh}{dt} = h_\infty(V) - h$$

Donde $m_\infty(V)$ y $h_\infty(V)$ son **funciones sigmoides** (formas en S) que describen el valor de equilibrio a cada voltaje:

$$m_\infty(V) = \frac{1}{1 + \exp[s(V_{1/2} - V)]}$$

Y $\tau_m(V)$ y $\tau_h(V)$ son **constantes de tiempo** que determinan la velocidad de cambio:

$$\tau_m(V) = A + \frac{B}{1 + \exp[s(V_{1/2} - V)]}$$

**Interpretación cualitativa**:
- $m_\infty(V)$: valor "objetivo" de compuerta a cada voltaje (entre 0 y 1).
- $\tau_m(V)$: tiempo que tarda la variable en aproximarse a ese objetivo.
- Matemáticamente, cada variable de compuerta se relaja exponencialmente hacia su valor de equilibrio.

Los parámetros $V_{1/2}$ (voltaje de media-activación), $s$ (pendiente), $A$ y $B$ para cada canal están dados en la **Tabla 2** del artículo original.

### 5.4. El Mecanismo Central: Regulación Dependiente de la Actividad

**Este es el mecanismo más importante del artículo.**

Las conductancias máximas de Ca²⁺ y K⁺ en el soma **no son constantes**. En su lugar, dependen de una variable lenta $z$ que evoluciona según la actividad:

$$\bar{g}_{Ca} = \frac{G_{Ca}}{2} [1 + \tanh(z)]$$

$$\bar{g}_K = \frac{G_K}{2} [1 - \tanh(z)]$$

Donde:
- $G_{Ca}$ y $G_K$ son parámetros constantes que fijan los máximos posibles para cada neurona.
- **Restricción importante**: La suma es siempre constante:

$$\bar{g}_{Ca} + \bar{g}_K = \frac{G_{Ca} + G_K}{2}$$

El modelo **no crea conductancia total**, sino que la **redistribuye entre K⁺ y Ca²⁺**.

#### Dinámica de z

La variable $z$ evoluciona según una ecuación lenta:

$$\tau_z \frac{dz}{dt} = \tanh(I_{target} - I_{Ca})$$

Con $\tau_z = 5$ s (mucho más lenta que la dinámica eléctrica).

Donde:
- $I_{Ca} = \bar{g}_{Ca} m_{Ca}^3 h_{Ca} (V_s - E_{Ca})$: la corriente de Ca²⁺ actual de esa neurona.
- $I_{target}$: valor objetivo fijo para cada tipo de neurona (parámetro).

#### Lógica del Feedback

- **Si $I_{Ca} < I_{target}$** (poca actividad):
  - $(I_{target} - I_{Ca}) > 0$ → $\tanh(\cdot) > 0$ → $\frac{dz}{dt} > 0$.
  - $z$ aumenta → $\tanh(z)$ aumenta → **$\bar{g}_{Ca}$ sube, $\bar{g}_K$ baja**.
  - La neurona se vuelve **más excitable** (más corriente entrante).

- **Si $I_{Ca} > I_{target}$** (mucha actividad):
  - $(I_{target} - I_{Ca}) < 0$ → $\tanh(\cdot) < 0$ → $\frac{dz}{dt} < 0$.
  - $z$ disminuye → $\tanh(z)$ disminuye → **$\bar{g}_{Ca}$ baja, $\bar{g}_K$ sube**.
  - La neurona se vuelve **menos excitable** (más corriente saliente).

**Interpretación**: Es un mecanismo de **homeostasis**: intenta mantener la actividad cercana a un nivel objetivo.

#### Valores de Parámetros

Los valores utilizados en el modelo son:

- $G_{Ca} = 0.2$ µS y $G_K = 16$ µS.
- $I_{target} = 0.4$ nA para AB/PD.
- $I_{target} = 0.3$ nA para LP.
- $I_{target} = 0.5$ nA para PY.

Estos fueron elegidos de modo que, en equilibrio, el circuito genere un patrón de actividad similar al ritmo pilórico.

### 5.5. Sinapsis Rápidas

Las sinapsis rápidas (de LP y PY, y componente rápida de AB/PD) se modelan como:

$$I_{fast} = \frac{\bar{g}_{fast}}{1 + \exp[s_{fast}(V_{fast} - V_s^{pre})]} (V_s^{post} - E_{syn})$$

Donde:
- $\bar{g}_{fast}$: conductancia máxima de la sinapsis.
- $V_s^{pre}$: voltaje somático **presináptico** (el que envía).
- $V_s^{post}$: voltaje somático **postsináptico** (el que recibe).
- $E_{syn} = -75$ mV: potencial de reversión (inhibidor).
- La función sigmoide representa la apertura/cierre de canales sinápticos en función del voltaje presináptico.

### 5.6. Sinapsis Lentas (solo AB/PD)

Las sinapsis lentas tienen una dinámica adicional:

$$I_{slow} = \bar{g}_{slow} m_{slow} (V_s^{post} - E_{syn})$$

$$\frac{dm_{slow}}{dt} = \frac{k_1(1 - m_{slow})}{1 + \exp[s_{slow}(V_{slow} - V_s^{pre})]} - k_2 m_{slow}$$

Donde:
- $m_{slow}$ se activa (sube) cuando $V_s^{pre}$ está alto.
- Se desactiva (baja) con una constante de tiempo controlada por $k_2$.
- Esto introduce una **componente postsináptica más lenta** que prolonga la inhibición.

Los valores de $\bar{g}_{fast}$ para cada conexión están en la **Tabla 3** del artículo; los de $\bar{g}_{slow}$, $k_1$, $k_2$ en la **Tabla 4**.

---

## 6. Experimentos Biológicos

### 6.1. Procedimiento

Se utilizaron cangrejos adultos machos de *Cancer borealis*:

- Mantenidos a ~13°C en agua de mar artificial.
- Se disecó el sistema nervioso para dejar expuesto el STG y los nervios conectados, incluyendo el stn.
- Se registró la actividad **extracelularmente** (desde nervios motores con electrodos de pinza).
- También se hicieron registros **intracelulares** (microelectrodos en soma de neuronas específicas).

### 6.2. Bloqueo del stn

El bloqueo del nervio estomatogástrico se realizó de dos formas:

1. Transección física con tijeras.
2. Colocación de un pozo de vaselina alrededor del nervio con una solución que bloquea los potenciales de acción (sacarosa hiperosmótica + tetrodotoxina).

### 6.3. Resultados Experimentales (Figura 2)

#### Figura 2A – Trazas Extracelulares

Se registró actividad del nervio lvn (lateral ventricular nerve) que contiene axones de LP, PD y PY:

- **Antes del bloqueo ("control")**: Ritmo pilórico robusto y regular. Se observa una secuencia repetida de grupos de potenciales de acción que corresponden a los tres tipos de neurona, con una frecuencia de ~0.88 Hz en el ejemplo mostrado.

- **Inmediatamente después del bloqueo ("after stn block")**: El ritmo desaparece casi por completo. Las neuronas LP y PD dejan de disparar en ráfagas ordenadas. Las PY pueden mostrar actividad más continua pero desorganizada.

- **Tras ~24 horas de recuperación ("recovery")**: El ritmo reaparece, mostrando nuevamente una secuencia tripásica. Sin embargo, la frecuencia es notablemente más lenta (~0.31 Hz en el ejemplo).

#### Figura 2B – Histograma de Frecuencias

Se presenta un histograma que cuantifica la frecuencia del ritmo en diferentes condiciones:

| Condición | Frecuencia (Hz) | n |
|-----------|-----------------|---|
| Control (pre-bloqueo) | 1.14 ± 0.09 | 28 |
| Inmediatamente post-bloqueo | 0.01 ± 0.01 | 27 |
| 24 horas después | 0.31 ± 0.10 | 13 |
| 48 horas después | 0.37 ± 0.10 | 7 |

Las diferencias entre grupos son **estadísticamente significativas** (ANOVA en rangos: F(3,16) = 10.671, P < 0.001).

**Interpretación**: 
- El bloqueo del stn elimina la neuromodulación → cae la frecuencia casi a cero.
- Con el tiempo, la red se reorganiza → el ritmo reaparece a menor frecuencia.
- La nueva frecuencia (~0.3-0.4 Hz) se parece a la observada in vivo en animales no alimentados (0.26-0.5 Hz).

---

## 7. Resultados del Modelo

### 7.1. Autoensamblaje del Circuito (Figuras 3 y 4)

#### Figura 3A – Neuronas Aisladas en Equilibrio

Cuando cada neurona está completamente aislada (sin sinapsis entre ellas), pero con sus conductancias $\bar{g}_{Ca}$ y $\bar{g}_K$ ya autorreguladas:

- **AB/PD**: Dispara en ráfagas (bursts), con períodos de actividad seguidos de silencio.
- **LP**: También presenta bursting intrínseco, pero con su propio patrón característico.
- **PY**: Dispara de forma **tónica** (potenciales de acción continuos sin pausas).

Esta es la "personalidad intrínseca" de cada neurona cuando está sola.

#### Figura 3B – Inmediatamente Después de Activar las Conexiones Sinápticas

Se encienden las sinapsis según el diagrama de conectividad (Figura 1B):

- Los patrones de actividad de cada neurona cambian **bruscamente**.
- Las neuronas comienzan a inhibirse mutuamente.
- En este momento, $\bar{g}_{Ca}$ y $\bar{g}_K$ aún no tienen tiempo de reajustarse:
  - Las corrientes $I_{Ca}$ devían de sus valores objetivo $I_{target}$.
  - La variable $z$ comienza a moverse lentamente.
  - Pero aún estamos muy cerca de las conductancias iniciales del Panel A.

#### Figura 3C – Red en Estado Estacionario

Después de mucho tiempo (cuando $z$ ha convergido a nuevos valores de equilibrio):

- La red muestra un **ritmo tripásico estable**:
  - AB/PD se activa en bursts.
  - Luego LP entra en su fase de actividad.
  - Luego PY se activa.
  - El ciclo se repite.
- Este patrón es muy similar al ritmo pilórico biológico observado experimentalmente.
- Las conductancias $\bar{g}_{Ca}$ y $\bar{g}_K$ han alcanzado nuevos valores que permiten esta oscilación coordinada.

#### Figura 3D – Neuronas Inmediatamente Después de Desconectar

Se apagan las sinapsis después de que la red haya estado acoplada largo tiempo:

- Cada neurona vuelve a estar aislada.
- Pero sus propiedades intrínsecas **han cambiado** (porque $z$ fue modificado durante el acoplamiento).
- **Comparando A y D**: 
  - LP, por ejemplo, ahora dispara de forma más tónica (menos bursting) que en el Panel A.
  - PY dispara a un ritmo más rápido.
  - AB/PD muestra cambios menores.

**Conclusión**: La experiencia de estar acoplado en red ha "reprogramado" las corrientes intrínsecas de cada neurona mediante el mecanismo de regulación dependiente de la actividad.

### 7.2. Robustez del Autoensamblaje (Figura 4)

#### Figura 4A – Convergencia desde Múltiples Condiciones Iniciales

Se demuestra que el autoensamblaje es **robusto**:

- Se muestran dos **ensayos paralelos** con condiciones iniciales **muy diferentes** en términos de $\bar{g}_{Ca}$ y $\bar{g}_K$ de las tres neuronas.
- En la parte superior de cada columna se ve la actividad inicial: los patrones son distintos entre los dos ensayos.
- En la parte inferior se ve la actividad final: **ambos ensayos convergen al mismo ritmo tripásico**.

Esto demuestra que existe un **atractor** en el espacio de estados (el espacio de valores posibles de $z_{ABPD}, z_{LP}, z_{PY}$).

#### Figura 4B – Evolución Temporal de las Conductancias

Se muestra cómo evolucionan en el tiempo:

- **Panel superior**: Evolución de $\bar{g}_K$ para cada neurona (AB/PD, LP, PY) en ambos ensayos.
- **Panel inferior**: Evolución de $\bar{g}_{Ca}$ para cada neurona.

**Observaciones clave**:
- Las trayectorias temporales son **diferentes** entre los dos ensayos (no son idénticas).
- Algunas conductancias muestran cambios **no monótonos** (suben, bajan, vuelven a subir).
- A pesar de esto, **convergen a los mismos valores finales** para cada neurona.

**Conclusión**: El sistema tiene un comportamiento de atractor. Independientemente del punto de partida, el circuito busca y encuentra una configuración de conductancias que permite el ritmo pilórico estable.

### 7.3. Bloqueo de Neuromodulación y Recuperación (Figura 5)

Esta es la comparación más importante: el modelo reproduce el fenómeno biológico.

#### Figura 5A – Control (stn Intacto / Proctolina Activa)

**Modelo (izquierda)**:
- Se muestran registros simulados (como si fueran extracelulares).
- Voltajes de AB/PD, LP, PY en el tiempo.
- Se observa un ritmo tripásico ordenado.

**Experimento biológico (derecha)**:
- Extracelular desde los nervios lvn y lpg/pyn.
- Registros intracelulares de PD y LP.
- Se observa el mismo patrón tripásico ordenado.

**Concordancia**: El modelo captura correctamente el patrón de control.

#### Figura 5B – Inmediatamente Tras Bloqueo (Proctolina = 0 / stn Bloqueado)

**Modelo (izquierda)**:
- Al quitar la conductancia de proctolina ($\bar{g}_{Proc} = 0$), el ritmo se derrumba.
- AB/PD y LP quedan **silenciosas** (sin actividad).
- **PY queda muy despolarizada** (voltaje muy positivo) y dispara **tónico** a alta frecuencia.

**Experimento (derecha)**:
- Inmediatamente tras bloquear el stn, se observa:
  - LP y PD sin actividad en ráfagas ordenadas.
  - PY con actividad más continua.
  - La organización rítmica ha desaparecido.

**Concordancia**: El modelo reproduce cualitativamente el colapso inmediato del ritmo.

#### Figura 5C – Recuperación Tras Tiempo Prolongado

**Modelo (izquierda)**:
- Con el tiempo (sin la conductancia de proctolina), el mecanismo de regulación dependiente de la actividad actúa:
  - La caída de $I_{Ca}$ en AB/PD y LP (por falta de modulación) causa que $z$ aumente lentamente.
  - Esto eleva $\bar{g}_{Ca}$ y reduce $\bar{g}_K$ en esas neuronas.
  - Las hace de nuevo más excitables.
- Poco a poco, el circuito recupera un patrón rítmico:
  - **Tripásico** (mismo orden de fases que antes).
  - **Pero más lento** (menor frecuencia).
  - **Completamente independiente de la proctolina**.

**Experimento (derecha)**:
- Tras ~24 h de bloqueo del stn, el ritmo pilórico vuelve:
  - Con la misma secuencia de fases (tripásico).
  - A una frecuencia más lenta (~0.31-0.37 Hz vs ~1.14 Hz original).

**Concordancia**: El modelo reproduce de forma cualitativamente correcta cómo la regulación de conductancias puede restaurar el ritmo sin neuromodulación.

---

## 8. Interpretación Biológica y Conceptual

### 8.1. Dos Modos Operacionales

El hallazgo principal es que **la misma red puede operar en dos "modos" fundamentalmente diferentes**:

1. **Modo Dependiente de Neuromodulación** (normal in vivo con moduladores)
   - Ritmo rápido (~1 Hz).
   - Controlado por neuromoduladores externos (proctolina del stn).
   - Flexible, sensible a cambios en los moduladores.

2. **Modo Intrínseco Autónomo** (después de pérdida prolongada de moduladores)
   - Ritmo más lento (~0.3-0.4 Hz).
   - Completamente independiente de entrada moduladora.
   - Generado por propiedades intrínsecas modificadas de las neuronas.

### 8.2. Mecanismo de Transición

La transición entre modos ocurre a través de:

1. **Bloqueo de moduladores**: Caída de actividad neuronal inmediata.
2. **Regulación dependiente de la actividad**: Cambio gradual de $\bar{g}_{Ca}$ y $\bar{g}_K$ durante horas/días.
3. **Reorganización de la red**: Convergencia a una nueva configuración estable.

### 8.3. Plasticidad como Estabilidad

Paradójicamente:
- La **plasticidad** (cambio de conductancias en respuesta a actividad) es normalmente vista como un mecanismo de cambio adaptativo.
- En este modelo, la misma plasticidad actúa como un mecanismo de **estabilidad homeostática**.

La red no "se desmorona" ante la pérdida de moduladores, sino que se reorganiza para mantener función (aunque en un modo diferente).

### 8.4. Contexto Biológico

En la naturaleza:
- Los cangrejos pasan períodos alimentados e inanicionados.
- En estado inanicionado: actividad digestiva baja, el ritmo pilórico es lento (~0.26-0.5 Hz).
- En estado alimentado: actividad digestiva alta, el ritmo es rápido (~0.8-1.0 Hz).

El modelo sugiere que:
- El modo rápido (con moduladores) corresponde al estado alimentado.
- El modo lento (sin moduladores) corresponde al estado inanicionado.
- La red **se autorregula** para adaptar su frecuencia al estado metabólico.

---

## 9. Implicaciones y Limitaciones

### 9.1. Importancia del Modelo

1. **Unificación teórica**: Muestra cómo una sola clase de mecanismo (regulación dependiente de actividad a nivel neuronal) puede explicar:
   - Autorregulación de células individuales.
   - Estabilidad de circuitos completos.
   - Adaptación a cambios en entrada moduladora.

2. **Minimalismo**: Se logra con solo:
   - Tres neuronas modelo.
   - Dos conductancias reguladas por neurona ($\bar{g}_{Ca}$ y $\bar{g}_K$).
   - Sin necesidad de plasticidad sináptica (las sinapsis permanecen fijas).

3. **Predicción experimental**: El modelo predice que debe haber **cambios significativos en las corrientes iónicas** mientras el ganglio se recupera del bloqueo del stn.

### 9.2. Limitaciones

1. **Simplificación**: 
   - El modelo real tiene muchas más neuronas y canales.
   - Solo se regula un subconjunto de conductancias.

2. **Plasticidad sináptica ignorada**:
   - El modelo mantuvo las sinapsis fijas.
   - En la realidad, probablemente también hay cambios en fuerzas sinápticas.

3. **Múltiples segundos mensajeros**:
   - El modelo usa solo Ca²⁺ como señal de actividad.
   - La célula real probablemente usa múltiples cascadas de señalización.

4. **Diferencias entre especies**:
   - Los resultados de Golowasch et al. usan *Cancer borealis*.
   - Otros grupos usan *Panulirus argus* (langosta), con diferentes cinéticas de recuperación.

### 9.3. Direcciones Futuras

- Inclusión de plasticidad sináptica.
- Modelos más realistas con más neuronas y más canales.
- Validación experimental de cambios predichos en corrientes iónicas.
- Estudio de interacciones entre plasticidad intrínseca y sináptica.

---

## 10. Glosario de Términos Clave

| Término | Definición |
|---------|-----------|
| **Bursting** | Patrón de disparo con ráfagas de potenciales de acción intercaladas con períodos de silencio. |
| **Conductancia** | Medida de facilidad con la que los iones pueden fluir a través de un canal (inverso de resistencia). |
| **Equilibrio (potencial de)** | Voltaje al que la concentración de un ion está en equilibrio químico-eléctrico. |
| **Homeostasis** | Tendencia del sistema a mantener condiciones internas estables. |
| **Hodgkin-Huxley** | Modelo matemático clásico de comportamiento de canales iónicos basado en ecuaciones diferenciales. |
| **Inhibición** | Reducción de la probabilidad de que una neurona dispare; hiperpolarización de la membrana. |
| **Neuropéptido** | Moléccula de señalización pequeña (cadena de aminoácidos) liberada por neuronas. |
| **Potencial de acción** | Cambio rápido y grande en el voltaje de membrana que permite la transmisión de información a larga distancia. |
| **Proctolina** | Neuropéptido modulador que afecta la actividad del ganglio estomatogástrico. |
| **Somatogástrico** | Relativo al estómago; en este contexto, refiere al ganglio que lo controla. |
| **Tripásico** | Patrón que tiene tres fases distintas que se repiten cíclicamente. |