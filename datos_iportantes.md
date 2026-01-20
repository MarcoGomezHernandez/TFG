1. ¿Por qué usamos add_synaptic_input si la neurona receptora es un CSV?**

Es código muerto (inútil). Al ser una grabación, el valor del CSV sobrescribe el estado de la neurona en cada paso, ignorando cualquier corriente que le intentes inyectar.


2. ¿Qué distancia buscamos para que funcione bien en un circuito biohíbrido?

Buscas la **Mínima Distancia posible (cercana a 0).

Se calcula la distancia entre la corriente: en parametrize_chemicalSynapsis_HR.cpp y parametrize_chemicalSynapsis_HR-PD.cpp: Se utiliza el potencial de la neurona PRESINÁPTICA; en el otro, la POSTSINÁPTICA (En los modelos donde se usa la presináptica (HR, HR-PD), se busca caracterizar la dinámica de emisión del modelo químico asegurando que la corriente se genere sin retardos espurios respecto al disparo de origen ($V_{pre} \approx I_{fast}$); mientras que en el modelo que usa la postsináptica (PD-HR), se busca garantizar la eficacia de excitación, seleccionando parámetros que aseguren que la corriente inyectada sea suficiente para detonar causalmente un potencial de acción en la neurona destino ($I_{fast} \approx V_{post}$).).

En ambos, se compara esta con ifast de la sinapsis (Se utiliza ifast porque su dinámica temporal rápida es topológicamente comparable a los picos de potencial de acción ($V$) porque depende directamente de ellos, a diferencia de la corriente total que incluye componentes lentas de meseta ($I_{slow}$), que no dependen directamente de los picos, lo que permite aislar la respuesta sináptica inmediata para minimizar eficazmente el error geométrico durante el ajuste de los parámetros cinéticos sin la interferencia de formas de onda incompatibles)


3. ¿Por qué comparar $V_{pre}$ vs $I_{syn}$ y no contra la real ($V_{post}$)?

Para calibrar la herramienta (la sinapsis). Estás comprobando que tu modelo matemático convierte voltaje en corriente correctamente.Comparar con la neurona real ($V_{post}$) es absurdo porque es una grabación del pasado que no reacciona a tu simulación; no te diría si tu sinapsis funciona.