**1. ¿Por qué usamos `add_synaptic_input` si la neurona receptora es un CSV?**

* **Respuesta:** Es **código muerto (inútil)**. Al ser una grabación, el valor del CSV sobrescribe el estado de la neurona en cada paso, ignorando cualquier corriente que le intentes inyectar.

**2. ¿Qué distancia buscamos para que funcione bien en un circuito biohíbrido?**

* **Respuesta:** Buscas la **Mínima Distancia posible (cercana a 0)**.
* En un circuito biohíbrido (lazo cerrado), necesitas **sincronización perfecta**: cuando la neurona artificial dispara, la corriente debe inyectarse en la biológica *al instante*.
* Una distancia alta indicaría "lag" (retraso) o distorsión, lo que desincronizaría la comunicación entre la parte viva y la artificial.


**3. ¿Por qué comparar  vs  y no contra la real ()?**

* **Respuesta:** Para **calibrar la herramienta (la sinapsis)**.
* Estás comprobando que tu modelo matemático convierte voltaje en corriente correctamente.
* Comparar con la neurona real () es absurdo porque es una grabación del pasado que no reacciona a tu simulación; no te diría si tu sinapsis funciona.