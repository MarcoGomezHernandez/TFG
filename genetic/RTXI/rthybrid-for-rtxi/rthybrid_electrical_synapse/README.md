### RTHybrid electrical synapse model

**Requirements:** None  
**Limitations:** None  

![RTHybrid electrical synapse model GUI](rthybrid_electrical_synapse.png)

<!--start-->
<p><b>RTHybrid electrical synapse model</b><br>RTHybrid module for RTXI that implements an electrical synapse model.</p>
<!--end-->

#### Input
1. input(0) - Post-synaptic Voltage (V)
2. input(1) - Pre-synaptic Voltage (V)
3. input(2) - Scale (Pre to Post) : Scale from the pre-synaptic neuron to the post-synaptic one, as given by the RTHybridAmplitudeScale module.
4. input(3) - Offset (Pre to Post) : Offset from the pre-synaptic neuron to the post-synaptic one, as given by the RTHybridAmplitudeScale module.

#### Output
1. output(0) - Current (nA) : Synaptic current (in nA)

#### Parameters
1. g (uS) - Conductance (uS)

#### States
1. Current (nA) 
2. offset
3. scale
