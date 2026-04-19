### RTHybrid Bidirectional Chemical Synapse BO

**Requirements:** libraries KFR DSP (compiled with -DCMAKE_POSITION_INDEPENDENT_CODE=ON) and Limbo (set the environment variable LIMBO_DIR to its the downloaded repo)
**Limitations:** The user should change the -lkfr_dsp_... flag from ard2 to the appropriate one for their system in the Makefile.

![RTHybrid Bidirectional Chemical Synapse BO GUI](bidirectional_chemical_synapse_BO.png)

<!--start-->
<p><b>RTHybrid Bidirectional Chemical Synapse BO</b><br>RTHybrid module for RTXI that implements a bidirectional chemical synapse model and runs Bayesian Optimization (BO) online to fit synaptic parameters based on captured voltage/current signals.</p>
<!--end-->

#### Input
1. input(0) - Voltage 1 (V) : Membrane potential 1
2. input(1) - Voltage 2 (V) : Membrane potential 2
3. input(2) - Voltage min 1 (V) : Dynamic min 1
4. input(3) - Voltage max 1 (V) : Dynamic max 1
5. input(4) - Voltage min 2 (V) : Dynamic min 2
6. input(5) - Voltage max 2 (V) : Dynamic max 2

#### Output
1. output(0) - Current 1->2 (nA) : Total synaptic current 1->2
2. output(1) - Current 2->1 (nA) : Total synaptic current 2->1

#### Parameters
1. BO initial samples - Number of initialization samples for BO
2. BO iterations - Number of BO iterations after initial sampling
3. BO evaluation time (ms) - Time to record signals per evaluation
4. BO stabilization time (ms) - Wait time after setting params before recording
5. BO search phase (1/0) - 1 = Enable, 0 = Disable
6. BO current min to achieve 1->2 (nA) - Target minimum current for direction 1->2
7. BO current max to achieve 1->2 (nA) - Target maximum current for direction 1->2
8. BO current min to achieve 2->1 (nA) - Target minimum current for direction 2->1
9. BO current max to achieve 2->1 (nA) - Target maximum current for direction 2->1
10. BO cutoff frequency 1 (kHz) - To separate the I_fast and I_slow for BO in synapse 1->2
11. BO cutoff frequency 2 (kHz) - To separate the I_fast and I_slow for BO in synapse 2->1
12. Dynamic voltage min and max 1 (1/0) - 1 = Enable, 0 = Disable; necessary for BO
13. Voltage min 1 (V) - Necessary for BO
14. Voltage max 1 (V) - Necessary for BO
15. Dynamic voltage min and max 2 (1/0) - 1 = Enable, 0 = Disable; necessary for BO
16. Voltage min 2 (V) - Necessary for BO
17. Voltage max 2 (V) - Necessary for BO
18. Current min 1->2 (nA) - Fixed output clamp min for current 1->2
19. Current max 1->2 (nA) - Fixed output clamp max for current 1->2
20. Current min 2->1 (nA) - Fixed output clamp min for current 2->1
21. Current max 2->1 (nA) - Fixed output clamp max for current 2->1
22. factor in dt (ms) = period (ms) * factor - Factor for calculating dt form the period; dt in ms
23. Use I_fast 1->2 (1/0) - 1 = Enable, 0 = Disable
24. Use I_slow 1->2 (1/0) - 1 = Enable, 0 = Disable
25. Use I_fast 2->1 (1/0) - 1 = Enable, 0 = Disable
26. Use I_slow 2->1 (1/0) - 1 = Enable, 0 = Disable
27. E_syn 1->2 (V) - Synaptic reversal potential (1->2)
28. g_fast 1->2 (nS) - Fast conductance (1->2)
29. s_fast 1->2 (1/V) - Fast sigmoid slope (1->2)
30. V_fast 1->2 (V) - Fast sigmoid threshold (1->2)
31. g_slow 1->2 (nS) - Slow conductance (1->2)
32. k1 1->2 (1/ms) - Slow gating opening rate (1->2)
33. k2 1->2 (1/ms) - Slow gating closing rate (1->2)
34. s_slow 1->2 (1/V) - Slow sigmoid slope (1->2)
35. V_slow 1->2 (V) - Slow sigmoid threshold (1->2)
36. E_syn 2->1 (V) - Synaptic reversal potential (2->1)
37. g_fast 2->1 (nS) - Fast conductance (2->1)
38. s_fast 2->1 (1/V) - Fast sigmoid slope (2->1)
39. V_fast 2->1 (V) - Fast sigmoid threshold (2->1)
40. g_slow 2->1 (nS) - Slow conductance (2->1)
41. k1 2->1 (1/ms) - Slow gating opening rate (2->1)
42. k2 2->1 (1/ms) - Slow gating closing rate (2->1)
43. s_slow 2->1 (1/V) - Slow sigmoid slope (2->1)
44. V_slow 2->1 (V) - Slow sigmoid threshold (2->1)

#### States
1. BO evaluations completed - Finishes when this is initial samples + iterations

