from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors


from component import GasFlow, HeatExchangerAdvanced
from design import optimise_hx_anypd, make_heat_exchanger_anypd, make_heat_exchanger_pd

massflow = 0.05
input_T = 250
input_p = 4000


gasflow = GasFlow(massflow, input_T, input_p)

deltah = gasflow.cp * gasflow.mass_flow * (330 - 250)
inlet_area = gasflow.mass_flow / (30 * gasflow.density())

coolant_speed = 1

gap_ratio_range = np.linspace(0.2, 30, 50*2)
area_ratio_range = np.linspace(1, 15, 20*2)
PD_range = np.linspace(0.01, 0.2, 20)

ESM_constrained = np.zeros((gap_ratio_range.size, area_ratio_range.size))
ESM_real = np.zeros((gap_ratio_range.size, area_ratio_range.size))
ESM_best = np.zeros((gap_ratio_range.size, area_ratio_range.size))

for GR_i, GR in enumerate(gap_ratio_range):
    for AR_i, AR in enumerate(area_ratio_range):
        HX_weight_constrained = make_heat_exchanger_anypd((AR, coolant_speed, GR), deltah, gasflow, inlet_area)
        HX_weight_real = make_heat_exchanger_anypd((AR, coolant_speed, GR), deltah, gasflow, inlet_area, output=True)
        ESM_constrained[GR_i, AR_i] = HX_weight_constrained
        ESM_real[GR_i, AR_i] = HX_weight_real

        sweep_ESM_constrained = []
        for PD in PD_range:
            HX_weight_constrained = make_heat_exchanger_pd((AR, coolant_speed, GR), 
                deltah, gasflow, inlet_area, PD, output=True)
            sweep_ESM_constrained.append(HX_weight_constrained + 8000*PD)
        ESM_best[GR_i, AR_i] = min(sweep_ESM_constrained)


print(np.amin(ESM_constrained))
print(np.amin(ESM_best))
fig, ax = plt.subplots()
plt.contourf(gap_ratio_range+1, area_ratio_range, np.transpose(ESM_constrained), locator=ticker.LogLocator())
plt.xlabel("Ratio of tube pitch to tube diameter")
plt.ylabel("Ratio of inlet area to bulk flow area")
plt.title("ESM field seen by optimiser, including constraints")
plt.colorbar()

fig, ax = plt.subplots()
plt.contourf(gap_ratio_range+1, area_ratio_range, np.transpose(ESM_best), locator=ticker.LogLocator())
plt.xlabel("Ratio of tube pitch to tube diameter")
plt.ylabel("Ratio of inlet area to bulk flow area")
plt.title("ESM field seen by optimiser, best-of-20 pressure drops")
plt.colorbar()

plt.show()


