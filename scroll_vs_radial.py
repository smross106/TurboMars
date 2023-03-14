import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from component import ScrollCompressor, GasFlow
from design import optimise_radial

inlet_pressure_range = np.logspace(3, 5.5, 55)
pressure_rise_range = np.linspace(1.1, 10.5, 35)

mass_flow = 0.075


inlet_temperature = 250
print((inlet_pressure_range.shape, pressure_rise_range.shape))
scroll_weight = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))
scroll_power = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))
scroll_outlet_T = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))

radial_weight = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))
radial_power = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))
radial_outlet_T = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))

for ip_i, inlet_pressure in enumerate(inlet_pressure_range):
    for pr_i, pressure_rise in enumerate(pressure_rise_range):
        if pressure_rise < 2.2:
            gasflow_1 = GasFlow(mass_flow, inlet_temperature, inlet_pressure)
            rad1 = optimise_radial(gasflow_1, pressure_rise, 0.8)

            temp_2 = gasflow_1.temperature * ((np.power(pressure_rise, gasflow_1.gm1/gasflow_1.gamma) - 1)/rad1.efficiency + 1)
            
            radial_weight[ip_i, pr_i] = rad1.weight
            radial_power[ip_i, pr_i] = rad1.power
            radial_outlet_T[ip_i, pr_i] = temp_2

        elif pressure_rise < 2.2**2:
            stage_pressure_rise = np.sqrt(pressure_rise)

            gasflow_1 = GasFlow(mass_flow, inlet_temperature, inlet_pressure)
            rad1 = optimise_radial(gasflow_1, stage_pressure_rise, 0.8)

            temp_2 = gasflow_1.temperature * ((np.power(stage_pressure_rise, gasflow_1.gm1/gasflow_1.gamma) - 1)/rad1.efficiency + 1)
            gasflow_2 = GasFlow(mass_flow, temp_2, inlet_pressure*stage_pressure_rise)
            rad2 = optimise_radial(gasflow_2, stage_pressure_rise, 0.8)

            temp_3 = gasflow_2.temperature * ((np.power(stage_pressure_rise, gasflow_2.gm1/gasflow_2.gamma) - 1)/rad2.efficiency + 1)

            radial_weight[ip_i, pr_i] = rad1.weight + rad2.weight
            radial_power[ip_i, pr_i] = rad1.power + rad2.power
            radial_outlet_T[ip_i, pr_i] = temp_3
        elif pressure_rise < 2.2**3:
            stage_pressure_rise = np.power(pressure_rise, 1/3)

            gasflow_1 = GasFlow(mass_flow, inlet_temperature, inlet_pressure)
            rad1 = optimise_radial(gasflow_1, stage_pressure_rise, 0.8)

            temp_2 = gasflow_1.temperature * ((np.power(stage_pressure_rise, gasflow_1.gm1/gasflow_1.gamma) - 1)/rad1.efficiency + 1)
            gasflow_2 = GasFlow(mass_flow, temp_2, inlet_pressure*stage_pressure_rise)
            rad2 = optimise_radial(gasflow_2, stage_pressure_rise, 0.8)

            temp_3 = gasflow_2.temperature * ((np.power(stage_pressure_rise, gasflow_2.gm1/gasflow_2.gamma) - 1)/rad2.efficiency + 1)
            gasflow_3 = GasFlow(mass_flow, temp_3, inlet_pressure*stage_pressure_rise*stage_pressure_rise)
            rad3 = optimise_radial(gasflow_3, stage_pressure_rise, 0.8)

            temp_4 = gasflow_3.temperature * ((np.power(stage_pressure_rise, gasflow_3.gm1/gasflow_3.gamma) - 1)/rad3.efficiency + 1)

            radial_weight[ip_i, pr_i] = rad1.weight + rad2.weight + rad3.weight
            radial_power[ip_i, pr_i] = rad1.power + rad2.power + rad3.power
            radial_outlet_T[ip_i, pr_i] = temp_4
        
        scroll = ScrollCompressor(gasflow_1, pressure_rise)
        scroll.estimate_efficiency()
        scroll.estimate_weight()

        scroll_weight[ip_i, pr_i] = scroll.weight
        scroll_power[ip_i, pr_i] = scroll.work_stage/scroll.efficiency
        temp_2 = gasflow_1.temperature * ((np.power(pressure_rise, gasflow_1.gm1/gasflow_1.gamma) - 1)/scroll.efficiency + 1)
        scroll_outlet_T[ip_i, pr_i] = temp_2
    
        print(100*(ip_i*len(pressure_rise_range)+pr_i)/(len(pressure_rise_range)*len(inlet_pressure_range)))


ESM_radial = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))
ESM_scroll = np.zeros((len(inlet_pressure_range), len(pressure_rise_range)))
for i, ip in enumerate(inlet_pressure_range):
    for j, pr in enumerate(pressure_rise_range):
        ESM_radial[i, j] = radial_weight[i,j] + .149*radial_power[i,j]
        ESM_scroll[i,j] = scroll_weight[i,j] + .149*scroll_power[i,j]

fig, ax = plt.subplots()

plt.contourf(inlet_pressure_range, pressure_rise_range, np.transpose(ESM_radial) / np.transpose(ESM_scroll), levels=np.logspace(-1.5, 1.5, 13), locator=ticker.LogLocator())
plt.colorbar()
plt.contour(inlet_pressure_range, pressure_rise_range, np.transpose(ESM_radial) / np.transpose(ESM_scroll), levels=[1], colors="k")

plt.xlabel("Inlet pressure, kPa")
plt.ylabel("Pressure ratio")
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
ax.xaxis.set_major_formatter(ticks)
plt.xscale("log")
plt.title("Ratio of equivalent system mass of radial stage(s) and scroll compressor \n 75g/s, power at 149kg/kW")

plt.show()