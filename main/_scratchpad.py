import matplotlib.pyplot as plt
import numpy as np



### Plotting state of the art 

fig, ax = plt.subplots(figsize=(2.96, 2.22))

fig.tight_layout()
#fig.subplots_adjust(left=0.15, bottom=0.15, right=0.96, top=0.9)
fig.subplots_adjust(left=0.2, bottom=0.21, right=0.96, top=0.9)
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular'
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

flow_rates_cryo = [0.0624, 0.458333333, 0.06, 0.0125, 0.5, 0.088, 0.02, 0.619]
specific_works_cryo = [28846.15385, 3173.236364, 24000, 3135.6, 4435.2, 19636.36364, 7200, 12969.30533]

flow_rates_other = [0.041666667, 0.458333333, 0.08, 0.5, 0.06732, 0.083]
specific_works_other = [4320, 5498.181818, 5175, 12888, 12834.2246, 5204.819277]

flow_rates_thesis = [1000/24, 1000/12, 1000/8, 1000/6, 1000/3]
specific_works_thesis = [517.67, 493.22, 490.47, 478.02, 473.18]

plt.scatter(flow_rates_cryo, specific_works_cryo, label="Cryofreezer")
plt.scatter(flow_rates_other, specific_works_other, label="Other")
plt.scatter(flow_rates_thesis, specific_works_thesis, label="This work", c="k",s=25)
plt.scatter(flow_rates_thesis[1], specific_works_thesis[1], c="k", marker="*", s=125)

plt.legend(loc="lower left", fontsize=10, borderpad=0.1, columnspacing=0, handletextpad=0.4, handlelength=1)
#plt.grid(visible=True, which="minor", linewidth=1)
plt.grid(visible=True, which="major", c="k")

#plt.xlim(0, 0.75)
#plt.xticks(np.arange(0, 0.8, 0.1), ["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7"], usetex=True, fontsize=10)
#plt.ylim(0, 30000)
#plt.yticks(np.arange(0, 35000, 5000), np.arange(0, 35000, 5000), usetex=True, fontsize=10)

ax.set_xlim(0.01, 1000)
#plt.xticks(np.logspace(-2, 3, 6), [0.01, 0.1, 1, 10, 100, 1000], usetex=True, fontsize=10)

plt.ylim(35, 35000)
#plt.yticks(np.logspace(2, 4, 3), [r"$10^2$", r"$10^3$", r"$10^4$"], usetex=True, fontsize=10)
plt.xscale("log")
plt.yscale("log")

#plt.vlines(1000/24, 41, 5000, colors="tab:green")
#plt.hlines([41, 5000], 1000/24, 1000, colors="tab:green")
plt.fill_between([1000/24, 3000], [41, 41], [5000, 5000], alpha=0.5)

plt.hlines(2990, 0.005, 3000, colors="tab:blue", linestyles="--")
plt.hlines(41, 0.005, 3000, colors="tab:red", linestyles="--")
plt.text(1.5, 2990, "Cryofreezer floor", c="tab:blue", fontsize=10, usetex=True)
plt.text(0.8, 41, "Isothermal process floor", c="tab:red", fontsize=10, usetex=True)
plt.text(1000/20, 600, "Starship", fontsize=10, usetex=True)

plt.xlabel(r"Flow rate, kg/hr", usetex=True, fontsize=10)
plt.ylabel(r"Specific work, kJ/kg $\textrm{CO}_2$", usetex=True, fontsize=10)
plt.title(r"Existing Atmospheric ISRU Systems", usetex=True, fontsize=12)
plt.show()




"""
### Heat exchanger optimisation testing

import CoolProp.CoolProp as CP
import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors

from component import HeatExchangerAdvanced, GasFlow
from design import *


flowin = GasFlow(0.05, 250, 5000)

pd_range = np.linspace(0.01, 0.2, 20)#[0.02, 0.05, 0.1, 0.15, 0.2]
weight_range = []
length_range = []
power_range = []
ESM_range = []

for pd in pd_range:
    hx1 = optimise_hx_pd(flowin, 870*50, 1e-2, pd)
    #print(hx1)
    #print(hx1.coolant_velocity)
    weight_range.append(hx1.weight)
    length_range.append(hx1.duct_length + hx1.diffuser_length + hx1.nozzle_length)
    power_range.append(hx1.pumping_power)
    ESM_range.append(hx1.ESM)

print(min(ESM_range), min(weight_range), min(power_range), min(length_range))

figs, axs = plt.subplots(2, 2)
axs[0,0].plot(pd_range, weight_range)
axs[0,0].set_ylabel("Weight (kg)")
axs[0,0].set_ylim(0, 3000)

axs[0,1].plot(pd_range, length_range)
axs[0,1].set_ylabel("Length (m)")
axs[0,1].set_ylim(0, 50)

axs[1,0].plot(pd_range, weight_range)
axs[1,0].set_ylabel("Power (W)")
axs[1,0].set_ylim(0, 1500)

axs[1,1].plot(pd_range, ESM_range)
axs[1,1].set_ylabel("ESM (kg-e)")
axs[1,1].set_ylim(0, 4000)
plt.show()

"""

"""
### Heat exchanger testing

p_co2 = [2695, 4727, 7934, 12805, 19948, 30106, 44156, 63111, 88119, 81317098, 0, 0, 25069347, 71847802, 128805935, 168665455]
s_co2 = [3.283, 3.165, 3.056, 2.955, 2.861, 2.774, 2.693, 2.616, 2.545, 2.477, 2.413, 2.352, 2.292, 2.234, 2.177, 2.138]
T_co2 = [143.15, 148.15, 153.15, 158.15, 163.15, 168.15, 173.15, 178.15, 183.15, 188.15, 193.15, 198.15, 203.15, 208.15, 213.15, 216.55]


T_start = 190
T_end = 260

s_start = PropsSI("Smass", "T|gas", T_start, "P", 900, "CO2")
s_end = PropsSI("Smass", "T|gas", T_end, "P", 550e3, "CO2")

n_thermo_points = 200
n_velocities = 200
dt = 0.1
D = 1

s_range = np.linspace(s_start, s_end, n_thermo_points)
T_range = np.linspace(T_start, T_end, n_thermo_points)
p_range = np.array([PropsSI("P", "T|gas", T_range[i], "Smass", s_range[i], "CO2") for i in range(n_thermo_points)])
d_range = np.array([PropsSI("D", "T|gas", T_range[i], "Smass", s_range[i], "CO2") for i in range(n_thermo_points)])

v_range = np.linspace(0.1, 50, n_velocities)

hx_HTC = np.zeros((T_range.size, v_range.size))
L_power_q_D = np.zeros((T_range.size, v_range.size))
L_D_0_1 = np.zeros((T_range.size, v_range.size))

for T_index, T in enumerate(T_range):
    gasflow = GasFlow(0.05, T, p_range[T_index])
    volumetric_flow = 0.05 / d_range[T_index]
    for v_index, V in enumerate(v_range):
        flow_area = volumetric_flow / V
        HX = HeatExchangerAdvanced(1, gasflow, flow_area, flow_area, neglect_peak_velocity=True)

        HX.duct_length = 1
        HX.n_rows = np.ceil(HX.duct_length / HX.tube_diameter)
        HX.n_tubes = HX.n_rows * HX.n_tubes_wide

        HTC = HX.heat_transfer_per_unit_area

        hx_HTC[T_index, v_index] = HTC
        
        L_power = np.power(V, gasflow.gm1/gasflow.gamma) * T * gasflow.cp * gasflow.density() / (
            4 * HTC * np.power(dt, -gasflow.gamma)
        )
        
        L_real = np.power(L_power*D, gasflow.gamma/gasflow.gm1)
        
        L_power_q_D[T_index, v_index] = L_power
        L_D_0_1[T_index, v_index] = L_real
"""

"""fig, ax = plt.subplots()
plt.plot(T_range, p_range/1000)
plt.title("Optimal compression path")
plt.yscale("log")
plt.ylabel("Pressure (kPa)")
plt.xlabel("Gas temperature (K)")

fig, ax = plt.subplots()
plt.contourf(T_range, v_range, np.transpose(hx_HTC))
plt.title("Heat transfer W/m2")
plt.xlabel("Gas temperature (K)")
plt.ylabel("Gas velocity in cylinder (m/s)")
plt.colorbar()

fig, ax = plt.subplots()
plt.contourf(T_range, v_range, np.transpose(L_power_q_D), locator=ticker.LogLocator())
plt.title("$L^{\\frac{\gamma - 1}{\gamma}}/D$")
plt.xlabel("Gas temperature (K)")
plt.ylabel("Gas velocity in cylinder (m/s)")
plt.colorbar()"""

"""


fig, ax = plt.subplots()
print(np.amax(L_D_0_1))
print(np.amin(L_D_0_1))
cs = plt.contourf(T_range, v_range, np.transpose(L_D_0_1),  locator=ticker.LogLocator(numticks=27))
plt.colorbar()
plt.xlabel("Gas temperature (K)")
plt.ylabel("Gas velocity in cylinder (m/s)")
plt.title("Piston L, D=10mm")"""

"""p = cs.collections[3].get_paths()[0]

v = p.vertices
x = v[:,0]
y = v[:,1]

print(x)
print(y)

fig, ax = plt.subplots()
plt.plot(x,y)"""

"""fig, ax = plt.subplots()
plt.plot(T_range, d_range)

plt.show()"""







### Optimisation simplex plotting for axial


"""



from component import *
from design import *
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl


gas = GasFlow(0.05, 250, 1200)

flow_coeff_range =  np.linspace(0.05, 2.25, 500)
work_coeff_range = np.linspace(0, 0.5, 550)

weights = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))
effs = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))
ESM_norm = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))
ESM_cons = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))

optimal = optimise_axial(gas, 1.1, 0.85, [1, .149, .121, 5.0])
optimal_speed = optimal.speed
#print(optimal)
#print(optimal.flow_coeff, optimal.work_coeff)



# Labelled points
flow_x = [optimal.flow_coeff, optimal.flow_coeff, 0.7, optimal.flow_coeff, 1.5]
work_y = [optimal.work_coeff, 0.1, optimal.work_coeff, 0.02, 0.25]
label_n = ["I", "II", "III", "IV", "V"]


bad = AxialStage(gas, 1.1, optimal_speed, work_y[4], flow_x[4], 0.85)
bad.weight_estimate()
bad.estimate_efficiency()
bad.estimate_ESM()

print(bad.R_mean_in * bad.speed_rad * bad.flow_coeff)
print(bad.flow_coeff, bad.work_coeff)



for flow_i, flow in enumerate(flow_coeff_range):
    for work_i, work in enumerate(work_coeff_range):
        esm_cons = make_axial((optimal_speed, work, flow), gas, 1.1, 0.85, [1, .149, .121, 5.0])
        esm_nocons = make_axial((optimal_speed, work, flow), gas, 1.1, 0.85, [1, .149, .121, 5.0], output=True)

        ESM_norm[flow_i, work_i] = esm_nocons
        ESM_cons[flow_i, work_i] = esm_cons

        ax = AxialStage(gas, 1.1, optimal_speed, work, flow, 0.85)
        ax.weight_estimate()
        ax.estimate_efficiency()
        ax.estimate_ESM()

        weights[flow_i, work_i] = ax.weight
        effs[flow_i, work_i] = ax.efficiency
    if flow_i%10:
        print(flow_i/5,"%")


specific_speed = 2.0
specific_diameter = 2.0

deltah = gas.delta_h_PR(1.1) / 0.85

speed_guess = specific_speed / (np.sqrt(gas.mass_flow / gas.density_deltah(deltah)) * np.power(deltah, -0.75)) 
work_coeff_guess = 1 / np.power(specific_speed * specific_diameter, 2)
flow_coeff_guess = 1 / (specific_speed * np.power(specific_diameter, 3))

print(flow_coeff_guess, work_coeff_guess)


fig, axs = plt.subplots(2, 2, figsize=(6*1.2,4*1.2), sharex=True, sharey=True)
fig.tight_layout(w_pad=1.0, h_pad=1.5)
fig.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95)
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular'
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

mapw = axs[0,0].contourf(flow_coeff_range, work_coeff_range, np.transpose(weights), locator=ticker.LogLocator(), levels=np.logspace(-2, 2, 9))
divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mapw, cax=cax)

axs[0,0].scatter(flow_x, work_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,0].annotate(txt, (flow_x[i]+0.02, work_y[i]))


mape = axs[0,1].contourf(flow_coeff_range, work_coeff_range, np.transpose(effs))
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mape, cax=cax)

axs[0,1].scatter(flow_x, work_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,1].annotate(txt, (flow_x[i]+0.02, work_y[i]))


map1 = axs[1,0].contourf(flow_coeff_range, work_coeff_range, np.transpose(ESM_norm), locator=ticker.LogLocator(), levels=np.logspace(0, 6, 13))
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map1, cax=cax)

axs[1,0].scatter(flow_x[1:], work_y[1:], marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,0].annotate(txt, (flow_x[i]+0.02, work_y[i]))

axs[1,0].plot([flow_coeff_guess, flow_coeff_guess*0.8, 2.0, flow_coeff_guess*0.8, flow_coeff_guess],
            [work_coeff_guess, work_coeff_guess*0.8, work_coeff_guess*0.8, 0.6, work_coeff_guess], 
             '--o', c="k", label="Initial simplex")

axs[1,0].scatter([optimal.flow_coeff], [optimal.work_coeff], color="k", marker="*", label="Optimised")

axs[1,0].contour(flow_coeff_range, work_coeff_range, np.transpose(ESM_norm), levels=[3e5], c="k")
#axs[1,0].text(0.2, 0.2, "Low efficiency", rotation=80)
#axs[1,0].text(1.3, 0.25, "Low efficiency", rotation=-45)

axs[1,0].legend()

map2 = axs[1,1].contourf(flow_coeff_range, work_coeff_range, np.transpose(ESM_cons), locator=ticker.LogLocator(), levels=np.logspace(0, 6, 13))
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map2, cax=cax)

axs[1,1].contour(flow_coeff_range, work_coeff_range, np.transpose(ESM_cons), levels=[1.1e6], c="k")

axs[1,1].scatter(flow_x, work_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,1].annotate(txt, (flow_x[i]+0.02, work_y[i]))


axs[1,1].text(1.45, 0.25, "Low efficiency", rotation=-30)
axs[1,1].text(0.3, 0.15, "Negative hub radius", rotation=40)
axs[1,1].text(1.2, 0.005, "Short blades", rotation=3)

axs[0,0].set_xlim(0, 2.0)
axs[0,1].set_xlim(0, 2.0)
axs[1,0].set_xlim(0, 2.0)
axs[1,1].set_xlim(0, 2.0)

#axs[0,0].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)
#axs[0,1].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)
axs[1,0].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)
axs[1,1].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)

axs[0,0].set_ylim(0, 0.4)
axs[0,1].set_ylim(0, 0.4)
axs[1,0].set_ylim(0, 0.4)
axs[1,1].set_ylim(0, 0.4)

axs[0,0].set_yticklabels(["0.0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"], usetex=True)
#axs[0,1].set_yticklabels(["0.0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"], usetex=True)
axs[1,0].set_yticklabels(["0.0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"], usetex=True)
#axs[1,1].set_yticklabels(["0.0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"], usetex=True)

#axs[0,0].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
axs[0,0].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[0,0].set_title(r"Weight at fixed speed", usetex=True)

#axs[0,1].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
#axs[0,1].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[0,1].set_title(r"Efficiency at fixed speed", usetex=True)


axs[1,0].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
axs[1,0].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[1,0].set_title(r"ESM at fixed speed", usetex=True)

axs[1,1].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
#axs[1,1].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[1,1].set_title(r"ESM at fixed speed with limits", usetex=True)

#plt.savefig("test.svg")
print("A")
plt.show()


exit()

"""


### Optimisation simplex plotting for centrifugal, flow vs work coefficient


"""



from component import *
from design import *
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl


gas = GasFlow(0.05, 250, 30000)

flow_coeff_range =  np.linspace(0.01, 2.0, 100)
work_coeff_range = np.linspace(0, 0.85, 120)

weights = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))
effs = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))
ESM_norm = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))
ESM_cons = np.zeros(shape=(flow_coeff_range.size, work_coeff_range.size))

optimal = optimise_radial(gas, 1.502, 0.85, [1, .149, .121, 5.0])
optimal_speed = optimal.speed
print(optimal)
print(optimal.flow_coeff, optimal.work_coeff, optimal.R_hub_inlet, optimal.diffusion_ratio)
print(make_radial((optimal_speed, optimal.work_coeff, optimal.flow_coeff, optimal.R_hub_inlet, optimal.diffusion_ratio), gas, 1.502, 0.85, [1, .149, .121, 5.0]))



# Labelled points

flow_x = [optimal.flow_coeff, optimal.flow_coeff, 1.50]
work_y = [optimal.work_coeff, 0.25, optimal.work_coeff] 
label_n = ["I", "II", "III"]


bad = RadialStage(gas, 1.502, optimal_speed, work_y[2], flow_x[2], 0.7, optimal.R_hub_inlet, optimal.diffusion_ratio)
bad.weight_estimate(381e6, 2810)
bad.estimate_efficiency()
bad.estimate_ESM()
#print(bad)

for flow_i, flow in enumerate(flow_coeff_range):
    for work_i, work in enumerate(work_coeff_range):
        esm_cons = make_radial((optimal_speed, work, flow, optimal.R_hub_inlet, optimal.diffusion_ratio), gas, 1.502, 0.7, [1, .149, .121, 5.0])
        esm_nocons = make_radial((optimal_speed, work, flow, optimal.R_hub_inlet, optimal.diffusion_ratio), gas, 1.502, 0.7, [1, .149, .121, 5.0], output=True)

        ESM_norm[flow_i, work_i] = esm_nocons
        ESM_cons[flow_i, work_i] = esm_cons

        rad = RadialStage(gas, 1.502, optimal_speed, work, flow, 0.85, optimal.R_hub_inlet, optimal.diffusion_ratio)
        rad.weight_estimate(381e6, 2810)
        rad.estimate_efficiency()
        rad.estimate_ESM()

        weights[flow_i, work_i] = rad.weight
        effs[flow_i, work_i] = rad.efficiency
    if flow_i%10==0:
        print(100*flow_i/flow_coeff_range.size,"%")


specific_speed = 0.5
specific_diameter = 5.0

deltah = gas.delta_h_PR(1.1) / 0.85

speed_guess = specific_speed / (np.sqrt(gas.mass_flow / gas.density_deltah(deltah)) * np.power(deltah, -0.75)) 
work_coeff_guess = 1 / np.power(specific_speed * specific_diameter, 2)
flow_coeff_guess = 1 / (specific_speed * np.power(specific_diameter, 3))

print(flow_coeff_guess, work_coeff_guess)


fig, axs = plt.subplots(2, 2, figsize=(6*1.2,4*1.2), sharex=True, sharey=True)
fig.tight_layout(w_pad=1.0, h_pad=1.5)
fig.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95)
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular'
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

mapw = axs[0,0].contourf(flow_coeff_range, work_coeff_range, np.transpose(weights), locator=ticker.LogLocator(), levels=[10**-1, 10**-0.5, 1, 10**0.5, 10, 10**1.5, 10**2, 10**3, 10**4])
divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mapw, cax=cax)


axs[0,0].scatter(flow_x, work_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,0].annotate(txt, (flow_x[i]+0.02, work_y[i]-0.05))



mape = axs[0,1].contourf(flow_coeff_range, work_coeff_range, np.transpose(effs))
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mape, cax=cax)


axs[0,1].scatter(flow_x, work_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,1].annotate(txt, (flow_x[i]+0.02, work_y[i]-0.05))
    


map1 = axs[1,0].contourf(flow_coeff_range, work_coeff_range, np.transpose(ESM_norm), locator=ticker.LogLocator(), levels=np.logspace(0, 6, 13))
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map1, cax=cax)


axs[1,0].scatter(flow_x[1:], work_y[1:], marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,0].annotate(txt, (flow_x[i]+0.02, work_y[i]-0.05))


axs[1,0].plot([flow_coeff_guess, flow_coeff_guess*0.8, 2.0, flow_coeff_guess*0.8, flow_coeff_guess],
            [work_coeff_guess, work_coeff_guess*0.8, work_coeff_guess*0.8, 0.85, work_coeff_guess], 
             '--o', c="k", label="Starting simplex")

axs[1,0].scatter([optimal.flow_coeff], [optimal.work_coeff], color="k", marker="*", label="Optimised result")
#axs[1,0].plot([flow_coeff_guess], [work_coeff_guess], "-o", color="tab:red", label="Specific speed/diameter guess")

axs[1,0].contour(flow_coeff_range, work_coeff_range, np.transpose(ESM_norm), levels=[3e5], c="k")
#axs[1,0].text(0.2, 0.2, "Low efficiency", rotation=80)
#axs[1,0].text(1.3, 0.25, "Low efficiency", rotation=-45)

axs[1,0].legend(loc="lower right")

map2 = axs[1,1].contourf(flow_coeff_range, work_coeff_range, np.transpose(ESM_cons), locator=ticker.LogLocator(), levels=np.logspace(0, 6, 13))
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map2, cax=cax)

#axs[1,1].contour(flow_coeff_range, work_coeff_range, np.transpose(ESM_cons), levels=[1.1e6], c="k")


axs[1,1].scatter(flow_x, work_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,1].annotate(txt, (flow_x[i]+0.02, work_y[i]-0.05))
    


#axs[1,1].text(0.75, 0.8, "Negative radius increase", rotation=0)
axs[1,1].text(0.75, 0.05, "Short blades", rotation=0)
#axs[1,1].text(1.2, 0.005, "Short blades", rotation=3)

axs[0,0].set_xlim(0, 2.0)
axs[0,1].set_xlim(0, 2.0)
axs[1,0].set_xlim(0, 2.0)
axs[1,1].set_xlim(0, 2.0)

#axs[0,0].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)
#axs[0,1].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)
axs[1,0].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)
axs[1,1].set_xticklabels(np.arange(0, 2.25, 0.25), usetex=True)

axs[0,0].set_ylim(0, 0.88)
axs[0,1].set_ylim(0, 0.88)
axs[1,0].set_ylim(0, 0.88)
axs[1,1].set_ylim(0, 0.88)

axs[0,0].set_yticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"], usetex=True)
#axs[0,1].set_yticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"], usetex=True)
#axs[1,0].set_yticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"], usetex=True)
axs[1,1].set_yticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"], usetex=True)

#axs[0,0].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
axs[0,0].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[0,0].set_title(r"Weight at fixed speed/diffusion ratio", usetex=True)

#axs[0,1].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
#axs[0,1].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[0,1].set_title(r"Efficiency at fixed speed/diffusion ratio", usetex=True)


axs[1,0].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
axs[1,0].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[1,0].set_title(r"ESM at fixed speed/diffusion ratio", usetex=True)

axs[1,1].set_xlabel(r"Flow coefficient $\Phi$", usetex=True)
#axs[1,1].set_ylabel(r"Work coefficient $\Psi$", usetex=True)
axs[1,1].set_title(r"ESM at fixed speed/diffusion ratio with limits", usetex=True)

#plt.savefig("test2.svg")
print("A")
plt.show()


"""

### Optimisation simplex plotting for centrifugal, speed vs diffusion coefficient






"""

from component import *
from design import *
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl


gas = GasFlow(0.05, 250, 30000)

diff_coeff_range =  np.linspace(1, 15, 500)
speed_range = np.logspace(1, 5, 520)

weights = np.zeros(shape=(diff_coeff_range.size, speed_range.size))
effs = np.zeros(shape=(diff_coeff_range.size, speed_range.size))
ESM_norm = np.zeros(shape=(diff_coeff_range.size, speed_range.size))
ESM_cons = np.zeros(shape=(diff_coeff_range.size, speed_range.size))

optimal = optimise_radial(gas, 1.502, 0.85, [1, .149, .121, 5.0])
optimal_speed = optimal.speed
print(optimal)
print(optimal.flow_coeff, optimal.work_coeff, optimal.R_hub_inlet, optimal.diffusion_ratio)
print(make_radial((optimal_speed, optimal.work_coeff, optimal.flow_coeff, optimal.R_hub_inlet, optimal.diffusion_ratio), gas, 1.502, 0.85, [1, .149, .121, 5.0]))



# Labelled points
diff_x = [optimal.diffusion_ratio, 5, 12, optimal.diffusion_ratio]
speed_y = [optimal.speed, optimal.speed, optimal.speed, 25000]
label_n = ["I", "IV", "V", "VI"]


bad = RadialStage(gas, 1.502, speed_y[3], optimal.work_coeff, optimal.flow_coeff, 0.7, optimal.R_hub_inlet, diff_x[3])
#bad = RadialStage(gas, 1.502, 100000, optimal.work_coeff, optimal.flow_coeff, 0.7, optimal.R_hub_inlet, diff_x[2])
bad.weight_estimate(381e6, 2810)
bad.estimate_efficiency()
bad.estimate_ESM()
print(bad)


for diff_i, diff in enumerate(diff_coeff_range):
    for speed_i, speed in enumerate(speed_range):
        esm_cons = make_radial((speed, optimal.work_coeff, optimal.flow_coeff, optimal.R_hub_inlet, diff), gas, 1.502, 0.7, [1, .149, .121, 5.0])
        esm_nocons = make_radial((speed, optimal.work_coeff, optimal.flow_coeff, optimal.R_hub_inlet, diff), gas, 1.502, 0.7, [1, .149, .121, 5.0], output=True)

        ESM_norm[diff_i, speed_i] = esm_nocons
        ESM_cons[diff_i, speed_i] = esm_cons

        rad = RadialStage(gas, 1.502, speed, optimal.work_coeff, optimal.flow_coeff, 0.7, optimal.R_hub_inlet, diff)
        rad.weight_estimate(381e6, 2810)
        rad.estimate_efficiency()
        rad.estimate_ESM()

        weights[diff_i, speed_i] = rad.weight
        effs[diff_i, speed_i] = rad.efficiency
    
    if diff_i%10==0:
        print(100*diff_i/diff_coeff_range.size, "%")


specific_speed = 0.5
specific_diameter = 5.0

deltah = gas.delta_h_PR(1.1) / 0.85

speed_guess = specific_speed / (np.sqrt(gas.mass_flow / gas.density_deltah(deltah)) * np.power(deltah, -0.75)) 
work_coeff_guess = 1 / np.power(specific_speed * specific_diameter, 2)
flow_coeff_guess = 1 / (specific_speed * np.power(specific_diameter, 3))

print(flow_coeff_guess, work_coeff_guess)


fig, axs = plt.subplots(2, 2, figsize=(6*1.2,4*1.2), sharex=True, sharey=True)
fig.tight_layout(w_pad=1.0, h_pad=1.5)
fig.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95)
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular'
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

mapw = axs[0,0].contourf(diff_coeff_range, speed_range, np.transpose(weights), locator=ticker.LogLocator(), levels=[10**-1, 10**-0.5, 1, 10**0.5, 10, 10**1.5, 10**2, 10**3, 10**4, 10**5, 10**6, 10**9])
divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mapw, cax=cax)

axs[0,0].scatter(diff_x, speed_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,0].annotate(txt, (diff_x[i]+0.02, speed_y[i]*0.7))


mape = axs[0,1].contourf(diff_coeff_range, speed_range, np.transpose(effs))
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mape, cax=cax)

axs[0,1].scatter(diff_x, speed_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,1].annotate(txt, (diff_x[i]+0.02, speed_y[i]*0.7))


map1 = axs[1,0].contourf(diff_coeff_range, speed_range, np.transpose(ESM_norm), locator=ticker.LogLocator(), levels=np.logspace(0, 6, 13))
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map1, cax=cax)

axs[1,0].scatter(diff_x, speed_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,0].annotate(txt, (diff_x[i]+0.02, speed_y[i]*0.7))

axs[1,0].plot([3, 3*0.8, 2.0, 15*0.8, 3],
            [speed_guess, speed_guess*0.8, speed_guess*0.8, 99e3, speed_guess], 
             '--o', c="k", label="Starting simplex")

axs[1,0].scatter([optimal.diffusion_ratio], [optimal.speed], color="k", marker="*", label="Optimised result")
#axs[1,0].plot([3], [speed_guess], "-o", color="tab:red", label="Specific speed/diameter guess")

#axs[1,0].contour(diff_coeff_range, speed_range, np.transpose(ESM_norm), levels=[3e5], c="k")
#axs[1,0].text(0.2, 0.2, "Low efficiency", rotation=80)
#axs[1,0].text(1.3, 0.25, "Low efficiency", rotation=-45)

axs[1,0].legend()

map2 = axs[1,1].contourf(diff_coeff_range, speed_range, np.transpose(ESM_cons), locator=ticker.LogLocator(), levels=np.logspace(0, 6, 13))
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map2, cax=cax)

#axs[1,1].contour(flow_coeff_range, work_coeff_range, np.transpose(ESM_cons), levels=[1.1e6], c="k")

axs[1,1].scatter(diff_x, speed_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,1].annotate(txt, (diff_x[i]+0.02, speed_y[i]*0.7))


axs[1,1].text(6, 100, "Low efficiency, large impeller", rotation=5)
axs[1,1].text(6, 70e3, "Negative radius increase")
#axs[1,1].text(1.2, 0.005, "Short blades", rotation=3)



axs[0,0].set_xticklabels(np.arange(0, 16, 1), usetex=True)
axs[0,1].set_xticklabels(np.arange(0, 16, 1), usetex=True)
axs[1,0].set_xticklabels(np.arange(0, 16, 1), usetex=True)
axs[1,1].set_xticklabels(np.arange(0, 16, 1), usetex=True)


axs[0,0].set_yticklabels([r"$10^1$", r"$10^2", r"$10^3", r"$10^4", r"$10^5"], usetex=True)
#axs[0,1].set_yticklabels([r"$10^1$", r"$10^2", r"$10^3", r"$10^4", r"$10^5"], usetex=True)
axs[1,0].set_yticklabels([r"$10^1$", r"$10^2", r"$10^3", r"$10^4", r"$10^5"], usetex=True)
#axs[1,1].set_yticklabels([r"$10^1$", r"$10^2", r"$10^3", r"$10^4", r"$10^5"], usetex=True)

axs[0,0].set_yscale("log")
axs[0,1].set_yscale("log")
axs[1,0].set_yscale("log")
axs[1,1].set_yscale("log")

#axs[0,0].set_xlabel(r"Diffusion ratio", usetex=True)
axs[0,0].set_ylabel(r"Speed \textit{RPM}", usetex=True)
axs[0,0].set_title(r"Weight at fixed flow/work coefficient", usetex=True)

#axs[0,1].set_xlabel(r"Diffusion ratio", usetex=True)
#axs[0,1].set_ylabel(r"Speed \textit{RPM}", usetex=True)
axs[0,1].set_title(r"Efficiency at fixed flow/work coefficient", usetex=True)


axs[1,0].set_xlabel(r"Diffusion ratio", usetex=True)
axs[1,0].set_ylabel(r"Speed \textit{RPM}", usetex=True)
axs[1,0].set_title(r"ESM at fixed flow/work coefficient", usetex=True)

axs[1,1].set_xlabel(r"Diffusion ratio", usetex=True)
#axs[1,1].set_ylabel(r"Speed \textit{RPM}", usetex=True)
axs[1,1].set_title(r"ESM at fixed flow/work coefficient with limits", usetex=True)

#plt.savefig("test.svg")
print("A")
plt.show()
"""
"""
### Convergence plotter

tolerance = [10000,1000,100,10,1,0.1,0.01,0.001,0.0001,0.00001]
error = [0.52533716,0.002422513,0.001171182,0.001164568,0.002026748,0.002026775,0.002026774,0.002026774,0.002026774,0.002026774]
steps = [18,56,125,155,708,758,805,827,842,862]

fig, ax1 = plt.subplots(figsize=(3.78, 2.83))
fig.subplots_adjust(left=0.15, right=0.85, bottom=0.15)

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":10})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

ax2 = ax1.twinx()

ax1.plot(tolerance, steps)
ax2.plot(tolerance, [i*100 for i in error], c="tab:orange")

ax1.set_xscale("log")
ax2.set_xscale("log")
ax1.set_yscale("log")
ax2.set_yscale("log")

ax1.set_xlim(tolerance[0], tolerance[-1])
ax2.set_xlim(tolerance[0], tolerance[-1])

ax1.tick_params(axis="y", labelcolor="tab:blue")
ax1.set_ylabel("Number of evaluations", color="tab:blue", usetex=True, fontsize=10)
ax2.tick_params(axis="y", labelcolor="tab:orange")

ax1.set_xlabel("Nelder-Mead tolerance", usetex=True, fontsize=10)
ax1.set_xticks([10**3, 10**1, 10**(-1), 10**(-3), 10**(-5)])
ax1.set_xticklabels([r"$10^3$", r"$10^1$", r"$10^{-1}$", r"$10^{-3}$", r"$10^{-5}$"], usetex=True, fontsize=10)
#ax2.set_ylabel("% error from lowest in search", color="tab:orange")

ax1.set_title("Nelder-Mead performance for axial compressor", usetex=True)

plt.show()

"""

"""

from component import *

testflow = GasFlow(29.71, 409.1, 287300, carbon_dioxide=False)
teststage = AxialStage(testflow, 1.581, 16042.3, 0.350, 0.484, 0.85)
teststage.estimate_efficiency()
teststage.estimate_ESM()
teststage.weight_estimate()

#print(teststage)
#print(teststage.R_hub_in)
#print(teststage.R_tip_in)

testflowR = GasFlow(13.61, 685.3, 1470700, carbon_dioxide=False)
testradial = RadialStage(testflowR, 2.042, 32650, 0.776, 0.35757, 0.7, 0.06985, 2.009)
testradial.estimate_efficiency()
testradial.weight_estimate(381e6, 2810)
testradial.estimate_ESM()

#print(testradial)
#print(testradial.R_mean_imp_exit)
#print(testradial.efficiency_impeller)


# Case from Cheng and Geld
#testflowH = GasFlow(0.4374, 343.15, 101325, carbon_dioxide=False)
#testHX = HeatExchangerAdvanced(6428.8, testflowH, 0.0768, 1, 2e-3, pitch_diameter_ratio=2, wall_thickness=2.65e-4, wall_thermal_conductivity=0.185, coolant="water", coolant_velocity=0.343)

#testHX.n_rows = 48
#testHX.n_tubes_wide = 42
#testHX.bulk_area = 0.0768
#testHX.inlet_area = 0.0768
#testHX.duct_length = 0.0793
#testHX.coolant_dT = 36.73

testflowH = GasFlow(0.15092175, 313.15, 101325, carbon_dioxide=False)
testHX = HeatExchangerAdvanced(2210, testflowH, 0.032467, 1, 3.2e-3, 2, 0.6e-3, 15, coolant="water", coolant_velocity=0.035)

testHX.n_rows = 38
testHX.n_tubes_wide = 39
testHX.duct_length = 0.124
testHX.coolant_dT = 27.81


h = testHX.heat_transfer_row() / testHX.coolant_dT
Nu_total = h * 6.4e-3 / testflowH.k
factor = 2* Nu_total / np.power(testHX.gasflow_in.Pr, 1/3)
print(factor)
print(testHX.Re_tubes)

testHX.gas_pressure_drop()
f = 2 * testflowH.density() * testHX.pressure_drop * 6.4e-3 / (0.124 * np.power(4*testflowH.mass_flow / testHX.bulk_area, 2))
print(f)
#print(2*testHX.gasflow_in.mass_flow / (testHX.bulk_area * testHX.gasflow_in.density()) )
#
"""

"""
### Performance map plotter for HX


from component import *
from design import *
from design2 import *
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl


gas = GasFlow(0.05, 300, 45000)

pitch_range =  np.linspace(1.1, 10, 500)
area_ratio_range = np.linspace(1, 10, 520)

weights = np.zeros(shape=(pitch_range.size, area_ratio_range.size))
pressure_drops = np.zeros(shape=(pitch_range.size, area_ratio_range.size))
ESM_norm = np.zeros(shape=(pitch_range.size, area_ratio_range.size))
ESM_cons = np.zeros(shape=(pitch_range.size, area_ratio_range.size))


optimal = generate_heat_exchanger(gas, 35000, 6.3e-3, 0.01, [1, .149, .121, 5.0])
print(optimal)
print(optimal.pitch_diameter_ratio)
print(optimal.bulk_area_ratio)
print(optimal.coolant_velocity)
print(optimal.pressure_drop)
print(optimal.n_tubes_wide)
print(optimal.inlet_size*1000)
print(optimal.duct_size*1000)
print(optimal.tube_diameter*optimal.n_tubes_wide*optimal.pitch_diameter_ratio*1000)
print(optimal.diffuser_length*1000)
print(optimal.duct_length*1000)
print(optimal.nozzle_length*1000)
print(optimal.duct_length*1000)
print(optimal)


# Labelled points
pitch_x = [optimal.pitch_diameter_ratio, 4, optimal.pitch_diameter_ratio, optimal.pitch_diameter_ratio]
area_ratio_y = [optimal.bulk_area_ratio, optimal.bulk_area_ratio, 9, 1.5]
label_n = ["I", "II", "III", "IV"]


bad = HeatExchangerAdvanced(35000*0.05, gas, 6.3e-3, area_ratio_y[3], pitch_diameter_ratio=pitch_x[3], coolant_velocity=optimal.coolant_velocity)
bad.gas_pressure_drop()
bad.coolant_pumping_power()
bad.weight_estimate()
bad.estimate_ESM()
print(bad)
print(bad.pitch_diameter_ratio)
print(bad.bulk_area_ratio)
print(bad.coolant_velocity)
print(bad.pressure_drop)
print(bad.n_tubes_wide)
print(bad.inlet_size*1000)
print(bad.duct_size*1000)
print(bad.diffuser_length*1000)
print(bad.duct_length*1000)
print(bad.nozzle_length*1000)


for pitch_i, pitch in enumerate(pitch_range):
    for area_ratio_i, area_ratio in enumerate(area_ratio_range):
        esm_cons = make_heat_exchanger_pd((area_ratio, optimal.coolant_velocity, pitch-1), 2075, gas, 6.3e-3, 0.01, [1, .149, .121, 5.0], "ammonia", output=False)
        esm_nocons = make_heat_exchanger_pd((area_ratio, optimal.coolant_velocity, pitch-1), 2075, gas, 6.3e-3, 0.01, [1, .149, .121, 5.0], "ammonia", output=True)

        ESM_norm[pitch_i, area_ratio_i] = esm_nocons
        ESM_cons[pitch_i, area_ratio_i] = esm_cons

        HX = HeatExchangerAdvanced(35000*0.05, gas, 6.3e-3, area_ratio, pitch_diameter_ratio=pitch, coolant_velocity=optimal.coolant_velocity)

        HX.gas_pressure_drop()
        HX.coolant_pumping_power()
        HX.weight_estimate()
        HX.estimate_ESM()

        weights[pitch_i, area_ratio_i] = HX.weight
        pressure_drops[pitch_i, area_ratio_i] = HX.pressure_drop
    
    if pitch_i%10==0:
        print(100*pitch_i/pitch_range.size, "%")


pitch_guess = 1.5
area_guess = 1.5


fig, axs = plt.subplots(2, 2, figsize=(6*1.2,4*1.2), sharex=True, sharey=True)
fig.tight_layout(w_pad=1.0, h_pad=1.5)
fig.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95)
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular'
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

mapw = axs[0,0].contourf(pitch_range, area_ratio_range, np.transpose(weights))#, locator=ticker.LogLocator(), levels=[10**-1, 10**-0.5, 1, 10**0.5, 10, 10**1.5, 10**2, 10**3, 10**4, 10**5, 10**6, 10**9])
divider = make_axes_locatable(axs[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mapw, cax=cax)

axs[0,0].scatter(pitch_x, area_ratio_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,0].annotate(txt, (pitch_x[i], area_ratio_y[i]-0.5))

print(np.amin(weights), np.amax(weights))



mappd = axs[0,1].contourf(pitch_range, area_ratio_range, np.transpose(pressure_drops), locator=ticker.LogLocator(), levels=np.logspace(-1, 5, 13))
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(mappd, cax=cax)

axs[0,1].scatter(pitch_x, area_ratio_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[0,1].annotate(txt, (pitch_x[i], area_ratio_y[i]-0.5))




map1 = axs[1,0].contourf(pitch_range, area_ratio_range, np.transpose(ESM_norm), levels=np.linspace(200, 400, 9))
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map1, cax=cax)

axs[1,0].scatter(pitch_x, area_ratio_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,0].annotate(txt, (pitch_x[i], area_ratio_y[i]-0.5))

axs[1,0].plot([pitch_guess, pitch_guess*0.8, 50, pitch_guess*0.8, pitch_guess],
            [area_guess, area_guess*0.8, area_guess*0.8, 10, area_guess], 
             '--o', c="k", label="Starting simplex")

axs[1,0].scatter([optimal.pitch_diameter_ratio], [optimal.bulk_area_ratio], color="k", marker="*", label="Optimised result")
#axs[1,0].plot([pitch_guess], [area_guess], "-o", color="tab:red", label="Initial guess")

#axs[1,0].contour(diff_coeff_range, speed_range, np.transpose(ESM_norm), levels=[3e5], c="k")

#axs[1,0].text(1.3, 0.25, "Low efficiency", rotation=-45)

axs[1,0].legend()




map2 = axs[1,1].contourf(pitch_range, area_ratio_range, np.transpose(ESM_cons), levels=np.linspace(200, 400, 9))
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(map2, cax=cax)

axs[1,1].contour(pitch_range, area_ratio_range, np.transpose(ESM_cons), levels=[1000], c="k")

axs[1,1].scatter(pitch_x, area_ratio_y, marker="+", c="k")
for i, txt in enumerate(label_n):
    axs[1,1].annotate(txt, (pitch_x[i], area_ratio_y[i]-0.5))



axs[1,1].text(1.1, 1, "High pressure drop", rotation=-30)
axs[1,1].text(6, 1.1, "Long duct")
##axs[1,1].text(1.2, 0.005, "Short blades", rotation=3)

axs[0,0].set_xlim(1, 10)
axs[0,1].set_xlim(1, 10)
axs[1,0].set_xlim(1, 10)
axs[1,1].set_xlim(1, 10)


#axs[0,0].set_xticklabels(np.arange(1, 11, 1), usetex=True)
#axs[0,1].set_xticklabels(np.arange(1, 11, 1), usetex=True)
axs[1,0].set_xticklabels(np.arange(1, 11, 1), usetex=True)
axs[1,1].set_xticklabels(np.arange(1, 11, 1), usetex=True)

axs[0,0].set_yticklabels(np.arange(1, 11, 1), usetex=True)
#axs[0,1].set_yticklabels(np.arange(1, 11, 1), usetex=True)
axs[1,0].set_yticklabels(np.arange(1, 11, 1), usetex=True)
#axs[1,1].set_yticklabels(np.arange(1, 11, 1), usetex=True)



#axs[0,0].set_xlabel(r"Pitch-diameter ratio", usetex=True)
axs[0,0].set_ylabel(r"Inlet-duct area ratio", usetex=True)
axs[0,0].set_title(r"Weight at fixed coolant flow", usetex=True)

#axs[0,1].set_xlabel(r"Pitch-diameter ratio", usetex=True)
#axs[0,1].set_ylabel(r"Inlet-duct area ratio", usetex=True)
axs[0,1].set_title(r"Pressure drop (Pa) at fixed coolant flow", usetex=True)


axs[1,0].set_xlabel(r"Pitch-diameter ratio", usetex=True)
axs[1,0].set_ylabel(r"Inlet-duct area ratio", usetex=True)
axs[1,0].set_title(r"ESM at fixed coolant flow", usetex=True)

axs[1,1].set_xlabel(r"Pitch-diameter ratio", usetex=True)
#axs[1,1].set_ylabel(r"Inlet-duct area ratio", usetex=True)
axs[1,1].set_title(r"ESM at fixed coolant flow with limits", usetex=True)

#plt.savefig("test.svg")
print("A")
plt.show()
"""

"""
import matplotlib.pyplot as plt

number_of_intercoolers = [10, 11, 12, 13, 14, 15, 20, 25, 30]
mass = [
    253.9, 233.8, 245.1, 292.3, 277.4, 228.7, 269.69, 342.69, 378.59]
ESM_cooling_power = [
    3879.6, 3839.0, 3800.5, 3800.9, 3816.8, 3799.2, 4057.1, 4076.0, 4034.8 
]

num_centris = [4, 5, 6, 6, 7, 8, 12, 14, 14]
num_scrolls = [6, 6, 6, 7, 7, 7, 8, 11, 10]

total_ESM = []
for i in range(len(number_of_intercoolers)):
    total_ESM.append(mass[i]+ESM_cooling_power[i])

fig, ax1 = plt.subplots(figsize=(2.96, 2.22))
fig.subplots_adjust(left=0.175, right=0.825, bottom=0.175)

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":8
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

ax2 = ax1.twinx()
l1 = ax2.scatter(number_of_intercoolers, mass, color="tab:green", label="Mass")
l2 = ax1.scatter(number_of_intercoolers, ESM_cooling_power, color="tab:blue", label="Power/cooling ESM")
l3 = ax1.scatter(number_of_intercoolers, total_ESM, color="tab:blue", label="Total ESM", marker="s", edgecolors='black')

ax1.set_xlabel('Number of intercoolers ', usetex=True, fontsize=8)

ax1.set_ylabel('ESM (kg)', color="tab:blue", usetex=True, fontsize=8)
ax2.set_ylabel('Mass (kg)', color="tab:green", usetex=True, fontsize=8)

ax1.set_yticklabels(np.arange(0, 4500, 500), color="tab:blue", usetex=True, fontsize=8)
ax2.set_yticklabels(np.arange(0, 450, 50), color="tab:green", usetex=True, fontsize=8)

ax1.set_xticks([10, 15, 20, 25, 30])
ax1.set_xticklabels(["10","15","20","25","30"], usetex=True, fontsize=8)
ax2.set_ylim(0, 400)

ax1.legend(handles=[l1,l2,l3],handlelength=1, fontsize=8)
plt.title("ESM variation with intercooler number", usetex=True, fontsize=12)



fig, ax1 = plt.subplots(figsize=(2.96, 2.22))
ax2 = ax1.twinx()
fig.subplots_adjust(left=0.175, right=0.825, bottom=0.175)

l2 = ax1.scatter(number_of_intercoolers, num_scrolls, label="\# scrolls", color="tab:blue", marker="s", edgecolors='black')
l1 = ax1.scatter(number_of_intercoolers, num_centris, label="\# centrifugals", color="tab:blue")


l3 = ax2.scatter(number_of_intercoolers, mass, color="tab:green", label="Mass")

ax1.set_xlabel('Number of intercoolers ', usetex=True, fontsize=8)
ax1.set_ylabel('Number of components', color="tab:blue", usetex=True, fontsize=8)
ax2.set_ylabel('Mass (kg)', color="tab:green", usetex=True, fontsize=8)

ax1.set_yticklabels(np.arange(4, 17, 2), color="tab:blue", usetex=True, fontsize=8)
ax2.set_yticklabels(np.arange(0, 450, 50), color="tab:green", usetex=True, fontsize=8)

ax1.set_xticks([10, 15, 20, 25, 30])
ax1.set_xticklabels(["10","15","20","25","30"], usetex=True, fontsize=8)


ax1.legend(handles=[l1,l2,l3],handlelength=1, fontsize=8)
plt.title("Compressor type with intercooler number", usetex=True, fontsize=12)

plt.show()

"""

"""


import matplotlib.pyplot as plt

mass_flow = [11.5, 22.5, 33.8, 45.1, 90.2]

mass = [
    116.8, 228.691, 400.691, 371.791, 675.891]
ESM_power = [
    1851.321241, 1798.706984, 2038.030351, 1730.807713, 1710.800243]
ESM_cooling = [
    760.5482588, 2310.462754, 3413.744021, 4430.514777, 8740.08197]


total_ESM = []
for i in range(len(mass_flow)):
    total_ESM.append(mass[i]+ESM_power[i]+ESM_cooling[i])

fig, ax1 = plt.subplots(figsize=(2.96, 2.22))
fig.subplots_adjust(left=0.19, right=0.84, bottom=0.175)

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":8
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

ax2 = ax1.twinx()

l1, = ax1.plot(mass_flow, total_ESM, color="tab:blue", label="Total ESM")
l2, = ax1.plot(mass_flow, ESM_power, color="tab:blue", label="Power ESM", linestyle="--")
l3, = ax1.plot(mass_flow, ESM_cooling, color="tab:blue", label="Cooling ESM", linestyle=":")

l4, = ax2.plot(mass_flow, mass, color="tab:green", label="Mass")

ax1.set_xlabel('Mass flow (g/s) ', usetex=True, fontsize=8)

ax1.set_ylabel('ESM (kg)', color="tab:blue", usetex=True, fontsize=8)
ax2.set_ylabel('Mass (kg)', color="tab:green", usetex=True, fontsize=8)

ax1.set_yticklabels(np.arange(0, 13000, 2000), color="tab:blue", usetex=True, fontsize=8)
ax1.set_ylim(0, 12000)

ax2.set_yticklabels(np.arange(0, 1000, 100), color="tab:green", usetex=True, fontsize=8)
ax2.set_ylim(0, 900)

ax1.set_xticks([0, 20, 40, 60, 80, 100])
ax1.set_xticklabels(["0","20","40","60","80", "100"], usetex=True, fontsize=8)
#ax2.set_ylim(0, 400)

ax1.legend(handles=[l1,l2,l3, l4])
plt.title("ESM variation with mass flow", usetex=True, fontsize=12)
plt.show()

"""


"""
from component import *
from design import *
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl


plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":10
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')
fig, ax = plt.subplots(figsize=(3.8, 2.8))
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.925)

gas = GasFlow(0.05, 220, 2000)


optimal = optimise_axial(gas, 1.1, 0.85, [1, .149, .121, 5.0])
optimal_speed = optimal.speed
print(optimal)
print(optimal.flow_coeff, optimal.work_coeff, optimal.speed)

bad = AxialStage(gas, 1.1, 1e3, optimal.work_coeff, optimal.flow_coeff, 0.85)
bad.estimate_efficiency()
bad.weight_estimate()
bad.estimate_ESM()
print(bad)


speed_range = np.logspace(2.5, 5, 1500)
ESM_cons_range = np.zeros(speed_range.size)
ESM_norms_range = np.zeros(speed_range.size)

for speed_i, speed in enumerate(speed_range):
    ESM_cons_range[speed_i] = make_axial((speed, optimal.work_coeff, optimal.flow_coeff), gas, 1.1, 0.85, [1, .149, .121, 5.0])
    ESM_norms_range[speed_i] = make_axial((speed, optimal.work_coeff, optimal.flow_coeff), gas, 1.1, 0.85, [1, .149, .121, 5.0], output=True)


plt.plot(speed_range, ESM_norms_range, label="ESM")
plt.plot(speed_range, ESM_cons_range, label="ESM with limits")
plt.scatter([optimal.speed], [optimal.ESM], c="k", marker="*", s=100, label="Optimum")

plt.annotate(r"Bladeheight  $\le 3\,$mm", [2.1e3, 150], [2.5e3, 175], 
             arrowprops=dict(arrowstyle="->", connectionstyle="angle, angleA=-90, angleB=180, rad=5"), usetex=True)

plt.annotate(r"Hub radius  $\le 0$", [3.5e4, 100], [5e3, 125], 
             arrowprops=dict(arrowstyle="->", connectionstyle="angle, angleA=-90, angleB=180, rad=5"), usetex=True)

plt.ylim(0, 200)
plt.xscale("log")
plt.legend(loc="lower left")
plt.xlabel("Speed (RPM)", usetex=True)
plt.ylabel("ESM (kg)", usetex=True)
plt.title("Axial compressor ESM function with limits", usetex=True, fontsize=12)
plt.show()
"""

"""
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(figsize=(2.83, 3.77))
fig.subplots_adjust(left=0.2, right=0.95, bottom=0.125, top=0.925)

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":10})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

plt.xlim(0,8)
plt.ylim(0.2, 4.0)
plt.yscale("log")

plt.xticks([0,1,2,3,4,5,6,7,8])
plt.yticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0], labels=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0])

plt.grid(visible=True)
plt.xlabel(r"Flow coefficient $\Phi$")
plt.ylabel(r"Work coefficient $\Psi$")
plt.title("Smyth and Miller contours")

plt.show()
"""
"""

import numpy as np

areas = np.linspace(1, 15)

L_W_diffuser = 0.6643 * np.power(areas, 2.7023) 
efficiency_diffuser = 1 - (0.7328 * np.power(areas, -0.997))

fig, ax1 = plt.subplots(figsize=(2.96, 2.22))
fig.subplots_adjust(left=0.2, right=0.84, bottom=0.2)

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":8
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

ax2 = ax1.twinx()

l1 = ax1.plot(areas, L_W_diffuser, c="tab:blue")
l2 = ax2.plot(areas, efficiency_diffuser, c="tab:orange")


ax1.set_xlabel('Diffusion ratio', usetex=True, fontsize=10)

ax1.set_ylabel('Length ratio', color="tab:blue", usetex=True, fontsize=10)
ax2.set_ylabel('Efficiency', color="tab:orange", usetex=True, fontsize=10)

ax1.set_yticklabels(np.arange(0, 1100, 200), color="tab:blue", usetex=True, fontsize=10)
ax1.set_ylim(0, 1100)

ax2.set_yticks(np.arange(0.3, 1.1, 0.1))
ax2.set_yticklabels(labels=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], color="tab:orange", usetex=True, fontsize=10)
ax2.set_ylim(0.3, 1.0)

#ax1.set_xticks([0, 20, 40, 60, 80, 100])
#ax1.set_xticklabels(["0","20","40","60","80", "100"], usetex=True, fontsize=8)
#ax2.set_ylim(0, 400)

plt.title("Diffuser performance", usetex=True, fontsize=12)
plt.show()

"""