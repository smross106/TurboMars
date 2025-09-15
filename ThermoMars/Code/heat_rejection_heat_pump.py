import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

"""
References:
 - [1] - A solar azimuth formula that renders circumstantial treatment unnecessary without compromising mathematical rigor: 
        Mathematical setup, application and extension of a formula based on the subsolar point and atan2 function, Zhang et al 2021
 - [2] - Über die Extinktion des Lichtes in der Erdatmosphäre, Schoenberg 1929
 - [3] - Mars Fact Sheet, https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
 - [4] - Thermal control of MSL Rover "Curiosity" using an Active Fluid Loop, Birur 2013 
 - [5] - Advanced Active Thermal Control Systems Architecture Study, Hanford and Ewert
 - [6] - Electric feed systems for liquid propellant rockets, Rachov et al
 - [7] - Mass Optimization and Experimental Validation of a Forced-Convection Heat Exchanger for Mars Surface Waste 
        Heat Rejection, Colgan 2023

"""

heat_exchanger_power_density = 2.72/1000 # kg/W
heat_pump_efficiency = 0.50
# From [5]
REL_heat_exchanger_power_density = heat_exchanger_power_density/100

def compressor_mass(power_kw):
    mass = 31.83 * np.power(power_kw, 0.476)

    return(mass)


def ESM_radiator(T_reject):
    # kg/W
    delta_T = T_reject - 235
    if delta_T < 0:
        return(np.nan)
    one_over_ESM = (5.858E-6 * np.power(delta_T, 2)) - (1.004E-4 * delta_T) + (1.498E-2)
    ESM = 1/one_over_ESM

    return(ESM/1000)

def ESM_FCHX(T_reject):
    # kg/W
    delta_T = T_reject - 235
    if delta_T < 0:
        return(np.nan)
    one_over_ESM = 2.853E-2 * np.log(delta_T) + 7.468E-2
    #one_over_ESM = (1.7225E-8 * np.power(T_reject, 3)) - (2.141E-5 * np.power(T_reject, 2)) + (8.863E-3 * T_reject) - 1.0003E0
    ESM = 1/one_over_ESM

    return(ESM/1000)

def increase_T(T_input, T_output, heat_in, power_ESM = 0.1, use_REL=False, eta=0.5):
    carnot = T_output / (T_output-T_input)
    
    CoP = carnot * eta

    work_in = heat_in / CoP
    heat_out = heat_in + work_in

    compressor_weight = compressor_mass(work_in/1000)
    if not use_REL:
        hx_in_weight = heat_in * heat_exchanger_power_density
        hx_out_weight = heat_out * heat_exchanger_power_density
    else:
        hx_in_weight = heat_in * REL_heat_exchanger_power_density
        hx_out_weight = heat_out * REL_heat_exchanger_power_density
    
    return(compressor_weight+hx_in_weight+hx_out_weight+work_in*power_ESM, heat_out, work_in, CoP)

print(ESM_radiator(290)*10000)
exit()
print(increase_T(260, 280, 100e3, 0, True, eta=0.4))
print(ESM_FCHX(280)*179e3)

exit()

T_start = 240
T_end = 350
Q_in = 100E3

W_start = ESM_radiator(T_start)*Q_in

plt.scatter(T_start, W_start, c="tab:blue", label="Radiator")

W_start = ESM_FCHX(T_start)*Q_in

plt.scatter(T_start, W_start, c="tab:orange", label="FCHX")

etas = [0.5, 0.6, 0.7]
radiators = [[],[],[]]
fchxs = [[],[],[]]
temps = []

for T_end in range(245, 505, 5):
    temps.append(T_end)
    for eta_i, hp_eta in enumerate(etas):
        W_pump, Q_out = increase_T(T_start, T_end, Q_in, eta=hp_eta, use_REL=True)
        W_end = W_pump + ESM_radiator(T_end)*Q_out
        radiators[eta_i].append(W_end)


        W_pump, Q_out = increase_T(T_start, T_end, Q_in, eta=hp_eta, use_REL=True)
        W_end = W_pump + ESM_FCHX(T_end)*Q_out
        fchxs[eta_i].append(W_end)

width=2


plt.plot(temps, radiators[0], c="tab:blue", linestyle="solid", lw=width)
plt.plot(temps, fchxs[0], c="tab:orange", linestyle="solid", lw=width)

plt.plot(temps, radiators[1], c="tab:blue", linestyle="dashed", lw=width)
plt.plot(temps, fchxs[1], c="tab:orange", linestyle="dashed", lw=width)

plt.plot(temps, radiators[2], c="tab:blue", linestyle="dotted", lw=width)
plt.plot(temps, fchxs[2], c="tab:orange", linestyle="dotted", lw=width)

legend_elements = [Line2D([0], [0], marker='o', color='w',markerfacecolor='tab:blue', label='Radiator', markersize=10),
                   Line2D([0], [0], marker='o', color='w',markerfacecolor='tab:orange', label='FCHX', markersize=10),
                   Line2D([0], [0], color='k', lw=2, label='$\eta$ = 50%', linestyle="solid"),
                   Line2D([0], [0], color='k', lw=2, label='$\eta$ = 60%', linestyle="dashed"),
                   Line2D([0], [0], color='k', lw=2, label='$\eta$ = 70%', linestyle="dotted")
]

plt.ylabel("Heat rejection system mass, kg")
plt.xlabel("Rejection temperature, K")
plt.title("Effect of heat pumping on heat rejection system mass \n 100kWth, power 100kg/kWe, REL HX weight, 300K cold side")
plt.legend(handles=legend_elements)


plt.show()