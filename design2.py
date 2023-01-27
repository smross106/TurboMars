import component

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import time
import design

T_in = 222.88
p_in = 1594.4


HX_weight = []
HX_length = []
HX_duct_length = []
HX_PR_range = np.linspace(0.01, 0.2, 20)

for i in range(10):
    flow_cold = component.GasFlow(0.05, T_in, p_in)

    PR = 2.0
    # ax1 = design.optimise_axial(flow, flow.delta_h_PR(1.1), 0.85)

    rad1 = design.optimise_radial(flow_cold, flow_cold.delta_h_PR(PR), 0.85)

    print("Radial",str(i+1))
    #print(2000*rad1.R_mean_inlet, 2000*rad1.R_mean_imp_exit, 2000*rad1.R_stator_exit, 1000*rad1.bladeheight_imp_exit, rad1.speed, rad1.efficiency, rad1.weight, rad1.power)

    
    eff = 0.88

    hot_T = flow_cold.temperature * ((np.power(PR, flow_cold.gm1/flow_cold.gamma) - 1)/rad1.efficiency + 1)

    flow_hot = component.GasFlow(0.05, hot_T, flow_cold.pressure*PR)

    best_HX_ESM = 1e100
    HX_PR = 0
    best_hx = None
    for PR in HX_PR_range:
        hx1 = design.optimise_hx_anypd(flow_hot, flow_hot.delta_h_delta_T(abs(250-flow_hot.temperature)), rad1.A_stator_exit, PD)
        if hx1.ESM < best_HX_ESM:
            best_HX_ESM = hx1.ESM
            best_hx = hx1
            HX_PR = PR
    print("HX",str(i+1))
    print(hx1)

    #HX_PR = hx1.pressure_drop / flow_hot.pressure
    if HX_PR >= 1:
        HX_PR = 0.99
    HX_length.append((hx1.duct_length+hx1.diffuser_length+hx1.nozzle_length))
    HX_weight.append(hx1.weight)
    HX_duct_length.append(hx1.duct_length)
    #print(hx1.cooling_power, (hx1.duct_length+hx1.diffuser_length+hx1.nozzle_length), hx1.duct_length, hx1.pumping_power, hx1.n_tubes, 1000*hx1.duct_size, hx1.weight)
    
    T_in = 250
    p_in = p_in * PR * (1 - HX_PR)

"""fig, ax = plt.subplots()
plt.plot(HX_duct_length, label="Length of HX duct")
plt.plot(HX_length, label="Total HX length")
plt.ylabel("Length (m)")
plt.xlabel("Heat exchanger number")
plt.title("Heat exchanger lengths with 50 bar liquid ammonia")
plt.legend()
"""
fig, ax = plt.subplots()
plt.plot(HX_weight)
plt.ylabel("Weight (kg)")
plt.xlabel("Heat exchanger number")
plt.title("Heat exchanger weight with VLT")
plt.show()
"""
print("total weight", str(sum(HX_weight)))
print("max weight", str(max(HX_weight)))

print("total length", str(sum(HX_length)))
print("max length", str(max(HX_length)))

plt.show()"""

"""print(p_in)
p_compress_out = p_in
T_compress_out = 250
compressed_flow = component.GasFlow(0.05, T_compress_out, p_compress_out)
hx_cond = design.optimise_hx(compressed_flow, 72.4e3, hx1.inlet_area, compressed_flow.pressure*0.06)
print(hx_cond)
print(hx_cond.duct_length)"""

"""
print(hx1.pitch_diameter_ratio*hx1.tube_diameter*1000)
print(hx1.bulk_area, , hx1.pitch_diameter_ratio)
print(hx1.required_area)
print(hx1.n_tubes_wide)"""

