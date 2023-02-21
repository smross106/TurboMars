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
HX_ESM = []

HX_PR_limit = 0.2
HX_PR_increase_ratio = 1.2
HX_PR_start = 0.01
HX_PR_increases_limit = 18

total_compression_work = 0
total_compressor_weight= 0
compressor_efficiency = []

total_ESM = 0

mass_flow = 0.1
HX_max_volume = max(2, 2*(mass_flow/0.05)**3)

N_stages = 25
rad_PR = 1.35

start = time.time()

for i in range(N_stages):
    flow_cold = component.GasFlow(mass_flow, T_in, p_in)

    
    # ax1 = design.optimise_axial(flow, flow.delta_h_PR(1.1), 0.85)

    rad1 = design.optimise_radial(flow_cold, rad_PR, 0.85)

    print("Radial",str(i+1))
    #print(rad1)
    #print(rad1.efficiency_delta_mach(0.85), rad1.W_mean_in, rad1.V_theta_imp_exit )
    #print(2000*rad1.R_mean_inlet, 2000*rad1.R_mean_imp_exit, 2000*rad1.R_stator_exit, 1000*rad1.bladeheight_imp_exit, rad1.speed, rad1.efficiency, rad1.weight, rad1.power)

    hot_T = flow_cold.temperature * ((np.power(rad_PR, flow_cold.gm1/flow_cold.gamma) - 1)/rad1.efficiency + 1)


    flow_hot = component.GasFlow(mass_flow, hot_T, flow_cold.pressure*rad_PR)
    Q = flow_hot.delta_h_delta_T(abs(250-flow_hot.temperature)) 

    total_compression_work += Q * flow_hot.mass_flow
    total_compressor_weight += rad1.weight
    compressor_efficiency.append(rad1.efficiency)
    total_ESM += rad1.ESM

    A_HX_in = flow_hot.mass_flow / (flow_hot.density() * 10)


    HX_PR_target = HX_PR_start
    HX_PR_increases = 0
    hx1_PR_actual = 1
    HX_approved = False

    while ((hx1_PR_actual > HX_PR_target) or (not HX_approved)) and (HX_PR_increases < HX_PR_increases_limit):
        if i < 3 or rad1.A_stator_exit < 6.25e-4:
            hx1 = design.optimise_hx_pd(flow_hot, Q, A_HX_in, HX_PR_target)
        else:
            hx1 = design.optimise_hx_pd(flow_hot, Q, rad1.A_stator_exit, HX_PR_target)
        hx1_PR_actual = hx1.pressure_drop / flow_hot.pressure

        HX_volume = hx1.total_length * (hx1.duct_size**2)
        if HX_volume < HX_max_volume:
            HX_approved = True

        HX_PR_target *= HX_PR_increase_ratio
        HX_PR_increases += 1
        print("HX number ", str(i+1),", iteration ",str(HX_PR_increases))

    #print(hx1)
    #print(hx1.coolant_velocity)
    HX_ESM.append(hx1.weight)

    total_ESM += hx1.ESM

    HX_PR = hx1.pressure_drop / flow_hot.pressure

    print("HX",str(i+1), HX_PR)


    if HX_PR > 0.99:
        HX_PR = 0.99

    """for PR in HX_PR_range:
        hx1 = design.optimise_hx_pd(flow_hot, flow_hot.delta_h_delta_T(abs(250-flow_hot.temperature)), rad1.A_stator_exit, PR)
        real_HX_PR = hx1.pressure_drop / flow_hot.pressure
        if hx1.ESM < best_HX_ESM:
            best_HX_ESM = hx1.ESM + (real_HX_PR * 160000 * flow_hot.mass_flow)
            best_hx = hx1
            HX_PR = PR
    
    print(best_hx)
    HX_ESM.append(best_hx.ESM)"""

    #HX_PR = hx1.pressure_drop / flow_hot.pressure
    """if HX_PR >= 1:
        HX_PR = 0.99
    HX_length.append((hx1.duct_length+hx1.diffuser_length+hx1.nozzle_length))
    HX_weight.append(hx1.weight)
    HX_duct_length.append(hx1.duct_length)"""
    #print(1-HX_PR, hx1.cooling_power, hx1.total_length*1000, hx1.duct_length*1000, hx1.pumping_power, hx1.n_tubes, 1000*hx1.duct_size, hx1.weight)
    
    T_in = 250
    p_in = flow_cold.pressure * rad_PR * (1 - HX_PR)


#print(sum(HX_ESM))
print("")
print("power", total_compression_work)
print("weight", total_compressor_weight)
print("average efficiency", sum(compressor_efficiency)/len(compressor_efficiency))
print("total ESM", total_ESM)

print(time.time() - start)

print(compressor_efficiency)
print(p_in/1000)

"""fig, ax = plt.subplots()
plt.plot(HX_duct_length, label="Length of HX duct")
plt.plot(HX_length, label="Total HX length")
plt.ylabel("Length (m)")
plt.xlabel("Heat exchanger number")
plt.title("Heat exchanger lengths with 50 bar liquid ammonia")
plt.legend()
"""
"""fig, ax = plt.subplots()
plt.plot(HX_weight)
plt.ylabel("Weight (kg)")
plt.xlabel("Heat exchanger number")
plt.title("Heat exchanger weight with VLT")
plt.show()"""
"""
print("total weight", str(sum(HX_weight)))
print("max weight", str(max(HX_weight)))

print("total length", str(sum(HX_length)))
print("max length", str(max(HX_length)))

plt.show()"""

#print(p_in)
"""p_compress_out = p_in
T_compress_out = 250
compressed_flow = component.GasFlow(0.05, T_compress_out, p_compress_out)
hx_cond = design.optimise_hx_pd(compressed_flow, 72.4e3, hx1.inlet_area, 0.02)
print(hx_cond)
print(hx_cond.duct_length)"""

"""
print(hx1.pitch_diameter_ratio*hx1.tube_diameter*1000)
print(hx1.bulk_area, , hx1.pitch_diameter_ratio)
print(hx1.required_area)
print(hx1.n_tubes_wide)"""

