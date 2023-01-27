from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

from component import GasFlow, HeatExchangerAdvanced
from design import optimise_hx_anypd



T_start = 190
T_end = 260

s_start = PropsSI("Smass", "T|gas", T_start, "P", 900, "CO2")
s_end = PropsSI("Smass", "T|gas", T_end, "P", 550e3, "CO2")

n_thermo_points = 5

s_range = np.linspace(s_start, s_end, n_thermo_points)
T_range = np.linspace(T_start, T_end, n_thermo_points)
p_range = np.array([PropsSI("P", "T|gas", T_range[i], "Smass", s_range[i], "CO2") for i in range(n_thermo_points)])
d_range = np.array([PropsSI("D", "T|gas", T_range[i], "Smass", s_range[i], "CO2") for i in range(n_thermo_points)])

m_dot_range = np.logspace(-3, -1, 4)*5

dT_range = np.linspace(5, 50, 10)

PD_range = np.linspace(0.01, 0.1, 5)

all_dTs = []
all_PDs = []


for T_index, T in enumerate(T_range):
    
    for m_index, m_dot in enumerate(m_dot_range):
        gasflow = GasFlow(m_dot, T, p_range[T_index])
        flow_velocity_in = 15
        flow_area_in = m_dot / (gasflow.density() * flow_velocity_in)
        for dT_index, dT in enumerate(dT_range):
            q = m_dot * dT * 870
            for PD_index, PD in enumerate(PD_range):
                HX = optimise_hx_anypd(gasflow, q, flow_area_in, PD*p_range[T_index])

                if (HX.duct_length + HX.nozzle_length + HX.diffuser_length) < 3:
                    all_dTs.append(dT)
                    all_PDs.append((HX.pressure_drop / p_range[T_index]))
        
    print(T, m_dot)

plt.scatter(all_dTs, all_PDs)
plt.show()