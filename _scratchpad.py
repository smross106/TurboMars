"""import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6,4))

flow_rates_cryo = [0.0624, 0.458333333, 0.06, 0.0125, 0.5, 0.088, 0.02, 0.619]
specific_works_cryo = [28846.15385, 3173.236364, 24000, 3135.6, 4435.2, 19636.36364, 7200, 12969.30533]

flow_rates_other = [0.041666667, 0.458333333, 0.08, 0.5, 0.06732, 0.083]
specific_works_other = [4320, 5498.181818, 5175, 12888, 12834.2246, 5204.819277]

plt.scatter(flow_rates_cryo, specific_works_cryo, label="Cryocooler")
plt.scatter(flow_rates_other, specific_works_other, label="Other")

plt.legend()
plt.grid(visible=True)
plt.xlim(0, 0.75)
plt.xlabel("Flow rate, kg/hr")
plt.ylabel("Specific work, kJ/kg CO2")
plt.title("Existing ISRU Atmosphere Extraction Systems")
plt.show()"""


import CoolProp.CoolProp as CP
import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
from component import HeatExchangerAdvanced, GasFlow

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

fig, ax = plt.subplots()
print(np.amax(L_D_0_1))
print(np.amin(L_D_0_1))
cs = plt.contourf(T_range, v_range, np.transpose(L_D_0_1),  locator=ticker.LogLocator(numticks=27))
plt.colorbar()
plt.xlabel("Gas temperature (K)")
plt.ylabel("Gas velocity in cylinder (m/s)")
plt.title("Piston L, D=10mm")

"""p = cs.collections[3].get_paths()[0]

v = p.vertices
x = v[:,0]
y = v[:,1]

print(x)
print(y)

fig, ax = plt.subplots()
plt.plot(x,y)"""

fig, ax = plt.subplots()
plt.plot(T_range, d_range)

plt.show()