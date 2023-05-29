import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, interpolate

sigma = 5.67e-8

#def T_sky(T_air):
#    T_sky = 0.0552 * np.power(T_air, 1.5)
#    return(T_sky)

# Material constants
emissivity = 0.94
absorbivity = 0.15

T_air = 185
T_sky = 0.0552 * np.power(T_air, 1.5)
h_pipe = 100
start_fluid_temps = [350]
end_fluid_temp = 190

# Thermal inertia = mdot c_p
# Value is for 2kJ/kgK, 0.5kg/s
thermal_inertia = 1000

fluid_temps = []

def radiative_loss(T_wall, T_fluid):
    loss = sigma * (emissivity * (T_wall**4) - absorbivity * (0.5 * T_air**4 + 0.5 * T_sky**4))
    return(loss) 

def wall_flux_root(T_wall, T_fluid):
    # Passed into scipy rootfinder
    val = (radiative_loss(T_wall, T_fluid) - (h_pipe * (T_fluid - T_wall)))
    return(val)

for sft_index, start_fluid_temp in enumerate(start_fluid_temps):
    fluid_temps.append([start_fluid_temp])

    while fluid_temps[sft_index][-1] >= end_fluid_temp:
    #for i in range(1):
        T_fluid = fluid_temps[sft_index][-1]
        roots = optimize.root(wall_flux_root, T_fluid, args=(T_fluid))
        T_wall = roots.x[0]
        Q = radiative_loss(T_wall, T_fluid)
        fluid_temps[sft_index].append(T_fluid - Q/thermal_inertia)


ft0 = fluid_temps[0]
maxlen = len(ft0)
xft = np.linspace(0, maxlen-1, maxlen)



interp_ft0 = interpolate.interp1d(xft, ft0, "quadratic")

temps_in = np.linspace(350, 190, 81)
temps_out = np.linspace(350, 190, 81)

merits = np.zeros((len(temps_in), len(temps_out)))

for ti_index, temp_in in enumerate(temps_in):
    for to_index, temp_out in enumerate(temps_out):
        if temp_in > temp_out:
            in_index = optimize.root(lambda x : interp_ft0(x)-temp_in, 0).x[0]
            out_index = optimize.root(lambda x : interp_ft0(x)-temp_out, 0).x[0]
            merits[ti_index, to_index] = (temp_in - temp_out) / (out_index - in_index) 
        else:
            merits[ti_index, to_index] = np.nan
#ft_330 = optimize.root(lambda x : interp_ft0(x)-349, 0).x[0]
#print(ft_330)

print(np.nanmax(merits))

rad_merit = plt.contourf(temps_out, temps_in, merits, levels=np.linspace(0, 1, 21))
plt.colorbar(rad_merit)
plt.ylabel("Inlet temperature")
plt.xlabel("Outlet temperature")
plt.title("Radiator merit")



#plt.plot(xft, ft0)
plt.show()