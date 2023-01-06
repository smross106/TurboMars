import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

## Global variables

# Peak efficiency of each type 
eta_nom_axi = 0.89
eta_nom_rad = 0.83
Re_nom = 1.5e5

# Assumed all efficiencies scale as Re^scaling factor
# TODO - use Casey and Dietmann to improve https://www.euroturbo.eu/paper/ETC2013-052.pdf
Re_scaling_axi = 0.01
Re_scaling_rad = 0.037

# Non-dim parameters
psi_axi = 0.4
psi_rad = 0.4
Ns_axi = 2.0
Ns_rad = 0.6
Ds_axi = 2.0
Ds_rad = 4.5

# Speed of sound, assuming T=270K
a = np.sqrt(270 * 1.33 * 189)
mu = 14e-6

mass_flows = np.linspace(1000/(24*3600), 1000/(4*3600), 50)
densities = np.logspace(-3, 2, 50)

# Testing set
#mass_flows = np.linspace(1000/(24*3600), 1000/(4*3600), 5)
#densities = np.logspace(-3, 2, 6)

efficiencies_axi = np.zeros((len(mass_flows), len(densities)))
efficiencies_rad = np.zeros((len(mass_flows), len(densities)))

Res_axi = np.zeros((len(mass_flows), len(densities)))
Res_rad = np.zeros((len(mass_flows), len(densities)))

for mf_index, mass_flow in enumerate(mass_flows):
    for den_index, density in enumerate(densities):
        # Axial fan
        # Get delta_h using stage loading
        delta_h_axi = psi_axi * np.power(a * 0.8, 2)
        # Get speed using specific speed, and diameter using specific diameter
        omega_axi = Ns_axi / (np.power(mass_flow / density, 0.5) * np.power(delta_h_axi, -0.75))
        D_axi = Ds_axi * mass_flow / (np.power(delta_h_axi, 0.25) * density)
        Re_axi = density * omega_axi * np.power(D_axi,2) / mu
        eta_axi = eta_nom_axi * np.power(Re_axi/Re_nom, Re_scaling_axi)

        Res_axi[mf_index][den_index] = Re_axi
        efficiencies_axi[mf_index][den_index] = eta_axi
        
        # Radial machine
        # Get delta_h using stage loading
        delta_h_rad = psi_rad * np.power(a * 0.8, 2)
        # Get speed using specific speed, and diameter using specific diameter
        omega_rad = Ns_rad / (np.power(mass_flow / density, 0.5) * np.power(delta_h_rad, -0.75))
        D_rad = Ds_rad * mass_flow / (np.power(delta_h_rad, 0.25) * density)
        Re_rad = density * omega_rad * np.power(D_rad,2) / mu
        eta_rad = eta_nom_rad * np.power(Re_rad/Re_nom, Re_scaling_rad)

        Res_rad[mf_index][den_index] = Re_rad
        efficiencies_rad[mf_index][den_index] = eta_rad

        #print(int(mass_flow*1000), density, eta_axi, eta_rad)

fig, axs  = plt.subplots(2,2)

axi_re = axs[0,0].contourf(densities, mass_flows*1000, Res_axi, locator=ticker.LogLocator(subs=(0.5, 1.0)))#levels=np.logspace(2, 7, 11))
rad_re = axs[0,1].contourf(densities, mass_flows*1000, Res_rad, locator=ticker.LogLocator(subs=(0.5, 1.0)))#levels=np.logspace(2, 7, 11))
axs[0,0].contour(densities, mass_flows*1000, Res_axi, levels=[1.5e5], colors="k")
axs[0,1].contour(densities, mass_flows*1000, Res_rad, levels=[1.5e5], colors="k")

axi_cf = axs[1,0].contourf(densities, mass_flows*1000, efficiencies_axi, levels=np.linspace(0.65, 0.95, 13))
rad_cf = axs[1,1].contourf(densities, mass_flows*1000, efficiencies_rad, levels=np.linspace(0.65, 0.95, 13))


print(np.amin(Res_axi))
print(np.amin(Res_rad))
print(np.amax(Res_axi))
print(np.amax(Res_rad))

plt.colorbar(rad_re, ax=axs[0,:])
plt.colorbar(axi_cf, ax=axs[1,:])


for i in range(0,2):
    for j in range(0,2):
        axs[i,j].set_xscale("log")
        axs[i,j].set_xlabel("Exit density (kg/m3)")
        axs[i,j].set_ylabel("Mass flow (g/s)")

axs[0,0].set_title("Axial compressor Reynolds number")
axs[0,1].set_title("Radial compressor Reynolds number")
axs[1,0].set_title("Axial compressor efficiency")
axs[1,1].set_title("Radial compressor efficiency")

plt.show()