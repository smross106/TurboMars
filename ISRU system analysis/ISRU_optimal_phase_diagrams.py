import numpy as np
import matplotlib.pyplot as plt

option_names = []

options = {
    "CO2 24h, water 24h, reactor 24h":{
    "Weight":     18244,
    "Power":    1512861,
    "Cooling":   470403,
    "Energy":      9077,
    },
    "CO2 24h, water 24h, reactor day":{
    "Weight":     32399,
    "Power":    2240181,
    "Cooling":  1395210,
    "Energy":       349
    },
    "CO2 24h, water day, reactor 24h":{
    "Weight":     25306,
    "Power":    1533972,
    "Cooling":   470403,
    "Energy":      8823
    },
    "CO2 24h, water day, reactor day":{
    "Weight":     39461,
    "Power":    2261292,
    "Cooling":  1395210,
    "Energy":        96
    },
    "CO2 night, water 24h, reactor 24h":{
    "Weight":     18522,
    "Power":    1290023,
    "Cooling":   378458,
    "Energy":      7830
    },
    "CO2 night, water 24h, reactor day":{
    "Weight":     32676,
    "Power":    1906418,
    "Cooling":  1054436,
    "Energy":       433
    },
    "CO2 night, water day, reactor 24h":{
    "Weight":     25584,
    "Power":    1311134,
    "Cooling":   378458,
    "Energy":      7576
    },
    "CO2 night, water day, reactor day":{
    "Weight":     39738,
    "Power":    1927529,
    "Cooling":  1054436,
    "Energy":       180
    }
}

power_range = np.linspace(.008, .151*1.4, 150)
cooling_range = np.linspace(.05, .121*3, 100)
#power_range = np.logspace(-8, 2, 150)
#cooling_range = np.logspace(-8, 2, 100)
battery = 5

ESM = np.zeros((power_range.size, cooling_range.size))
optimal_system = np.zeros((power_range.size, cooling_range.size))
for power_i, power in enumerate(power_range):
    for cooling_i, cooling in enumerate(cooling_range):
        ESM_options = []
        for key in options.keys():
            ESM_options.append(options[key]["Weight"] * (4-key.count("24h")) + 
                               options[key]["Power"] * power +
                               options[key]["Cooling"] * cooling +
                               options[key]["Energy"] * power)
        ESM[power_i, cooling_i] = min(ESM_options)
        optimal_system[power_i, cooling_i] = ESM_options.index(min(ESM_options))

plt.contourf(cooling_range, power_range, ESM, cmap="viridis")
plt.colorbar()

plt.contour(cooling_range, power_range, optimal_system, colors="r", levels=[2])

plt.hlines([.012, .0225, .104, .149], cooling_range[0], cooling_range[-1], colors="k", linestyles="--")
plt.annotate("Rankine SP-100 580kWe nuclear", [cooling_range[-1], .012],horizontalalignment='right')
plt.annotate("ATK MegaFlex solar - orbital", [cooling_range[-1], .0225],horizontalalignment='right')
plt.annotate("BVAD solar, 40% efficiency", [cooling_range[-1], .104],horizontalalignment='right')
plt.annotate("BVAD solar, 28% efficiency", [cooling_range[-1], .149],horizontalalignment='right')


plt.vlines([.121, .149, .324], power_range[0], power_range[-1], colors="k", linestyles="--")
plt.annotate("BVAD lightweight radiators", [.121, power_range[-1]], rotation_mode = 'anchor',rotation=270)
plt.annotate("BVAD radiators", [.149, power_range[-1]], rotation_mode = 'anchor',rotation=270)
plt.annotate("ISS radiators", [.324, power_range[-1]], rotation_mode = 'anchor',rotation=270)

plt.xlabel("Mass-cost of cooling, kg/W")
plt.ylabel("Mass-cost of power generation, kg/W")
plt.title("Equivalent System Mass (kg-e) of entire fuel production stack \n Above red line, CO2 acqusition runs at night only \n Below red line, all systems run continuously")

plt.show()

