import numpy as np
import matplotlib.pyplot as plt
from math import log

import component, Ts_cycle, design, co2_load


# STEP 1: design a cycle
mass_flow = 0.05
pressure_ratios, efficiencies, IC_temps, IC_pressure_drops = Ts_cycle.generate_radiator_cycle(5, 190, 900, 600e3, 250, 5e-4)
print("Cycle designed")

# STEP 2: get carbon dioxide properties
fig, ax, Trange, srange, hgrid, pgrid, dgrid = Ts_cycle.main()
print("Properties loaded")

# Step 3: get the properties of each component in the cycle
machine_params = Ts_cycle.plot_ideal_process(ax, [190, 900], pressure_ratios, efficiencies, IC_temps, IC_pressure_drops, Trange, srange, hgrid, pgrid, dgrid, label="1", return_machine_params=True)
print("Machine parameters created")

# Step 4: decide if each machine should be an axial or radial

machines = []
max_axial_PR = 1.2
max_radial_PR = 2.4

for machine in machine_params[-2:]:
    if machine["type"] == "c":
        pressure_ratio = machine["pressure ratio"]
        flow = component.GasFlow(mass_flow, machine["T in"], machine["p in"])
        if pressure_ratio < max_axial_PR:
            # Single stage axial, or radial
            axial = design.optimise_axial(flow, machine["delta h"], 0.7)[2]
            print(axial)

            radial = design.optimise_radial(flow, machine["delta h"], 0.7)[2]
            print(radial)

            print("\n")
        
        elif pressure_ratio < max_radial_PR:
            # Single stage radial, multi-stage axial
            n_axials = int(np.ceil(np.log(pressure_ratio) / np.log(max_axial_PR)))
            axial_delta_h = machine["delta h"] / n_axials
            axials = []
            for n in range(n_axials):
                axial = design.optimise_axial(flow, machine["delta h"], 0.7)[2]
                print(axial)
                axials.append(axial)
                flow = axial.gasflow_out()

            radial = design.optimise_radial(flow, machine["delta h"], 0.7)[2]
            print(radial)

        else:
            
            n_axials = int(np.ceil(np.log(pressure_ratio) / np.log(max_axial_PR)))
            axial_delta_h = machine["delta h"] / n_axials
            axials = []
            for n in range(n_axials):
                axial = design.optimise_axial(flow,axial_delta_h, 0.7)[2]
                print(axial)
                axials.append(axial)
                flow = axial.gasflow_out()
            print("\n")
            
            flow = component.GasFlow(mass_flow, machine["T in"], machine["p in"])
            n_radials = int(np.ceil(np.log(pressure_ratio) / np.log(max_radial_PR)))
            radial_delta_h = machine["delta h"] / n_radials
            radials = []
            for n in range(n_radials):
                radial = design.optimise_radial(flow, radial_delta_h, 0.7)[2]
                print(radial)
                radials.append(radial)
                flow = radial.gasflow_out()

            radial = design.optimise_radial(flow, machine["delta h"], 0.7)[2]
            print(radial)
            print("# # \n")
            print(pressure_ratio, n_axials, n_radials)

            
            
            