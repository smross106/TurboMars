import numpy as np
import matplotlib.pyplot as plt

from component import HeatExchangerAdvanced, GasFlow


def velocities_sweep(flow, cooling_power, inlet_area, gasflow_range, coolantflow_range, wall=False):
    volumetric_flow = flow.mass_flow / flow.density()

    max_area = volumetric_flow / gasflow_range[1] 
    min_area = volumetric_flow / gasflow_range[0]

    areas = np.linspace(min_area, max_area, 100)
    gas_speeds = [volumetric_flow / area for area in areas]
    coolant_speeds = np.linspace(coolantflow_range[0], coolantflow_range[1], 100) 

    total_weights = np.zeros((len(areas), len(coolant_speeds)))
    dry_weights = np.zeros((len(areas), len(coolant_speeds)))
    pressure_drops = np.zeros((len(areas), len(coolant_speeds)))
    pressure_drops_percents = np.zeros((len(areas), len(coolant_speeds)))
    lengths = np.zeros((len(areas), len(coolant_speeds)))
    pumping_powers = np.zeros((len(areas), len(coolant_speeds)))
    n_tubes = np.zeros((len(areas), len(coolant_speeds)))
    n_rows = np.zeros((len(areas), len(coolant_speeds)))
    array_areas = np.zeros((len(areas), len(coolant_speeds)))
    array_gas_speeds = np.zeros((len(areas), len(coolant_speeds)))
    array_coolant_speeds = np.zeros((len(areas), len(coolant_speeds)))

    for area_index, area in enumerate(areas):
        for coolant_index, coolant in enumerate(coolant_speeds):
            if wall:
                hx = HeatExchangerAdvanced(cooling_power, flow, inlet_area, area, coolant_velocity=coolant, 
                glycol=False, pitch_diameter_ratio=5, tube_diameter=2e-3, wall_thickness=0.1e-3, geometry="wall")
            else:
                hx = HeatExchangerAdvanced(cooling_power, flow, inlet_area, area, coolant_velocity=coolant, 
                glycol=False, pitch_diameter_ratio=2, tube_diameter=2e-3, wall_thickness=0.1e-3, geometry="row")
            hx.weight_estimate()
            hx.coolant_pumping_power()
            hx.gas_pressure_drop()

            total_weights[area_index, coolant_index] = hx.weight + 0.121*hx.pumping_power
            dry_weights[area_index, coolant_index] = hx.weight
            pressure_drops[area_index, coolant_index] = hx.pressure_drop
            pressure_drops_percents[area_index, coolant_index] = 100 * hx.pressure_drop / flow.pressure
            lengths[area_index, coolant_index] = hx.duct_length
            pumping_powers[area_index, coolant_index] = hx.pumping_power
            n_tubes[area_index, coolant_index] = hx.n_tubes
            n_rows[area_index, coolant_index] = hx.n_rows

            array_areas[area_index, coolant_index] = area
            array_gas_speeds[area_index, coolant_index] = volumetric_flow / area 
            array_coolant_speeds[area_index, coolant_index] = coolant
    
    return({
        "Area": areas,
        "Gas speed": gas_speeds,
        "Coolant speed": coolant_speeds,
        "Total weight": total_weights,
        "Dry weight": dry_weights,
        "Pressure drop": pressure_drops,
        "Pressure drop percent": pressure_drops_percents,
        "Length": lengths,
        "Pumping power": pumping_powers,
        "Number of tubes": n_tubes,
        "Number of rows": n_rows,
        "Array area": array_areas,
        "Array gas speed": array_gas_speeds,
        "Array coolant speed": array_coolant_speeds
    })

def contourplot(x_range, y_range, data, unit="", log=[False, False]):
    datamax = np.amax(data)
    datamin = np.amin(data)
    datamax_log = np.ceil(2*np.log10(datamax))/2
    datamin_log = np.floor(2*np.log10(datamin))/2

    fig, ax = plt.subplots()
    wt = ax.contourf(x_range, y_range, data, levels=np.logspace(datamin_log, datamax_log, 15))
    plt.colorbar(wt, ax=ax)

    if log[0]:
        ax.set_xscale("log")
    if log[1]:
        ax.set_yscale("log")
    
    ax.text(.5,.95,('Minimum {:.1f}'+unit+', maximum {:.1f}'+unit).format(datamin, datamax),fontsize=10,ha='center', transform=ax.transAxes)


    return(fig, ax)



bc_flow = GasFlow(0.05, 300, 1900)
exit_flow = GasFlow(0.05, 300, 100e3)

lp_tubes = velocities_sweep(bc_flow, 2000, 8.8e-3, [1, 150], [0.5, 2], wall=False)
hp_tubes = velocities_sweep(exit_flow, 2000, 1.7e-4, [1, 150], [0.5, 2], wall=False)

lp_wall = velocities_sweep(bc_flow, 2000, 8.8e-3, [1, 150], [0.5, 2], wall=True)
hp_wall = velocities_sweep(exit_flow, 2000, 1.7e-4, [1, 150], [0.5, 2], wall=True)


"""fig, ax = contourplot(lp_tubes["Coolant speed"], lp_tubes["Gas speed"], lp_tubes["Length"], unit="m", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("HX length, 1.9kPa")

fig, ax = contourplot(lp_tubes["Coolant speed"], lp_tubes["Gas speed"], lp_tubes["Number of rows"], unit="", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Number of rows, 1.9kPa")

fig, ax = contourplot(lp_tubes["Coolant speed"], lp_tubes["Gas speed"], lp_tubes["Total weight"], unit="kg", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Total weight, 1.9kPa")

fig, ax = contourplot(lp_tubes["Coolant speed"], lp_tubes["Gas speed"], lp_tubes["Pumping power"], unit="W", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Coolant pumping power, 1.9kPa")"""


"""

fig, ax = contourplot(hp_tubes["Coolant speed"], hp_tubes["Gas speed"], hp_tubes["Total weight"], unit="kg", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Total weight, 100kPa")

fig, ax = contourplot(hp_tubes["Coolant speed"], hp_tubes["Gas speed"], hp_tubes["Pumping power"], unit="W", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Coolant pumping power, 100kPa")

fig, ax = contourplot(hp_tubes["Coolant speed"], hp_tubes["Gas speed"], hp_tubes["Number of tubes"], unit="", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Number of tubes, 100kPa")"""

fig, ax = contourplot(lp_tubes["Coolant speed"], lp_tubes["Gas speed"], lp_tubes["Total weight"], unit="kg", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Total weight, 1.9kPa, tube rows")

fig, ax = contourplot(lp_wall["Coolant speed"], lp_wall["Gas speed"], lp_wall["Total weight"], unit="kg", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Total weight, 1.9kPa, sintered tubes")

fig, ax = contourplot(lp_tubes["Coolant speed"], lp_tubes["Gas speed"], lp_tubes["Pumping power"], unit="W", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Pump power, 1.9kPa, tube rows")

fig, ax = contourplot(lp_wall["Coolant speed"], lp_wall["Gas speed"], lp_wall["Pumping power"], unit="W", log=[False, True])
ax.set_xlabel("Coolant speed, m/s")
ax.set_ylabel("Gas speed, m/s")
ax.set_title("Pump power, 1.9kPa, sintered tubes")

print("hp tubes", np.amin(hp_tubes["Pressure drop"]))
print("hp tubes", np.amax(hp_tubes["Pressure drop"]))
print("lp tubes", np.amin(lp_tubes["Pressure drop"]))
print("lp tubes", np.amax(lp_tubes["Pressure drop"]))

print("hp wall", np.amin(hp_wall["Pressure drop"]))
print("hp wall", np.amax(hp_wall["Pressure drop"]))
print("lp wall", np.amin(lp_wall["Pressure drop"]))
print("lp wall", np.amax(lp_wall["Pressure drop"]))

"""fig, ax = plt.subplots()
ax.plot([item for sublist in hp_tubes["Array gas speed"] for item in sublist], 
    [item for sublist in hp_tubes["Pressure drop percent"] for item in sublist], label="100kPa tubes")
ax.plot([item for sublist in lp_tubes["Array gas speed"] for item in sublist], 
    [item for sublist in lp_tubes["Pressure drop percent"] for item in sublist], label="1.9kPa tubes")

ax.plot([item for sublist in hp_wall["Array gas speed"] for item in sublist], 
    [item for sublist in hp_wall["Pressure drop percent"] for item in sublist], label="100kPa sintered")
ax.plot([item for sublist in lp_wall["Array gas speed"] for item in sublist], 
    [item for sublist in lp_wall["Pressure drop percent"] for item in sublist], label="1.9kPa sintered")

plt.vlines(7, 0, 1000, colors="r", linestyles="--", label="Diffuser maximum")
plt.vlines(100, 0, 1000, colors="g", linestyles="--", label="Compressor outlet")
plt.hlines(5, 0, 1000, colors="b", linestyles="--", label="Previous assumption")
plt.legend()
ax.set_xlabel("Flow velocity, m/s")
ax.set_ylabel("Pressure drop, %")
ax.set_title("Pressure drop in heat exchanger")
ax.set_ylim(0, 25)
ax.set_xscale("log")
ax.set_xlim(1, 150)"""

plt.show()

"""bulk_areas = np.logspace(-3, 1, 45)
velocities = np.linspace(0.5, 2, 34)
total_weights_lp = []
total_weights_hp = []
total_lengths_vlt = []

pressure_drops_lp = []
pressure_drops_hp = []
speeds_lp = []
speeds_hp = []

for bulk_area in bulk_areas:
    hp_weights = []
    lp_weights = []
    vlt_lengths = []
    for velocity in velocities:
        hx1 = HeatExchangerAdvanced(2000, exit_flow, 1e-3, bulk_area, coolant_velocity=velocity, glycol=False, wall_thickness=0.1e-3)
        hx1.weight_estimate()
        hx1.coolant_pumping_power()
        hp_weights.append(hx1.weight + 0.121 * hx1.pumping_power)

        hx1 = HeatExchangerAdvanced(2000, bc_flow, 1e-3, bulk_area, coolant_velocity=velocity, glycol=False, wall_thickness=0.1e-3)
        hx1.weight_estimate()
        hx1.coolant_pumping_power()
        lp_weights.append(hx1.weight + 0.121 * hx1.pumping_power)
        vlt_lengths.append(hx1.duct_length)
        #print(hx1.weight + 0.121 * hx1.pumping_power, hx1.duct_length)
    
    total_lengths_vlt.append(vlt_lengths)
    total_weights_lp.append(lp_weights)
    total_weights_hp.append(hp_weights)
    hx1.gas_pressure_drop()
    pressure_drops_lp.append(100 * hx1.pressure_drop / bc_flow.pressure)
    speeds_lp.append(bc_flow.mass_flow / (bc_flow.density() * bulk_area))

    hx1 = HeatExchangerAdvanced(2000, exit_flow, 5e-3, bulk_area, coolant_velocity=velocity, glycol=False)
    hx1.gas_pressure_drop()
    hx1.weight_estimate()
    pressure_drops_hp.append(100 * hx1.pressure_drop / exit_flow.pressure)
    speeds_hp.append(exit_flow.mass_flow / (exit_flow.density() * bulk_area))

print("max weight LP", np.amax(total_weights_lp))
print("min weight LP", np.amin(total_weights_lp))

print("max weight HP", np.amax(total_weights_hp))
print("min weight HP", np.amin(total_weights_hp))

print("max length", np.amax(total_lengths_vlt))
print("min length", np.amin(total_lengths_vlt))

hx1 = HeatExchangerAdvanced(2000, bc_flow, 1e-3, 0.5, coolant_velocity=0.5, glycol=False, wall_thickness=0.1e-3)
hx1.weight_estimate()
hx1.coolant_pumping_power()
print(hx1.duct_length, hx1.weight, hx1.pumping_power, hx1.n_tubes)

fig, ax = plt.subplots()
wt = ax.contourf(velocities, speeds_lp, total_weights_lp, levels=np.logspace(1, 3, 20))
#ax.contour(velocities, bulk_areas, total_weights_lp, levels=[100], colors="r")
ax.set_yscale("log")
ax.set_title("Heat exchanger weight at 1.9kPa inc. pump power")
ax.set_xlabel("Coolant flow speed, m/s")
ax.set_ylabel("Gas flow speed, m/s")
plt.colorbar(wt, ax=ax)

fig, ax = plt.subplots()
wt = ax.contourf(velocities, speeds_hp, total_weights_hp, levels=np.logspace(1, 3, 20))
#ax.contour(velocities, bulk_areas, total_weights_hp, levels=[100], colors="r")
ax.set_yscale("log")
ax.set_title("Heat exchanger weight at 100kPa inc. pump power")
ax.set_xlabel("Coolant flow speed, m/s")
ax.set_ylabel("Gas flow speed, m/s")
plt.colorbar(wt, ax=ax)


fig, ax = plt.subplots()
wt = ax.contourf(velocities, speeds_lp, total_lengths_vlt, levels=np.logspace(0, 2.5, 15))
ax.contour(velocities, speeds_lp, total_lengths_vlt, levels=[5], colors="r")
ax.set_yscale("log")
ax.set_title("Heat exchanger duct length")
ax.set_xlabel("Coolant flow speed, m/s")
ax.set_ylabel("Gas flow speed, m/s")
plt.colorbar(wt, ax=ax)

fig, ax = plt.subplots()
ax.scatter(speeds_lp, pressure_drops_lp, label="1.9kPa")
ax.scatter(speeds_hp, pressure_drops_hp, label="100kPa")
plt.hlines(5, 1e-3, 1e1, label="Previous assumption")
plt.vlines(0.5, 0, 20, colors="k", linestyles="--", label="Approx optimum")
ax.set_xscale("log")
ax.set_ylim(0, 20)
ax.set_title("Pressure drop with varying duct area")
ax.set_xlabel("Gas flow speed, m/s")
ax.set_ylabel("Pressure drop, %")
ax.legend()

plt.show()"""