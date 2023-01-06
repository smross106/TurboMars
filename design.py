import component

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import time

def radial_smyth_optimisation(input_args, *args):
    (speed, work_coeff, flow_coeff, radius_hub_inlet) = input_args
    (gasflow_in, work_stage, efficiency) = args

    speed_rad = speed * np.pi / 30
    
    #SINGLE-EQUATION APPROACH
    #DON'T WORK
    func = np.power(speed_rad, 2) / np.sqrt(np.power(flow_coeff, 2) + work_coeff)
    func *= 0.5 * (radius_hub_inlet + np.sqrt(
        (np.power(2, 1/3) * np.power(gasflow_in.mass_flow / (gasflow_in.density() * np.pi * speed_rad), 2/3))
        + np.power(radius_hub_inlet, 2)
    ))
    func -= 16.16 * work_stage
    func = func*func
    

    """flow_func = 1 / (np.power(speed_rad, 2/3) * np.sqrt(work_coeff))
    flow_func -= np.power(gasflow_in.mass_flow * np.power(np.pi, 2) / (2 * gasflow_in.density()), 1/3) * 5 * 2.2 / np.sqrt(work_stage)

    R_mean = 0.5 * (radius_hub_inlet + np.sqrt(
        np.power(2, 1/3) * np.power(gasflow_in.mass_flow / (gasflow_in.density() * np.pi * speed_rad), 2/3) + np.power(radius_hub_inlet, 2)
    ))

    work_func = np.sqrt((np.power(flow_coeff, 2) / work_coeff) + 1) / np.power(speed_rad, 8/3)
    work_func *= 4 / np.power(R_mean, 2)
    work_func -= np.power((gasflow_in.mass_flow * np.power(np.pi, 2)) / (2 * gasflow_in.density()), 1/3) * (5.5 / 0.85) * 0.25 * np.power(work_stage, -0.5)

    func = np.sqrt(flow_func**2 + work_func**2)"""
    
    return(func)

def make_axial(input_args, *args, output=False):
    (speed, work_coeff, flow_coeff) = input_args
    (gasflow_in, work_stage, efficiency) = args

    axial = component.AxialStage(gasflow_in, work_stage, speed, work_coeff, flow_coeff, efficiency)

    # Assume everything is 7075 aluminium
    yield_stress = 381e6 # Safety factor of 1.5
    density = 2810

    axial_weight = axial.weight_estimate()

    input_power = gasflow_in.mass_flow * work_stage / axial.estimate_efficiency()
    ESM_power_weight = input_power * .121

    constraint_weight = 0
    if axial.bladeheight_in < 3e-3:
        constraint_weight += abs(axial.bladeheight_in - 3e-3) * 5e8
    
    if speed < 1:
        constraint_weight += 1e15
    elif speed > 100000:
        constraint_weight += abs(flow_coeff - 100000)**2 * 1e4

    if axial.R_hub_in < 0:
        constraint_weight += 1e10
    elif axial.R_hub_in < 10e-3:
        constraint_weight += abs(axial.R_hub_in - 10e-3)**2 * 1e7
    elif axial.R_hub_in > 1:
        constraint_weight += abs(axial.R_hub_in - 1)**2 * 1e5

    if work_coeff > 0.85:
        constraint_weight += abs(work_coeff - 0.85) * 1e5
    elif work_coeff < 1e-3:
        constraint_weight += abs(1 - work_coeff) * 1e5

    if flow_coeff > 2:
        constraint_weight += abs(flow_coeff - 2)**2 * 1e5
    elif flow_coeff < 0.25:
        constraint_weight += abs(flow_coeff - 0.25)**2 * 1e5

    if output:
        return(axial_weight+ESM_power_weight)
    else:
        return(axial_weight+ESM_power_weight+constraint_weight)

def make_radial(input_args, *args, output=False):
    (speed, work_coeff, flow_coeff, radius_hub_inlet, diffusion_ratio) = input_args
    (gasflow_in, work_stage, efficiency) = args

    radial = component.RadialStage(gasflow_in, work_stage, speed, work_coeff, flow_coeff, efficiency, radius_hub_inlet, diffusion_ratio=diffusion_ratio)

    # Assume everything is 7075 aluminium
    yield_stress = 381e6 # Safety factor of 1.5
    density = 2810

    radial_weight = radial.weight_estimate(yield_stress, density)

    input_power = gasflow_in.mass_flow * work_stage / radial.estimate_efficiency()
    ESM_power_weight = input_power * .121

    constraint_weight = 0
    # Add weight factors for constraints
    
    if diffusion_ratio < 1:
        constraint_weight += 1e15
    elif diffusion_ratio > 15:
        constraint_weight += 1e15
    
    if radial.bladeheight_in < 3e-3:
        constraint_weight += abs(radial.bladeheight_in - 3e-3)**2 * 5e8
    
    if radial.bladeheight_imp_exit < 3e-3:
        constraint_weight += abs(radial.bladeheight_imp_exit - 3e-3)**2 * 5e8
    
    if radial.R_mean_imp_exit < radial.R_tip_inlet:
        constraint_weight += 1e15

    if work_coeff < 0.01:
        constraint_weight += 1e15
    elif work_coeff > 0.85:
        constraint_weight += abs(work_coeff - 0.85)**2 * 1e5
    elif work_coeff < 0.1:
        constraint_weight += abs(work_coeff - 0.1)**2 * 1e5

    if flow_coeff < 0.001:
        constraint_weight += 1e20
    elif flow_coeff > 2:
        constraint_weight += abs(flow_coeff - 2)**2 * 1e5
    elif flow_coeff < 0.01:
        constraint_weight += abs(flow_coeff - 0.01)**2 * 1e5
    
    if speed < 1:
        constraint_weight += 1e15
    elif speed > 100000:
        constraint_weight += abs(flow_coeff - 100000)**2 * 1e4
    
    if radius_hub_inlet < 0:
        constraint_weight += 1e15
    elif radius_hub_inlet < 10e-3:
        constraint_weight += abs(radius_hub_inlet - 10e-3)**2 * 1e7
    elif radius_hub_inlet > 1:
        constraint_weight += abs(radius_hub_inlet - 1)**2 * 1e5

    if output:
        return(radial_weight+ESM_power_weight)
    else:
        return(radial_weight + ESM_power_weight + constraint_weight)

def make_heat_exchanger(input_args, *args, output=False):
    (bulk_area, coolant_velocity, pitch_diameter_ratio) = input_args
    (cooling_power, gasflow, inlet_area, target_PD) = args

    if bulk_area < 1e-7:
        return(1e50)

    HX = component.HeatExchangerAdvanced(cooling_power, gasflow, inlet_area, 
        bulk_area, pitch_diameter_ratio=pitch_diameter_ratio, coolant_velocity=coolant_velocity)

    HX.gas_pressure_drop()
    HX.coolant_pumping_power()
    HX.weight_estimate()

    hx_weight = HX.weight

    ESM_power_weight = HX.pumping_power * .121

    ESM_volume_weight = (HX.duct_length + 2*HX.diffuser_length) * bulk_area * 100

    constraint_weight = 0

    constraint_weight += np.power(abs(target_PD - HX.pressure_drop), 2) * 1e3

    if pitch_diameter_ratio < 1.1:
        constraint_weight += 1e20
    elif pitch_diameter_ratio > 50:
        constraint_weight += abs(pitch_diameter_ratio - 50) **2 * 1e6
    
    if (pitch_diameter_ratio * HX.tube_diameter) > HX.duct_size:
        constraint_weight += 1e20
    
    if bulk_area/inlet_area < 1:
        constraint_weight += 1e15
    elif bulk_area/inlet_area > 15:
        constraint_weight += (bulk_area/inlet_area) * 1e6
    
    if coolant_velocity < 0.1:
        constraint_weight += 1e15
    
    if (HX.duct_length + 2*HX.diffuser_length) > 3:
        constraint_weight += np.power((HX.duct_length + 2*HX.diffuser_length - 3), 2) * 1e10

    if output:
        return(hx_weight + ESM_power_weight)
    else:
        return(hx_weight + ESM_power_weight + constraint_weight)


def optimise_axial(flow_in, deltah, efficiency):
    specific_speed = 2.0
    specific_diameter = 2.0

    speed_guess = specific_speed / (np.sqrt(flow_in.mass_flow / flow_in.density_deltah(deltah)) * np.power(deltah, -0.75)) 
    work_coeff_guess = 1 / np.power(specific_speed * specific_diameter, 2)
    flow_coeff_guess = 1 / (specific_speed * np.power(specific_diameter, 3))

    #crude_weight = make_axial([speed_guess, work_coeff_guess, flow_coeff_guess], flow_in, deltah, efficiency)
    x0 = np.array([speed_guess, work_coeff_guess, flow_coeff_guess])
    x0_upper = np.array([99.9e3, 0.6, 2])
    N = len(x0)
    zdelt = 0.00025
    sim = np.empty((N + 1, N), dtype=x0.dtype)
    sim[0] = x0
    for k in range(N):
        y = 0.8*np.array(x0, copy=True)
        if y[k] != 0:
            y[k] = x0_upper[k]
        else:
            y[k] = zdelt
        sim[k + 1] = y

    properties = optimize.minimize(make_axial, 
        x0, 
        args=(flow_in, deltah, efficiency),
        #bounds=[(0, 100e3), (0, 0.5), (0, 2)], 
        #constraints=cons, method="SLSQP"
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":1500}
        method="Nelder-Mead", options={"maxiter":10000, "initial_simplex":sim}
        )


    (speed, work_coeff, flow_coeff) = properties.x

    axial = component.AxialStage(flow_in, deltah, speed, work_coeff, flow_coeff, efficiency)
    axial.weight_estimate()
    axial.estimate_efficiency()
    axial.ESM_estimate()
    #opt_weight = make_axial([speed, work_coeff, flow_coeff], flow_in, deltah, efficiency)

    #return([crude_weight, opt_weight, axial])
    return(axial)

def optimise_radial(flow_in, deltah, efficiency):
    specific_speed = 0.5
    specific_diameter = 5.0

    speed_guess = min((specific_speed / ((np.sqrt(flow_in.mass_flow / flow_in.density_deltah(deltah)) * np.power(deltah, -0.75)) )) * 30 / np.pi, 99e3)
    work_coeff_guess = 1 / np.power(specific_speed * specific_diameter, 2)
    flow_coeff_guess = 1 / (specific_speed * np.power(specific_diameter, 3))
    inlet_hub_radius_guess = max(400e-3, 0.5*np.sqrt(flow_in.mass_flow / (np.pi * flow_in.density() * 100 * flow_coeff_guess)))
    diffusion_ratio_guess = 3

    #crude_weight = make_radial([speed_guess, work_coeff_guess, flow_coeff_guess, inlet_hub_radius_guess], flow_in, deltah, efficiency)
    x0 = np.array([speed_guess, work_coeff_guess, flow_coeff_guess, inlet_hub_radius_guess, diffusion_ratio_guess])
    x0_upper = np.array([99.9e3, 0.85, 2, 1, 15])
    N = len(x0)
    zdelt = 0.00025
    sim = np.empty((N + 1, N), dtype=x0.dtype)
    sim[0] = x0
    for k in range(N):
        y = 0.8*np.array(x0, copy=True)
        if y[k] != 0:
            y[k] = x0_upper[k]
        else:
            y[k] = zdelt
        sim[k + 1] = y
    #print(sim)

    properties = optimize.minimize(make_radial, 
        x0, 
        args=(flow_in, deltah, efficiency),
        bounds=[(1, 100e3), (0.1, 0.85), (0.2, 2), (10e-3, 1), (1, 15)],
        method="Nelder-Mead", options={"maxiter": 10000, "disp":False, "adaptive":False, "return_all":False, "xatol":1e-4, "initial_simplex":sim}
        #method='SLSQP'
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":2000}
        )
    """properties = optimize.basinhopping(make_radial, 
        x0=[speed_guess, work_coeff_guess, flow_coeff_guess, inlet_hub_radius_guess, diffusion_ratio_guess],
        niter=n_iter_try,
        minimizer_kwargs={
                "method":"Nelder-Mead",
                "args":(flow_in, deltah, efficiency)},
        )"""
    
    #print(properties)
    
    (speed, work_coeff, flow_coeff, inlet_hub_radius, diffusion_ratio) = properties.x

    radial = component.RadialStage(flow_in, deltah, speed, work_coeff, flow_coeff, efficiency, inlet_hub_radius, diffusion_ratio)
    radial.weight_estimate(381e6, 2810)
    radial.estimate_efficiency()
    radial.ESM_estimate()

    #opt_weight = make_radial([speed, work_coeff, flow_coeff, inlet_hub_radius], flow_in, deltah, efficiency, output=True)

    return(radial)


def optimise_hx(flow_in, cooling_deltah, inlet_area, target_PD):
    
    bulk_area_guess = inlet_area * 2
    coolant_velocity_guess = 1
    pitch_diameter_ratio_guess = 10

    cooling_power = flow_in.mass_flow * cooling_deltah


    x0 = np.array([bulk_area_guess, coolant_velocity_guess, pitch_diameter_ratio_guess])
    x0_upper = np.array([bulk_area_guess*15, 2, 50])
    N = len(x0)
    zdelt = 0.00025
    sim = np.empty((N + 1, N), dtype=x0.dtype)
    sim[0] = x0
    for k in range(N):
        y = 0.8*np.array(x0, copy=True)
        if y[k] != 0:
            y[k] = x0_upper[k]
        else:
            y[k] = zdelt
        sim[k + 1] = y
    #print(sim)

    properties = optimize.minimize(make_heat_exchanger, 
        x0, 
        args=(cooling_power, flow_in, inlet_area, target_PD) ,
        bounds=[(1e-7, 10), (0.1, 3), (1.1, 50)],
        method="Nelder-Mead", options={"maxiter": 10000, "disp":False, "adaptive":False, "return_all":False, "xatol":1e-4, "initial_simplex":sim}
        #method='SLSQP'
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":2000}
        )
    """properties = optimize.basinhopping(make_radial, 
        x0=[speed_guess, work_coeff_guess, flow_coeff_guess, inlet_hub_radius_guess, diffusion_ratio_guess],
        niter=n_iter_try,
        minimizer_kwargs={
                "method":"Nelder-Mead",
                "args":(flow_in, deltah, efficiency)},
        )"""
    
    #print(properties)
    
    (bulk_area, coolant_velocity, pitch_diameter_ratio) = properties.x

    HX = component.HeatExchangerAdvanced(cooling_power, flow_in, inlet_area, 
        bulk_area, pitch_diameter_ratio=pitch_diameter_ratio, coolant_velocity=coolant_velocity)
    HX.gas_pressure_drop()
    HX.coolant_pumping_power()
    HX.weight_estimate()
    HX.estimate_ESM()

    #opt_weight = make_radial([speed, work_coeff, flow_coeff, inlet_hub_radius], flow_in, deltah, efficiency, output=True)

    return(HX)

#DESIGN VALIDATION
"""
lp_flow = component.GasFlow(0.05, 190, 900)
mp_flow = component.GasFlow(0.05, 250, 1e5)

PR_axial = 1.1
PR_radial = 2.361


print("PR: ", PR_radial)"""





"""lp_radial = optimise_radial(lp_flow, lp_flow.delta_h_PR(PR_radial), 0.7)
print(lp_radial)
mp_radial = optimise_radial(mp_flow, mp_flow.delta_h_PR(PR_radial), 0.7)
print(mp_radial)
"""



# Axial
#delta_hs = np.linspace(0.5e3, 12e3, 24)

# Radial
"""delta_hs = np.linspace(3e3, 60e3, 20)
#delta_hs = [10e3]
print(delta_hs)
crudes_ax = []
optis_ax = []

optis_rad_ESM = []
optis_rad_physical = []
optis_rad_energy = []
rads = []

start = time.time()

for delta_h in delta_hs:
    radial_lp = optimise_radial(lp_flow, delta_h, 0.7)
    optis_rad_ESM.append(radial_lp.ESM)
    optis_rad_physical.append(radial_lp.weight)
    optis_rad_energy.append(radial_lp.ESM - radial_lp.weight)
    rads.append(radial_lp)"""
    #optis_rad_LP.append(weights_lp[1])
    #optis_rad_physical.append(weights_lp[2].weight_estimate(381e6, 2810))
    #optis_rad_energy.append(weights_lp[1] - weights_lp[2].weight_estimate(381e6, 2810))
    #rads.append(weights_lp[2])
    #optis_rad_energy.append(.121 * delta_h * lp_flow.mass_flow / weights_lp[2].estimate_efficiency())

    #print(delta_h, weights_lp[2].bladeheight_in, weights_lp[2].estimate_efficiency(verbose=False))

    #print(delta_h, weights[2].estimate_efficiency())

"""print((time.time() - start)/len(rads))

print("Design of radial at 30kJ/kg, PR=2.08")
print(rads[-1])
print(rads[-1].flow_coeff)
print(rads[-1].work_coeff)

#print("Min weight for optimised ",str(min(optis_rad)))

fig, ax1 = plt.subplots()

#ax1.scatter(delta_hs, optis_ax, label="axial")
#ax1.scatter(delta_hs, crudes_ax, label="crude guess")
ax1.plot(delta_hs, optis_rad_ESM, label="ESM",c="b")
ax1.plot(delta_hs, optis_rad_physical, label="Physical weight",c="r")
ax1.plot(delta_hs, optis_rad_energy, label="Energy weight", c="g")

#ax1.scatter(delta_hs, optis_rad, label="1 bar inlet")
#ax1.scatter(delta_hs, crudes_rad, label="crude guess")

ax1.set_xlabel('delta h')
ax1.set_ylabel('Weight')
#ax1.set_ylim(0, 10000)
fig.legend()
plt.show()
#print(make_axial((5000, 0.4, 1.0), bc_flow, 50e3, 0.85))

"""