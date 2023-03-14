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
    (gasflow_in, pressure_ratio_stage, efficiency, ESM_vector) = args

    axial = component.AxialStage(gasflow_in, pressure_ratio_stage, speed, work_coeff, flow_coeff, efficiency)

    # Assume everything is 7075 aluminium
    yield_stress = 381e6 # Safety factor of 1.5
    density = 2810

    axial_weight = axial.weight_estimate()

    input_power = gasflow_in.mass_flow * axial.work_stage / axial.estimate_efficiency()

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
        return((axial_weight*ESM_vector[0])+(input_power*ESM_vector[1]))
    else:
        return((axial_weight*ESM_vector[0])+(input_power*ESM_vector[1])+constraint_weight)

def make_radial(input_args, *args, output=False):
    (speed, work_coeff, flow_coeff, radius_hub_inlet, diffusion_ratio) = input_args
    (gasflow_in, pressure_ratio_stage, efficiency_guess, ESM_vector) = args

    radial = component.RadialStage(gasflow_in, pressure_ratio_stage, speed, work_coeff, flow_coeff, efficiency_guess, radius_hub_inlet, diffusion_ratio=diffusion_ratio)

    # Assume everything is 7075 aluminium
    yield_stress = 381e6 # Safety factor of 1.5
    density = 2810

    radial_weight = radial.weight_estimate(yield_stress, density)

    input_power = gasflow_in.mass_flow * radial.work_stage / radial.estimate_efficiency()

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
        constraint_weight += np.power(radial.R_tip_inlet / radial.R_mean_imp_exit, 2) * 1e40

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

    if np.isnan([radial_weight]):
        radial_weight = 1e50
    if np.isnan([input_power]):
        input_power = 1e50

    if output:
        return((radial_weight*ESM_vector[0])+(input_power*ESM_vector[1]))
    else:
        return((radial_weight*ESM_vector[0])+(input_power*ESM_vector[1]) + constraint_weight)

def make_heat_exchanger_anypd(input_args, *args, output=False):
    (bulk_area_ratio, coolant_velocity, gap_diameter_ratio) = input_args
    (cooling_power, gasflow, inlet_area, ESM_vector) = args

    coolant_velocity = abs(coolant_velocity)

    limit_weight = 0
    if bulk_area_ratio < 1:
        limit_weight += np.power(bulk_area_ratio-1, 2)*1e50
    elif bulk_area_ratio > 15:
        limit_weight += np.power(bulk_area_ratio-50, 2)*1e50
    
    if coolant_velocity < 0.1:
        limit_weight += np.power(coolant_velocity-0.1, 2)*1e50
    elif coolant_velocity > 1:
        limit_weight += np.power(coolant_velocity-1, 2)*1e50
    
    if gap_diameter_ratio < 0.1:
        limit_weight += np.power(gap_diameter_ratio-0.1, 2)*1e50
    elif gap_diameter_ratio >50:
        limit_weight += np.power(gap_diameter_ratio-50, 2)*1e50
    
    if limit_weight != 0:
        return(limit_weight)

    HX = component.HeatExchangerAdvanced(cooling_power, gasflow, inlet_area, 
        bulk_area_ratio, pitch_diameter_ratio=gap_diameter_ratio+1, coolant_velocity=coolant_velocity, coolant="ammonia")

    HX.gas_pressure_drop()
    HX.coolant_pumping_power()
    HX.weight_estimate()

    relative_pd = 100*HX.pressure_drop/gasflow.pressure

    hx_weight = HX.weight

    ESM_power_weight = HX.pumping_power * .149

    ESM_volume_weight = 0*(HX.duct_length + 2*HX.diffuser_length) * bulk_area_ratio * inlet_area * 100

    constraint_weight = 0

    constraint_weight += relative_pd * 80 * 5
    
    if ((gap_diameter_ratio+1) * HX.tube_diameter) > HX.duct_size:
        constraint_weight += 1e30
    
    if (HX.duct_length + HX.diffuser_length + HX.nozzle_length) > 10:
        constraint_weight += np.power((HX.duct_length + HX.diffuser_length + HX.nozzle_length - 10), 2) * 1e10
    
    if relative_pd > 20:
        constraint_weight += np.power(relative_pd-20, 2) * 1e20
    """if relative_pd <= 0:
        constraint_weight += 1e30
    elif relative_pd > 10:
        constraint_weight += np.power(relative_pd - 10, 4) * 1e11"""

    if output:
        return(hx_weight + ESM_power_weight + ESM_volume_weight)
    else:
        return(hx_weight + ESM_power_weight + ESM_volume_weight + constraint_weight)

def make_heat_exchanger_pd(input_args, *args, output=False):
    (bulk_area_ratio, coolant_velocity, gap_diameter_ratio) = input_args
    (cooling_power, gasflow, inlet_area, target_pressure_drop_relative, ESM_vector, coolant) = args

    limit_weight = 0
    if bulk_area_ratio < 1.01:
        limit_weight += np.power(bulk_area_ratio-1.01, 2)*1e50
    elif bulk_area_ratio > 15:
        limit_weight += np.power(bulk_area_ratio-50, 2)*1e50
    
    if coolant_velocity < 0.1:
        limit_weight += np.power(coolant_velocity-0.1, 2)*1e50
    elif coolant_velocity > 5:
        limit_weight += np.power(coolant_velocity-5, 2)*1e50
    
    if gap_diameter_ratio < 0.1:
        limit_weight += np.power(gap_diameter_ratio-0.1, 2)*1e50
    elif gap_diameter_ratio > 50:
        limit_weight += np.power(gap_diameter_ratio-50, 2)*1e50
    
    if limit_weight != 0:
        return(limit_weight)

    HX = component.HeatExchangerAdvanced(cooling_power, gasflow, inlet_area, 
        bulk_area_ratio, pitch_diameter_ratio=gap_diameter_ratio+1, coolant_velocity=coolant_velocity, coolant=coolant)

    HX.gas_pressure_drop()
    HX.coolant_pumping_power()
    HX.weight_estimate()

    target_PD = target_pressure_drop_relative * gasflow.pressure

    hx_weight = HX.weight

    hx_volume = (HX.duct_length + 2*HX.diffuser_length) * bulk_area_ratio * inlet_area 

    constraint_weight = 0

    if HX.pressure_drop > target_PD:
        constraint_weight += np.power(target_PD - HX.pressure_drop, 2) * 1e10
    
    if HX.pressure_drop > gasflow.pressure:
        constraint_weight += np.power(HX.pressure_drop - gasflow.pressure, 2) * 1e50

    #constraint_weight += 100*(HX.pressure_drop/gasflow.pressure) * 80 
    
    # if ((gap_diameter_ratio+1) * HX.tube_diameter) > HX.duct_size:
    #     constraint_weight += 1e30
    
    # if (HX.duct_length + HX.diffuser_length + HX.nozzle_length) > 3:
    #     constraint_weight += np.power((HX.duct_length + HX.diffuser_length + HX.nozzle_length - 3), 2) * 1e10
    
    if output:
        return((hx_weight*ESM_vector[0]) + (HX.pumping_power*ESM_vector[1]) + (cooling_power*ESM_vector[2]) +(hx_volume*ESM_vector[3]))
    else:
        return((hx_weight*ESM_vector[0]) + (HX.pumping_power*ESM_vector[1]) + (cooling_power*ESM_vector[2]) +(hx_volume*ESM_vector[3]) + constraint_weight)


def optimise_axial(flow_in, pressure_ratio_stage, efficiency_guess, ESM_vector):
    specific_speed = 2.0
    specific_diameter = 2.0

    deltah = flow_in.delta_h_PR(pressure_ratio_stage) / efficiency_guess

    speed_guess = specific_speed / (np.sqrt(flow_in.mass_flow / flow_in.density_deltah(deltah)) * np.power(deltah, -0.75)) 
    work_coeff_guess = 1 / np.power(specific_speed * specific_diameter, 2)
    flow_coeff_guess = 1 / (specific_speed * np.power(specific_diameter, 3))

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
        args=(flow_in, pressure_ratio_stage, efficiency_guess, ESM_vector),
        method="Nelder-Mead", options={"maxiter":10000, "initial_simplex":sim}
        )


    (speed, work_coeff, flow_coeff) = properties.x

    axial = component.AxialStage(flow_in, pressure_ratio_stage, speed, work_coeff, flow_coeff, efficiency_guess)
    axial.weight_estimate()
    axial.estimate_efficiency()
    axial.estimate_ESM()

    return(axial)

def optimise_radial(flow_in, pressure_ratio_stage, efficiency_guess, ESM_vector):
    specific_speed = 0.5
    specific_diameter = 5.0

    deltah = flow_in.delta_h_PR(pressure_ratio_stage) / efficiency_guess

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
        args=(flow_in, pressure_ratio_stage, efficiency_guess, ESM_vector),
        bounds=[(1, 100e3), (0.1, 0.85), (0.2, 2), (10e-3, 1), (1, 15)],
        method="Nelder-Mead", options={"maxiter": 10000, "disp":False, "adaptive":False, "return_all":False, "xatol":1e-4, "initial_simplex":sim}
        #method='SLSQP'
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":2000}
        )
    
    #print(properties)
    
    (speed, work_coeff, flow_coeff, inlet_hub_radius, diffusion_ratio) = properties.x

    radial = component.RadialStage(flow_in, pressure_ratio_stage, speed, work_coeff, flow_coeff, efficiency_guess, inlet_hub_radius, diffusion_ratio)
    radial.weight_estimate(381e6, 2810)
    radial.estimate_efficiency()
    radial.estimate_ESM()

    return(radial)


def optimise_hx_anypd(flow_in, cooling_deltah, inlet_area):
    
    bulk_area_ratio_guess = 2 
    coolant_velocity_guess = 0.2
    gap_diameter_ratio_guess = 2

    cooling_power = flow_in.mass_flow * cooling_deltah


    x0 = np.array([bulk_area_ratio_guess, coolant_velocity_guess, gap_diameter_ratio_guess])
    x0_upper = np.array([15, 1, 50])
    N = len(x0)
    zdelt = 0.00025
    sim = np.empty((N + 1, N), dtype=x0.dtype)
    sim[0] = x0
    for k in range(N):
        y = 0.9*np.array(x0, copy=True)
        if y[k] != 0:
            y[k] = x0_upper[k]
        else:
            y[k] = zdelt
        sim[k + 1] = y
    #print(sim)

    properties = optimize.minimize(make_heat_exchanger_anypd, 
        x0, 
        args=(cooling_power, flow_in, inlet_area),
        bounds=((1, 15), (0.1, 3), (0.1, 50)),
        method="Nelder-Mead", 
        options={"maxiter": 10000, "disp":False, "adaptive":False, "return_all":False, 
            "xatol":1e-4, "initial_simplex":sim, "bounds":((1, 15), (0.1, 3), (0.1, 50))}
        #method='SLSQP'
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":2000}
        )
    
    print(properties)
    
    (bulk_area, coolant_velocity, gap_diameter_ratio) = properties.x

    HX = component.HeatExchangerAdvanced(cooling_power, flow_in, inlet_area, 
        bulk_area, pitch_diameter_ratio=gap_diameter_ratio+1, coolant_velocity=coolant_velocity, coolant="ammonia")
    HX.gas_pressure_drop()
    HX.coolant_pumping_power()
    HX.weight_estimate()
    HX.estimate_ESM()

    return(HX)

def optimise_hx_pd(flow_in, cooling_deltah, inlet_area, target_PD, ESM_vector, coolant="ammonia"):
    
    bulk_area_guess = 1.5
    coolant_velocity_guess = 0.2
    gap_diameter_ratio_guess = 0.5

    cooling_power = flow_in.mass_flow * cooling_deltah


    x0 = np.array([bulk_area_guess, coolant_velocity_guess, gap_diameter_ratio_guess])
    x0_upper = np.array([10, 2, 50])
    N = len(x0)
    zdelt = 0.00025
    sim = np.empty((N + 1, N), dtype=x0.dtype)
    sim[0] = x0
    for k in range(N):
        y = np.array(x0, copy=True)
        if y[k] != 0:
            y[k] = x0_upper[k]
        else:
            y[k] = zdelt
        sim[k + 1] = y
    #print(sim)

    properties = optimize.minimize(make_heat_exchanger_pd, 
        x0, 
        args=(cooling_power, flow_in, inlet_area, target_PD, ESM_vector, coolant) ,
        method="Nelder-Mead", options={"maxiter": 20000, "disp":False, "adaptive":True, 
        "return_all":False, "initial_simplex":sim}
        #method='SLSQP'
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":2000}
        )

    #print(properties)
    
    
    (bulk_area, coolant_velocity, gap_diameter_ratio) = properties.x

    HX = component.HeatExchangerAdvanced(cooling_power, flow_in, inlet_area, 
        bulk_area, pitch_diameter_ratio=gap_diameter_ratio+1, coolant_velocity=coolant_velocity, coolant=coolant)
    HX.gas_pressure_drop()
    HX.coolant_pumping_power()
    HX.weight_estimate()
    HX.estimate_ESM()

    if properties.success == False:
        print(properties)
        print(HX.pressure_drop/HX.airflow.pressure, HX.ESM)
    
    HX.design_T_in = flow_in.temperature
    HX.design_delta_T = flow_in.delta_T_delta_h(cooling_deltah)

    return(HX)

