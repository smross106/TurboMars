from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import numpy as np
import math
import matplotlib.pyplot as plt

T_K = PropsSI("Hmass", "T", 260, "Q", 1, "CO2") - PropsSI("Hmass", "T", 260, "Q", 0, "CO2")
print(T_K)

exit()

CP.apply_simple_mixing_rule('CarbonDioxide', 'Water', 'linear')
#print(PropsSI("Hmass", "P", 50*1e5, "T", 800, "CarbonDioxide[0.6]&Water[0.4]"))

def temperature_line_water(p_bar):

    # Below critical pressure, use CoolProp
    T_K = PropsSI("T", "P", p_bar*1e5, "Q", 1, "water")
    
    return(T_K)


def T_safety_water(p_bar):
    limit_T_K = temperature_line_water(p_bar)
    mach_correction = 1 + 0.5*(1.4-1)
    T_safe_K = 5 + (limit_T_K * mach_correction)

    return(T_safe_K)

def temperature_line_co2(p_bar):
    if p_bar < 5.19:
        # Below the triple point, use an inverted Antoine equation
        # Curvefit from Vukalovich data
        # Takes data in bar, outputs K
        A = 15.73
        B = 3007.9
        C = -3.386

        T_K = (B / (A - math.log10(p_bar))) - C

    elif p_bar < 73.773:
        # Below critical pressure, use CoolProp
        T_K = PropsSI("T", "P", p_bar*1e5, "Q", 1, "CO2")
    
    else:
        # Above critical pressure, return critical temperature
        T_K = PropsSI("co2", "Tcrit")
    
    return(T_K)


def T_safety_co2(p_bar):
    limit_T_K = temperature_line_co2(p_bar)
    mach_correction = 1 + 0.5*(1.31-1)
    T_safe_K = 10 + (limit_T_K * mach_correction)

    return(T_safe_K)


def pressure_drop_hx(mass_flow, delta_h, pressure_in, high=False):
    Q = mass_flow * delta_h
    # W, so delta_h is J/kg

    if high:
        p_out_bar = pressure_in * (1 - Q * 1.4e-7)
    else:
        p_out_bar = pressure_in * (1 - Q * 7e-8 )
    
    return(p_out_bar)

def weight_turbine():
    return(10)

def weight_hx(mass_flow, delta_h, high=False):
    Q = mass_flow * delta_h
    # W, so delta_h is J/kg

    if high: 
        weight = 0.008 * Q
    else:
        weight = 0.004 * Q
    
    return(weight)

def mixture_h(p_bar, T_k, mass_ratio_water, mass_ratio_co2):
    
    co2_h = PropsSI("Hmass", "P", p_bar*1e5, "T", T_k, "co2")
    h20_h = PropsSI("Hmass", "P", p_bar*1e5, "T", T_k, "water")
    return(mass_ratio_co2*co2_h + mass_ratio_water*h20_h)

def mixture_gamma(p_bar, T_k, mass_ratio_water, mass_ratio_co2):
    co2_cp = PropsSI("Cpmass", "P", p_bar*1e5, "T", T_k, "co2")
    h20_cp = PropsSI("Cpmass", "P", p_bar*1e5, "T", T_k, "water")
    co2_cv = PropsSI("Cvmass", "P", p_bar*1e5, "T", T_k, "co2")
    h20_cv = PropsSI("Cvmass", "P", p_bar*1e5, "T", T_k, "water")
    cp = co2_cp * mass_ratio_co2 + h20_cp * mass_ratio_water
    cv = co2_cv * mass_ratio_co2 + h20_cv * mass_ratio_water
    return(cp/cv)

def turbine_string(reactant_mass_flow, combustor_p, combustor_T, exit_p=0.015, turbine_max_PR=4, turbine_eta=0.8):
    # ALL SPECIFIC WORK SCALED TO REACTANT MASS FLOW


    # Use carbon dioxide dilution to decrease maximum temperature
    combustor_feed_T = 500
    energy_release_fuel = reactant_mass_flow * 0.2 * 50e6 # MJ/kg
    enthalpy_product_heat = (reactant_mass_flow * 0.55*(
        PropsSI("Hmass", "P", combustor_p*1e5, "T", combustor_T, "co2") - 
        PropsSI("Hmass", "P", 101325, "T", 298, "co2")) + 
        reactant_mass_flow * 0.45 * 
        (PropsSI("Hmass", "P", combustor_p*1e5, "T", combustor_T, "water") - 
        PropsSI("Hmass", "P", 101325, "T", 298, "water")))
    
    co2_mass_flow = (energy_release_fuel - enthalpy_product_heat) / (
        PropsSI("Hmass", "P", combustor_p*1e5, "T", combustor_T, "co2") - 
        PropsSI("Hmass", "P", combustor_p*1e5, "T", combustor_feed_T, "co2")
    )

    total_mass_flow = co2_mass_flow + reactant_mass_flow
    reactant_moles = reactant_mass_flow/26.67
    co2_moles = co2_mass_flow/44
    water_mass_fraction = (0.45 * reactant_mass_flow) / total_mass_flow
    co2_mass_fraction = (0.55 * reactant_mass_flow + co2_mass_flow) / total_mass_flow

    co2_acquisition_weight = 1263.8 + 263.57*np.log(co2_mass_flow)
    co2_acquisition_specific_work_J = (416.93 - 21.43*np.log(co2_mass_flow))*1000  * (co2_mass_flow/reactant_mass_flow)
    co2_acquisition_cooling_W = 1000*(1.495 + 783.8*co2_mass_flow)


    # Assume k=0.8 scaling
    reactant_acquisition_weight = 1550 * np.power(reactant_mass_flow/8.44e-4, 0.8)
    reactant_acquisition_specific_work_J = 50e6
    reactant_acquisition_cooling_W = 200e3 * (reactant_mass_flow/0.0211)

    reactant_pump_specific_work = (0.2*(
        PropsSI("Hmass", "P", combustor_p*1e5, "T", 120, "methane") - PropsSI("Hmass", "P", 2e5, "Q", 0, "methane")) + 
        0.8*(PropsSI("Hmass", "P", combustor_p*1e5, "T", 105, "oxygen") - PropsSI("Hmass", "P", 2e5, "Q", 0, "oxygen")) 
    )
    co2_pump_specific_work = (co2_mass_flow/reactant_mass_flow)*(
        PropsSI("Hmass", "P", combustor_p*1e5, "T", 225, "co2") - PropsSI("Hmass", "P", 6e5, "Q", 0, "co2"))
    
    pump_power = (reactant_pump_specific_work + co2_pump_specific_work) * reactant_mass_flow

    total_turbine_weight = 0
    n_turbines = 0
    total_turbine_specific_work = 0
    total_turbine_power = 0

    p_in = combustor_p
    T_in = combustor_T
    h_in = mixture_h(p_in, T_in, water_mass_fraction, co2_mass_fraction)

    # Turbines operating on both CO2 and steam
    reset=True
    while p_in > 2*exit_p:
        if reset:
            guess_PR = turbine_max_PR
        p_turbine = p_in / guess_PR
        gamma = mixture_gamma(p_in, T_in, water_mass_fraction, co2_mass_fraction)

        T_turbine = T_in * (1 + turbine_eta*(np.power(guess_PR, (1-gamma)/gamma)) - 1)
        T_limit = T_safety_water(p_turbine)

        if T_turbine > T_limit:
            print("delta T mix", int(T_in), int(T_turbine))
            p_in = p_turbine
            T_in = T_turbine
            h_out = mixture_h(p_turbine, T_turbine, water_mass_fraction, co2_mass_fraction)
            n_turbines += 1
            total_turbine_weight += weight_turbine()
            total_turbine_specific_work += (h_in - h_out) * (total_mass_flow/reactant_mass_flow)
            total_turbine_power += (h_in - h_out) * total_mass_flow
            print("delta H mix", int((h_in - h_out)), "\n")
            h_in = h_out
            
            reset=True
        else:
            if guess_PR >= 1.6:
                guess_PR /= 1.25
                reset=False
            else:
                break
    
    # Heat exchanger to condense out the water

    delta_h = PropsSI("Hmass", "P", p_in*1e5, "T", T_in, "water") - PropsSI("Hmass", "P", p_in*1e5, "Q", 0, "water")
    condenser_weight = weight_hx(total_mass_flow*water_mass_fraction, delta_h, high=True)
    p_in = pressure_drop_hx(total_mass_flow*water_mass_fraction, delta_h, p_in, high=True)
    T_in = PropsSI("T", "P", p_in*1e5, "Q", 0, "water")
    h_in = PropsSI("Hmass", "P", p_in*1e5, "T", T_in, "co2")


    # Turbines operating on CO2 only
    reset=True
    while p_in > 2*exit_p:
        if reset:
            guess_PR = turbine_max_PR
        p_turbine = p_in / guess_PR
        gamma = PropsSI("Cpmass", "P", p_in*1e5, "T", T_in, "co2")/PropsSI("Cvmass", "P", p_in*1e5, "T", T_in, "co2")

        T_turbine = T_in * (1 + turbine_eta*(np.power(guess_PR, (1-gamma)/gamma)) - 1)
        T_limit = T_safety_co2(p_turbine)

        if T_turbine > T_limit:
            print("delta T co2", int(T_in), int(T_turbine))
            p_in = p_turbine
            T_in = T_turbine
            h_out = PropsSI("Hmass", "P", p_in*1e5, "T", T_in, "co2")
            n_turbines += 1
            total_turbine_weight += weight_turbine()
            total_turbine_specific_work += (h_in - h_out) * (total_mass_flow/reactant_mass_flow)
            total_turbine_power += (h_in - h_out) * total_mass_flow
            reset=True
            print("delta H co2", int((h_in - h_out)), "\n")
            h_in = h_out
        else:
            if guess_PR >= 1.6:
                guess_PR /= 1.25
                reset=False
            else:
                print(p_turbine, T_turbine, T_limit)
                break
            

    #print(p_in, T_in)
    print(co2_mass_flow, reactant_mass_flow)
    print(total_turbine_power-pump_power)
    #print(total_turbine_specific_work)
    #print(n_turbines)
    print(total_turbine_weight, condenser_weight, co2_acquisition_weight, reactant_acquisition_weight)
    print((co2_acquisition_cooling_W+reactant_acquisition_cooling_W))
    #print((total_turbine_weight+condenser_weight+co2_acquisition_weight+reactant_acquisition_weight)+(
    #    co2_acquisition_cooling_W+reactant_acquisition_cooling_W)*.0605)
    #print(total_turbine_specific_work / (reactant_acquisition_specific_work_J + co2_acquisition_specific_work_J+co2_pump_specific_work+reactant_pump_specific_work))

        

turbine_string(0.10824, 100, 1950, turbine_max_PR=4, turbine_eta=0.85)