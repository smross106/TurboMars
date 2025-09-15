from CoolProp.CoolProp import PropsSI
import numpy as np
import math
import matplotlib.pyplot as plt

GAMMA = 1.31

def pressure_line(temperature):
    if temperature < 216.593:
        # Below the triple point, use quartic fit to Vukalovich data
        p_bar = ((1.7095E-7 * np.power(temperature, 4)) -
                 (9.3886E-5 * np.power(temperature, 3)) +
                 (1.9028E-2 * np.power(temperature, 2)) -
                 (1.6831E0  * np.power(temperature, 1)) +
                  5.4740E1)
    elif temperature < 304.1282:
        p = PropsSI("P", "T", temperature, "Q", 1, "CO2")
        p_bar = p/1e5
    else:
        p_bar = PropsSI("co2", "pcrit")/1e5
    
    return(p_bar)

def temperature_line(p_bar):
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


def T_safety(p_bar):
    limit_T_K = temperature_line(p_bar)
    mach_correction = 1 + 0.5*(GAMMA-1)
    T_safe_K = 5 + (limit_T_K * mach_correction)

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

def generate_machine(mass_flow_in, mass_flow_out, pump_p, reheat_T, liquid_p=5.5, exit_p=0.015, turbine_max_PR=4, turbine_eta=0.8):
    compressor_weight_kg = 1263.8 + 263.57*np.log(mass_flow_in)
    compressor_specific_work_J = (416.93 - 21.43*np.log(mass_flow_in))*1000
    compressor_cooling_kW = 1.495 + 783.8*mass_flow_in

    liquid_enthalpy = PropsSI("Hmass", "P", liquid_p*1e5, "Q", 0, "co2")
    liquid_T = PropsSI("T", "P", liquid_p*1e5, "Q", 0, "co2")
    pump_enthalpy = PropsSI("Hmass", "P", pump_p*1e5, "T|liquid", liquid_T+5, "co2")

    pump_specific_work_J = pump_enthalpy - liquid_enthalpy

    # Boiler
    T_in = liquid_T+5
    p_in = pump_p
    boiler_T_out = reheat_T
    h_in = pump_enthalpy
    approx_h_out = PropsSI("Hmass", "P", p_in*1e5, "T", boiler_T_out, "co2")
    boiler_p_out = pressure_drop_hx(mass_flow_out, (approx_h_out-h_in), p_in, high=True)
    boiler_h_out = PropsSI("Hmass", "P", boiler_p_out*1e5, "T", boiler_T_out, "co2")
    boiler_specific_heat = (boiler_h_out - h_in)
    boiler_weight = weight_hx(mass_flow_out, (approx_h_out-h_in), high=True)

    optimal_PR = turbine_max_PR - 0.5
    

    for i in range(7):
        n_turbines = 0
        n_HX = 0
        total_turbine_work = 0
        total_turbine_specific_work = 0
        total_HX_specific_work = 0
        total_turbine_weight = 0
        total_HX_weight = 0
        average_heat_addition_T = 0

        # Loop to optimise optimal pressure ratio in turbines
        p_in = boiler_p_out
        T_in = boiler_T_out
        h_in = PropsSI("Hmass", "P", p_in*1e5, "T", T_in, "co2")
        while p_in > exit_p:
            # Create a turbine-HX pair

            # Trial a turbine with the optimal PR, if it causes an undertemperature then decrease it
            PR_guess = optimal_PR
            for PR_try in range(3):
                p_out = p_in/PR_guess
                isentropic_h_out = PropsSI("Hmass", "P", p_out*1e5, "Smass", PropsSI("Smass", "P", p_in*1e5, "T", T_in, "co2"), "co2")
                real_h_out = h_in - turbine_eta*(h_in - isentropic_h_out)
                T_out = PropsSI("T", "Hmass", real_h_out, "P", p_out*1e5, "co2")
                T_limit = T_safety(p_out)

                if T_out < T_limit and PR_guess > 1.1:
                    PR_guess *= 1 - (1/(PR_try+2))
                    PR_guess = max(PR_guess, 1.1)
                elif T_out > T_limit and PR_guess < optimal_PR:
                    PR_guess *= 1 + (1/(PR_try+2))
                    PR_guess = min(PR_guess, optimal_PR)
                elif T_out > T_limit and PR_guess >= optimal_PR:
                    break
                elif T_out < T_limit and abs(PR_guess - 1.1) < 0.05:
                    #print("Pressure too high and temperature too low, break")
                    raise ValueError("Reduce pump pressure, increase reheat temperature")
            
            total_turbine_work += mass_flow_out * (h_in - real_h_out)
            total_turbine_specific_work += (h_in - real_h_out)
            total_turbine_weight += weight_turbine()
            n_turbines += 1
            h_in = real_h_out
            p_in = p_out
            T_in = T_out

            # Add a heat exchanger if above the target pressure
            if p_in >= exit_p*1.5:
                T_out = reheat_T
                approx_h_out = PropsSI("Hmass", "P", p_in*1e5, "T", T_out, "co2")
                p_out = pressure_drop_hx(mass_flow_out, (approx_h_out-h_in), p_in, high=True)
                h_out = PropsSI("Hmass", "P", p_out*1e5, "T", T_out, "co2")
                total_HX_weight += weight_hx(mass_flow_out, h_out-h_in)
                total_HX_specific_work += (h_out - h_in)
                average_heat_addition_T += 0.5*(T_in + T_out) * (h_out - h_in)
                n_HX += 1

                h_in = h_out
                p_in = p_out
                T_in = T_out

            #print(PR_guess, Q_hx, T_in, T_out, p_out)

            

        # Adjust optimal PR to hit the ideal exit pressure
        target_overall_PR = pump_p / exit_p
        target_machine_PR = np.power(target_overall_PR, 1/n_turbines)
        if target_machine_PR > turbine_max_PR:
            target_machine_PR = np.power(target_overall_PR, 1/(n_turbines+1))
        
        optimal_PR = target_machine_PR

    # Average heat addition from boiler
    average_heat_addition_T += 0.5 * (liquid_T+5 + PropsSI("T", "P", pump_p*1e5, "Q", 0, "co2")) * (
        PropsSI("T", "P", pump_p*1e5, "Q", 0, "co2")-pump_enthalpy)
    average_heat_addition_T += PropsSI("T", "P", pump_p*1e5, "Q", 0, "co2") * (
        PropsSI("T", "P", pump_p*1e5, "Q", 1, "co2")-PropsSI("T", "P", pump_p*1e5, "Q", 0, "co2"))
    
    average_heat_addition_T += 0.5 * (reheat_T + PropsSI("T", "P", pump_p*1e5, "Q", 0, "co2")) * (
        boiler_h_out - PropsSI("T", "P", pump_p*1e5, "Q", 1, "co2")
    )
            
    system_weight = boiler_weight + total_turbine_weight + total_HX_weight
    system_power = total_turbine_specific_work * mass_flow_out
    net_system_power = system_power - (pump_specific_work_J*mass_flow_out)
    heat_in = (boiler_specific_heat+total_HX_specific_work)*mass_flow_out
    average_heat_addition_T /= (boiler_specific_heat+total_HX_specific_work)

    #print(system_power/system_weight)
    #print(system_weight, boiler_weight, total_HX_weight)
    #print(system_power)
    #print((boiler_specific_heat+total_HX_specific_work)*mass_flow_out)
    print(n_turbines)
    print(system_weight)
    print(compressor_weight_kg)
    print(compressor_cooling_kW*121)

    output_dict = {
        "Net power out": net_system_power,
        "Thermal power in": heat_in,
        "Thermal efficiency": (net_system_power-compressor_specific_work_J*mass_flow_out)/heat_in,
        "Round trip efficiency": (total_turbine_specific_work-pump_specific_work_J)/compressor_specific_work_J,
        "Specific power": net_system_power / (system_weight + compressor_weight_kg), 
        "Specific power with cooling": net_system_power / (system_weight + compressor_weight_kg + compressor_cooling_kW*60.5), 
        "Temperature of heat addition": average_heat_addition_T
    }
    print(output_dict)
    return(output_dict)
    #print(total_turbine_specific_work-pump_specific_work_J, compressor_specific_work_J, (total_turbine_specific_work-pump_specific_work_J)/compressor_specific_work_J)


#output = generate_machine(2.82, 2.82, 50, 340, turbine_max_PR=4)
#output = generate_machine(0.2083, 0.2083, 5.5, 340, turbine_max_PR=4)

print(PropsSI("T", "P", 2500e3, "Q", 0, "CO2"))
print(PropsSI("Hmass", "P", 600e3, "Q", 1, "CO2") - PropsSI("Hmass", "P", 600e3, "Q", 0, "CO2"))
print(PropsSI("Dmass", "P", 600e3, "Q", 0, "CO2"))
print(PropsSI("Hmass", "P", 100e3, "Q", 1, "oxygen") - PropsSI("Hmass", "P", 100e3, "Q", 0, "oxygen"))
print(PropsSI("Dmass", "P", 100e3, "Q", 0, "oxygen"))
print(PropsSI("Hmass", "P", 100e3, "Q", 1, "methane") - PropsSI("Hmass", "P", 100e3, "Q", 0, "methane"))
print(PropsSI("Dmass", "P", 100e3, "Q", 0, "methane"))

exit()
#reheat_T_range = np.linspace(273, 400, 50)
reheat_T_range = np.concatenate((np.arange(273, 360, 0.25), np.arange(360, 400, 0.5)))
pump_p_range = np.linspace(10, 70, 125)

RTE_range = np.zeros((reheat_T_range.size, pump_p_range.size))
specific_power_range = np.zeros((reheat_T_range.size, pump_p_range.size))
power_range = np.zeros((reheat_T_range.size, pump_p_range.size))
input_power_range = np.zeros((reheat_T_range.size, pump_p_range.size))

for rhT_i, reheat_T in enumerate(reheat_T_range):
    for p_i, pump_p in enumerate(pump_p_range):
        try:
            output = generate_machine(0.215, 0.215, pump_p, reheat_T, turbine_max_PR=4)
            RTE_range[rhT_i, p_i] = output["Round trip efficiency"]
            specific_power_range[rhT_i, p_i] = output["Specific power"]
            power_range[rhT_i, p_i] = output["Net power out"]
            input_power_range[rhT_i, p_i] = output["Thermal power in"]
        except:
            RTE_range[rhT_i, p_i] = np.nan
            specific_power_range[rhT_i, p_i] = np.nan
            power_range[rhT_i, p_i] = np.nan
            input_power_range[rhT_i, p_i] = np.nan
    
    print(int(rhT_i/reheat_T_range.size * 100),"%")

fig, ax = plt.subplots()       
cont = ax.contourf(reheat_T_range, pump_p_range, np.transpose(RTE_range))
fig.colorbar(cont, ax=ax)
cont = ax.contour(reheat_T_range, pump_p_range, np.transpose(RTE_range), levels=[1], c="r")
cont2 = ax.contour(reheat_T_range, pump_p_range, np.transpose(input_power_range), levels=[180e3], c="k", linestyles="dotted")
ax.set_xlabel("Reheat temperature (K)")
ax.set_ylabel("Pump pressure (bar)")
ax.set_title("Round-trip efficiency")

fig, ax = plt.subplots()       
cont = ax.contourf(reheat_T_range, pump_p_range, np.transpose(specific_power_range))
fig.colorbar(cont, ax=ax)
cont2 = ax.contour(reheat_T_range, pump_p_range, np.transpose(input_power_range), levels=[180e3], c="k", linestyles="dotted")
ax.set_xlabel("Reheat temperature (K)")
ax.set_ylabel("Pump pressure (bar)")
ax.set_title("Specific power (W/kg)")

fig, ax = plt.subplots()       
cont = ax.contourf(reheat_T_range, pump_p_range, np.transpose(power_range))
fig.colorbar(cont, ax=ax)
cont2 = ax.contour(reheat_T_range, pump_p_range, np.transpose(input_power_range), levels=[180e3], c="k", linestyles="dotted")
ax.set_xlabel("Reheat temperature (K)")
ax.set_ylabel("Pump pressure (bar)")
ax.set_title("Power output (W)")

fig, ax = plt.subplots()       
cont = ax.contourf(reheat_T_range, pump_p_range, np.transpose(input_power_range))
fig.colorbar(cont, ax=ax)
cont2 = ax.contour(reheat_T_range, pump_p_range, np.transpose(input_power_range), levels=[180e3], c="k", linestyles="dotted")
ax.set_xlabel("Reheat temperature (K)")
ax.set_ylabel("Pump pressure (bar)")
ax.set_title("Power input (W)")

plt.show()

# https://ntrs.nasa.gov/api/citations/20040010319/downloads/20040010319.pdf