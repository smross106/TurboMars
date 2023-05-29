import numpy as np
import time

def approximate_radial(inlet_pressure, pressure_ratio, mass_flow, conservative=True):
    # Correlated for a range 1.6kPa < inlet pressure < 600kPa
    # Correlated for a range 1.22 < PR < 1.82
    # Correlated for a range 0.015kg/s < mass flow < 0.1kg/s

    # Coefficients found using statsmodel.ols ANOVA
    # 474 observations
    # All p-values for coefficients <0.001
    # Third quartile factor is the increase needed to ensure 75% of training points are below the curve

    # F-statistic = 2.83e-71
    approximate_mass = np.exp(5.9520
                              - 1.101e-5 * inlet_pressure
                              + 15.1466 * mass_flow
                              + 1.025e-5 * inlet_pressure * pressure_ratio
                              - 0.4972 * np.log(inlet_pressure)
                              )
    
    # F-statistic = 2.54e-98
    approximate_efficiency = (0.4413 
                              + 2.02e-6 * inlet_pressure
                              - 2.4722 * mass_flow  
                              - 2.714e-6 * inlet_pressure * pressure_ratio
                              + 2.6277 * pressure_ratio * mass_flow
                              + 1.094e-5 * inlet_pressure * pressure_ratio * mass_flow
                              + 0.0255 * np.log(inlet_pressure)
                              )
    
    third_quartile_factor_mass = 1.6827 
    third_quartile_factor_efficiency = 1 - 0.0829
           
    if conservative:
        approximate_mass *= third_quartile_factor_mass
        approximate_efficiency *= third_quartile_factor_efficiency
    
    if approximate_efficiency < 0.01:
        approximate_efficiency = 0.01

    
    return(approximate_mass, approximate_efficiency)

def approximate_hx(inlet_pressure, cooling_power, mass_flow, conservative=True):
    # Correlated for a range 2kPa < inlet pressure < 700kPa
    # Correlated for a range 157W < cooling power < 4406W
    # Correlated for a range 0.015kg/s < mass flow < 0.1kg/s

    # Coefficients found using statsmodel.ols ANOVA
    # 98 observations
    # All p-values for coefficients <0.002
    # Third quartile factor is the increase needed to ensure 75% of training points are below the curve

    # F-statistic = 2.82e-51
    approximate_mass = np.exp(-2.2736
                              - 0.2220 * np.log(inlet_pressure)
                              + 0.7122 * np.log(cooling_power)
    )

    # F-statistic = 2.45e-39
    approximate_pressure_drop_ratio = np.exp(-8.2904
                                             + 0.2856 * np.log(cooling_power)
                                             + 2.7959 * mass_flow
                                             - 0.1639 * np.log(inlet_pressure))
    
    third_quartile_factor_mass = 1.1726 
    third_quartile_factor_pressure_drop_ratio = 1.1478
           
    if conservative:
        approximate_mass *= third_quartile_factor_mass
        approximate_pressure_drop_ratio *= third_quartile_factor_pressure_drop_ratio
    
    return(approximate_mass, approximate_pressure_drop_ratio)

mass_flows = [0.01, 0.03, 0.075, 0.12]
pressure_ratios = [2.115, 1.187]
N_stag = [8, 35]

times_list = []
for mass_flow in mass_flows:
    for pr_i, pressure_ratio in enumerate(pressure_ratios):
        N_stages = N_stag[pr_i]
#mass_flow = 0.075
#pressure_ratio = 1.272
#N_stages = 25

        total_mass = 0
        total_efficiency = 0
        total_Q = 0
        starting_pressure = 1600
        pressure = starting_pressure

        start = time.time()

        for i in range(N_stages):
            rad_mass, rad_eff = approximate_radial(pressure, pressure_ratio, mass_flow,conservative=True)
            total_mass += rad_mass
            total_efficiency += rad_eff

            if rad_eff < 0.5:
                rad_eff = 0.5

            Q = 250 * mass_flow * 830 * (np.power(pressure_ratio, 0.31/1.31)-1) / rad_eff
            total_Q += Q
            hx_mass, hx_pr = approximate_hx(pressure * pressure_ratio, Q, mass_flow, conservative=True)
            total_mass += hx_mass
            pressure *= (pressure_ratio * (1-hx_pr))

            #print(Q)
            
        times_list.append( time.time() - start)
        print(total_mass, "\t", total_Q, "\t", total_efficiency/N_stages, "\t", pressure/1000)

print("###")
for i in times_list:print(i)