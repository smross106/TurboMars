import component
import design


import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import time
from math import log
import copy

def run_machine(machines, T_start, p_start, mass_flow, label="", verbose=True):
    rejection_T = False
    
    running_machines = copy.deepcopy(machines)
    T_in = T_start
    p_in = p_start

    compressor_type_numbers = [1, 1, 1]

    machine_label = ["In"]
    compressor_label = []
    HX_label = []

    machines_power = []
    machines_cooling = []
    machines_reject_temperature_cooling = []
    pressures = [p_start]
    temperatures = [T_start]
    machines_efficiency = []

    
    for machine in running_machines:
        flow = component.GasFlow(mass_flow, T_in, p_in)

        machine.gasflow_in = flow
        
        if type(machine) == component.AxialStage or type(machine) == component.RadialStage or type(machine) == component.ScrollCompressor:
            machine.estimate_efficiency()
            machine.estimate_ESM()

            machines_power.append(machine.power)
            machines_efficiency.append(machine.efficiency)

            if type(machine) == component.AxialStage:
                machine_label.append(str("A"+str(compressor_type_numbers[0])))
                compressor_label.append(str("A"+str(compressor_type_numbers[0])))
                compressor_type_numbers[0] += 1
            elif type(machine) == component.RadialStage:
                machine_label.append(str("R"+str(compressor_type_numbers[1])))
                compressor_label.append(str("R"+str(compressor_type_numbers[1])))
                compressor_type_numbers[1] += 1
            elif type(machine) == component.ScrollCompressor:
                machine_label.append(str("S"+str(compressor_type_numbers[2])))
                compressor_label.append(str("S"+str(compressor_type_numbers[2])))
                compressor_type_numbers[2] += 1
            
            p_in = p_in * machine.pressure_ratio_stage
            T_in = flow.temperature * ((np.power(machine.pressure_ratio_stage, flow.gm1/flow.gamma) - 1)/max(machine.efficiency, 0.1) + 1)
        
        elif type(machine) == component.HeatExchangerAdvanced:
            machine.gas_pressure_drop()

            # Assume heat transfer scaled with LMTD
            machine_max_cooling_duty = (machine.coolant_max_T_out - machine.coolant_T_in) / machine.design_delta_T
            machine.cooling *= min(((T_in - machine.design_T_in) + machine.coolant_dT) / machine.coolant_dT, machine_max_cooling_duty)

            #print(min(((T_in - machine.design_T_in) + machine.coolant_dT) / machine.coolant_dT, machine_max_cooling_duty), ((T_in - machine.design_T_in) + machine.coolant_dT) / machine.coolant_dT, machine_max_cooling_duty)

            machines_power.append(machine.power)
            machines_cooling.append(machine.cooling)
            machines_reject_temperature_cooling.append(machine.cooling * 0.5 * (T_in + T_in - flow.delta_T_delta_h(machine.cooling/mass_flow)))

            machine_label.append(str("H"+str(len(HX_label)+1)))
            HX_label.append(str("H"+str(len(HX_label)+1)))

            T_in = T_in - flow.delta_T_delta_h(machine.cooling/mass_flow)
            p_in -= machine.pressure_drop
        
        elif type(machine) == component.CondensingHX:
            machine.estimate_pressure_drop()

            machines_power.append(machine.power)
            machines_cooling.append(machine.cooling)
            machines_reject_temperature_cooling.append(machine.cooling * 210)

            machine_label.append(str("COND1"))
            HX_label.append(str("COND1"))

            T_in = T_in
            p_in -= machine.pressure_drop
        
        elif type(machine) == component.BufferTank:
            machines_power.append(machine.power)

            machine_label.append(str("T1"))

        
        pressures.append(p_in)
        temperatures.append(T_in)
    
    #print([int(T) for T in temperatures])
    #print(T_start, "\n", sum(machines_power), "\n", sum(machines_cooling), "\n")
    plt.plot(machine_label, temperatures, label=label)
    if rejection_T:
        return(sum(machines_power), sum(machines_cooling), sum(machines_reject_temperature_cooling)/sum(machines_cooling))
    else:
        return(sum(machines_power), sum(machines_cooling))


def generate_heat_exchanger(flow_hot, deltaH, A_HX_in, maximum_HX_PR, ESM_vector, coolant="ammonia"):
    HX_PR_limit = maximum_HX_PR
    HX_PR_increase_ratio = 1.2
    HX_PR_start = 0.005
    HX_PR_increases_limit = int(np.log(HX_PR_limit/HX_PR_start)/np.log(HX_PR_increase_ratio))

    HX_PR_target = HX_PR_start
    HX_PR_increases = 0
    hx1_PR_actual = 1
    HX_approved = False

    while (not HX_approved) and (HX_PR_increases < HX_PR_increases_limit):
        hx1 = design.optimise_hx_pd(flow_hot, deltaH, A_HX_in, HX_PR_target, ESM_vector, coolant=coolant)
        hx1_PR_actual = hx1.pressure_drop / flow_hot.pressure

        if hx1_PR_actual < HX_PR_target:
            HX_approved = True

        HX_PR_target *= HX_PR_increase_ratio
        HX_PR_increases += 1

    if np.isnan(hx1.pressure_drop) or hx1_PR_actual > maximum_HX_PR:
        print("### \t ERROR \t ###")
        print(hx1)
        print(hx1.bulk_area_ratio, hx1.coolant_velocity, hx1.pitch_diameter_ratio)
    else:
        return(hx1)

def ESM(weight, power, cooling, ESM_vector, verbose=True):
    ESM = (weight*ESM_vector[0]) + (power*ESM_vector[1]) + (cooling*ESM_vector[2])

    return(ESM)


def design_machine(mass_flow, T_in, p_in, ESM_vector, num_intercoolers, rad_PR, p_out_target=555.6e3, condense=True, verbose=True):

    #ESM_vector = [1, 0.149, 0.121, 100]

    #T_in = 190
    #p_in = 781
    #p_out_target = 555.6e3

    total_compression_work = 0
    total_compressor_weight= 0
    total_cooling = 0
    compressor_efficiency = []
    machines_before_scroll = []

    #mass_flow = 0.0113

    N_stages = num_intercoolers

    #rad_PR = 1.222       # For 30 stages
    #rad_PR = 1.272      # for 25 stages
    #rad_PR = 1.3516       # for 20 stages
    #rad_PR = 1.5066       # For 15 stages
    #rad_PR = 1.8304       # For 10 stages

    after_scroll = False
    weight_after_scroll = 0
    power_after_scroll = 0
    cooling_after_scroll = 0
    machines_after_scroll = []
    p_inlet_scroll = 0
    T_inlet_scroll = 0

    start = time.time()

    scroll_compressor = []

    if verbose:print("\t Designing initial axials")
    for i in range(6):
        flow_cold = component.GasFlow(mass_flow, T_in, p_in)
        ax1 = design.optimise_axial(flow_cold, 1.1, 0.85, ESM_vector)
        #print("Axial", str(i+1), ax1.efficiency, ax1.weight)
        compressor_efficiency.append(ax1.efficiency)
        total_compressor_weight += ax1.weight
        total_compression_work += ax1.power

        hot_T = flow_cold.temperature * ((np.power(1.1, flow_cold.gm1/flow_cold.gamma) - 1)/ax1.efficiency + 1)

        T_in = hot_T
        p_in *= 1.1

        machines_before_scroll.append(ax1)

    if verbose:print("\t Designing intercooled radials")
    for i in range(N_stages):
        flow_cold = component.GasFlow(mass_flow, T_in, p_in)

        # ax1 = design.optimise_axial(flow, flow.delta_h_PR(1.1), 0.85)

        rad1 = design.optimise_radial(flow_cold, rad_PR, 0.7, ESM_vector)

        equivalent_scroll = component.ScrollCompressor(flow_cold, rad_PR)
        equivalent_scroll.estimate_ESM()

        equivalent_scroll_ESM_scale = equivalent_scroll.weight*ESM_vector[0] + equivalent_scroll.power*(ESM_vector[1]+ESM_vector[2])
        rad1_ESM_scale = rad1.weight*ESM_vector[0] + rad1.power*(ESM_vector[1]+ESM_vector[2])

        if equivalent_scroll_ESM_scale < rad1_ESM_scale:
            rad1 = equivalent_scroll
            use_scroll = True
            if after_scroll == False:
                after_scroll = True
                p_inlet_scroll = p_in
                T_inlet_scroll = T_in
        else:
            use_scroll = False

        compressor_efficiency.append(rad1.efficiency)
        scroll_compressor.append(use_scroll)

        #if use_scroll:
        #    print("Scroll",str(i+1), rad1.efficiency, rad1.weight)
        #else:
        #    print("Radial",str(i+1), rad1.efficiency, rad1.weight)

        hot_T = flow_cold.temperature * ((np.power(rad_PR, flow_cold.gm1/flow_cold.gamma) - 1)/rad1.efficiency + 1)

        
        flow_hot = component.GasFlow(mass_flow, hot_T, flow_cold.pressure*rad_PR)
        deltaH = flow_hot.delta_h_delta_T(abs(250-flow_hot.temperature)) 

        total_compression_work += rad1.power
        total_compressor_weight += rad1.weight

        if after_scroll:
            power_after_scroll += rad1.power
            weight_after_scroll += rad1.weight
            machines_after_scroll.append(rad1)
        else:
            machines_before_scroll.append(rad1)


        A_HX_in = flow_hot.mass_flow / (flow_hot.density() * 10)

        hx1 = generate_heat_exchanger(flow_hot, deltaH, A_HX_in, 1-(1/rad_PR), ESM_vector)

        total_cooling += hx1.cooling
        total_compression_work += hx1.power
        total_compressor_weight += hx1.weight

        if after_scroll:
            power_after_scroll += hx1.power
            weight_after_scroll += hx1.weight
            cooling_after_scroll += hx1.cooling
            machines_after_scroll.append(hx1)
        else:
            machines_before_scroll.append(hx1)


        HX_PR = hx1.pressure_drop / flow_hot.pressure
        
        #print("HX",str(i+1), HX_PR, hx1.weight)
        
        T_in = 250
        p_in = flow_hot.pressure * (1 - HX_PR)

    p_out_initial = p_in
    T_out_initial = T_in

    best_scroll_machines = machines_after_scroll

    if sum(scroll_compressor) > 1:
        if verbose:print("\t Multiple scrolls detected - running scroll number optimisation pass")
        nominal_scroll_compressors = sum(scroll_compressor)
        total_scroll_compression = np.power(rad_PR, nominal_scroll_compressors)

        best_scroll_N = nominal_scroll_compressors
        best_scroll_weight = weight_after_scroll
        best_scroll_power = power_after_scroll
        best_scroll_cooling = cooling_after_scroll
        best_scroll_ESM = ESM(best_scroll_weight, best_scroll_power, best_scroll_cooling, ESM_vector)
        best_scroll_machines = machines_after_scroll
        p_out_best = p_out_initial
        T_out_best = T_out_initial

        for test_num_scroll_compressors in reversed(range(1,nominal_scroll_compressors, 1)):
            scroll_PR = np.power(total_scroll_compression, 1/test_num_scroll_compressors)

            for i in range(3):
                T_in = T_inlet_scroll
                p_in = p_inlet_scroll
                test_scroll_weight = 0
                test_scroll_power = 0
                test_scroll_cooling = 0
                test_scroll_machines = []

                for stage in range(test_num_scroll_compressors):
                    flow_cold = component.GasFlow(mass_flow, T_in, p_in)
                    
                    scr1 = component.ScrollCompressor(flow_cold, scroll_PR)
                    scr1.estimate_ESM()
                    test_scroll_weight += scr1.weight
                    test_scroll_power += scr1.power
                    test_scroll_machines.append(scr1)

                    hot_T = flow_cold.temperature * ((np.power(scroll_PR, flow_cold.gm1/flow_cold.gamma) - 1)/scr1.efficiency + 1)

                    flow_hot = component.GasFlow(mass_flow, hot_T, flow_cold.pressure*scroll_PR)
                    Q = flow_hot.delta_h_delta_T(abs(250-flow_hot.temperature)) 

                    A_HX_in = flow_hot.mass_flow / (flow_hot.density() * 10)

                    hx1 = generate_heat_exchanger(flow_hot, Q, A_HX_in, 1-(1/scroll_PR), ESM_vector)

                    test_scroll_cooling += hx1.cooling
                    test_scroll_power += hx1.power
                    test_scroll_weight += hx1.weight
                    test_scroll_machines.append(hx1)

                    T_in = 250
                    HX_PR = hx1.pressure_drop / flow_hot.pressure
                    p_in = flow_hot.pressure * (1 - HX_PR)

                scroll_PR *= np.power(p_out_target/p_in, 1/test_num_scroll_compressors)
            
            p_out = p_in
            T_out = T_in

            test_scroll_ESM = ESM(test_scroll_weight, test_scroll_power, test_scroll_cooling, ESM_vector)

            if test_scroll_ESM < best_scroll_ESM and p_out>p_out_target:
                if verbose:print("\t \t", test_num_scroll_compressors, "scrolls beat", best_scroll_N,"scrolls")
                best_scroll_ESM = test_scroll_ESM
                best_scroll_weight = test_scroll_weight
                best_scroll_power = test_scroll_power
                best_scroll_cooling = test_scroll_cooling
                best_scroll_machines = test_scroll_machines
                best_scroll_N = test_num_scroll_compressors  
                #print(p_out_initial, p_out, scroll_PR)
                p_out_best = p_out
                T_out_best = T_out

        if verbose:print("\t \t", best_scroll_N,"scrolls are best")  

    all_machines = machines_before_scroll + best_scroll_machines

    if condense:
        flow = component.GasFlow(mass_flow, T_out_best, p_out_best)
        precondenser = generate_heat_exchanger(flow, flow.delta_h_delta_T(T_out - 220), flow.mass_flow / (flow.density() * 10), 0.1, ESM_vector)

        p_out_best -= precondenser.pressure_drop
        
        condenser_flow_in = component.GasFlow(mass_flow, 220, p_out_best)

        condenser = component.CondensingHX(condenser_flow_in)
        condenser.estimate_ESM()
        condenser.estimate_pressure_drop()
        p_out_best -= condenser.pressure_drop

        buffer_tank = component.BufferTank(1500/1160)
        buffer_tank.estimate_ESM()

        all_machines.append(precondenser)
        all_machines.append(condenser)
        all_machines.append(buffer_tank)

        if verbose:print("\t As-designed radial output pressure {:.1f}kPa".format((p_out_initial-precondenser.pressure_drop-condenser.pressure_drop)/1000))

    total_work = 0
    total_weight = 0
    total_cooling = 0

    num_axials = 0
    num_radials = 0
    num_scrolls = 0
    num_hxs = 0

    for machine in all_machines:
        total_work += machine.power
        total_cooling += machine.cooling
        total_weight += machine.weight
        if type(machine) == component.AxialStage: 
            num_axials += 1
        elif type(machine) == component.RadialStage: 
            num_radials += 1
        elif type(machine) == component.ScrollCompressor: 
            num_scrolls += 1
        elif type(machine) == component.HeatExchangerAdvanced: 
            num_hxs += 1

    exec_time = time.time() - start
    if verbose:
        print("")
        print("Runtime {:.1f}".format(exec_time))
        print("Weight {:.1f}kg".format(total_weight))
        print("Power {:.1f}W-e".format(total_work))
        print("Cooling {:.1f}W-th".format(total_cooling))
        print("Specific compression work {:.1f}kJ/kg".format((total_work/mass_flow)/1000))
        print("Total ESM {:.1f}kg-eq".format(ESM(total_weight, total_work, total_cooling, ESM_vector)))
        print("Output pressure {:.1f}kPa".format(p_out_best/1000))
        print("A {}-stage compressor with {} axial, {} radial and {} scroll compressor stages and {} intercoolers".format(
            num_axials+num_radials+num_scrolls, num_axials, num_radials, num_scrolls, num_hxs))
    
    data_dict = {"Runtime": exec_time,
                 "Weight": total_weight,
                 "Output pressure": p_out_best}
    
    return(all_machines, data_dict)

def run_machines_day(machines, mass_flow, hours, design_T_in, design_p_in, verbose=True):
    T_reject = False


    if verbose:print("Running daily variation")
    all_data = [[]]*8
    if hours[0]:
        all_data[0] = run_machine(machines, 184.83, 775.79, mass_flow, "H0")
    if hours[1]:
        all_data[1] = run_machine(machines, 180.18, 779.94, mass_flow, "H3")
    if hours[2]:
        all_data[2] = run_machine(machines, 184.71, 792.95, mass_flow, "H6")
    if hours[3]:
        all_data[3] = run_machine(machines, 205.68, 796.48, mass_flow, "H9")
    if hours[4]:
        all_data[4] = run_machine(machines, 221.67, 782.96, mass_flow, "H12")
    if hours[5]:
        all_data[5] = run_machine(machines, 222.98, 770.60, mass_flow, "H15")
    if hours[6]:
        all_data[6] = run_machine(machines, 206.66, 773.60, mass_flow, "H18")
    if hours[7]:
        all_data[7] = run_machine(machines, 192.11, 773.11, mass_flow, "H21")

    if T_reject:
        design_pow, design_cool, design_T = run_machine(machines, design_T_in, design_p_in, mass_flow, "Design")
    else:
        design_pow, design_cool = run_machine(machines, design_T_in, design_p_in, mass_flow, "Design")

    num_axials = 0
    num_radials = 0
    num_scrolls = 0
    num_hxs = 0

    total_weight = 0
    
    for machine in machines:
        total_weight += machine.weight
        if type(machine) == component.AxialStage: 
            num_axials += 1
        elif type(machine) == component.RadialStage: 
            num_radials += 1
        elif type(machine) == component.ScrollCompressor: 
            num_scrolls += 1
        elif type(machine) == component.HeatExchangerAdvanced: 
            num_hxs += 1
    
    p_string = ""
    c_string = ""
    if T_reject:
        T_string = ""

    p_string += "{} \t".format(num_hxs-1)
    c_string += " \t"
    p_string += "Power \t"
    c_string += "Cooling \t"
    if T_reject:
        T_string += "ICTreject \t"


    for i in all_data:
        if len(i) != 0:
            p_string += "{:.4f} \t".format(i[0])
            c_string += "{:.4f} \t".format(i[1])
            if T_reject:
                T_string += "{:.2f} \t".format(i[2])
        else:
            p_string += "0 \t"
            c_string += "0 \t"
            if T_reject:
                T_string += "0 \t"

    p_string += "{:.4f} \t".format(design_pow)
    c_string += "{:.4f} \t".format(design_cool)

    p_string += "{:.1f} \t {} \t {} \t {}".format(total_weight, num_axials, num_radials, num_scrolls)

    print(p_string)
    print(c_string)
    if T_reject:
        print(T_string)

def optimise_machine_design(mass_flow, design_T_in, design_p_in, scroll_target_p_out, PR_guess, N_rads, ESM_vector, condense=True, target_p_out=550e3, verbose=True):

    """mass_flow = 0.0113
    design_T_in = 199.85
    design_p_in = 770.6

    PR_guess = 1.84
    target_p_out = 550e3
    scroll_target_p_out = 550e3
    
    N_rads = 10"""
    rad_PR = PR_guess
    if verbose:print("Designing machine")
    for i in range(10):
        
        machines, data_dict = design_machine(mass_flow, design_T_in, design_p_in, ESM_vector, N_rads, rad_PR, condense=condense, p_out_target=scroll_target_p_out, verbose=verbose)
        if abs(data_dict["Output pressure"] -target_p_out) < 0.3e3:
            break
        else:
            if verbose:print("Iteration to improve output pressure")
            rad_PR *= np.power(target_p_out/data_dict["Output pressure"], 1/N_rads)
    
    return(machines)

def optimise_design_system():
    # Design space evaluation
    mass_flows = [0.0113, 0.0225, 0.0338, 0.0451, 0.0902]
    design_T_ins = [199.85, 185.46, 183.04, 182.45, 180.18]
    design_p_ins = [770.6, 775.79, 775.79, 779.94, 779.94]
    hours = [[True]*8,
             [True, True, True, False, False, False, False, True],
             [True, True, True, False, False, False, False, False],
             [False, True, True, False, False, False, False, False],
             [False, True, False, False, False, False, False, False]]

    condenses = [False, True, True, True, True]
    scroll_p_out = [550e3, 573e3, 573e3, 573e3, 573e3]
    
    N_radials = [10, 15, 20, 25, 30]
    PR_guesses = [1.833, 1.502, 1.358, 1.278, 1.2285]
    
    # Find the best number of radials in the determined optimal range
    """mass_flows = [0.0225]
    design_T_ins = [185.46]
    design_p_ins = [775.79]
    hours = [[True, True, True, False, False, False, False, True]]
    condenses = [True]
    scroll_p_out = [573e3]

    N_radials = [10, 11, 12, 13, 14, 15]
    PR_guesses = [1.833, 1.736, 1.659, 1.597, 1.540, 1.502]"""

    

    ESM_vector = [1, .149, .121, 100]

    for massflow_i, mass_flow in enumerate(mass_flows):
        design_T_in = design_T_ins[massflow_i]
        design_p_in = design_p_ins[massflow_i]

        for nstg_i, nstg in enumerate(N_radials):
            PR_guess = PR_guesses[nstg_i]

            all_machines = optimise_machine_design(
                mass_flow, design_T_in, design_p_in, scroll_p_out[massflow_i], PR_guess, nstg, ESM_vector, condenses[massflow_i], verbose=False)

            run_machine(all_machines, design_T_in, design_p_in, mass_flow, "Design", verbose=False)

            run_machines_day(all_machines, mass_flow, hours[massflow_i], design_T_in, design_p_in, verbose=False)
        
        print("")



if __name__ == "__main__":

    # Find the best compressor stack in the design family
    #optimise_design_system()


    # Best ESM when integrated with rest of system
    mass_flow = 0.0225
    design_T_in = 185.46
    design_p_in = 775.79
    hours = [True, True, True, False, False, False, False, True]
    condense = True
    scroll_p_out = 573e3

    # Best ESM when taken alone
    """mass_flow = 0.0113
    design_T_in = 199.9
    design_p_in = 770.60
    hours = [True, True, True, True, True, True, True, True]
    condense = False
    scroll_p_out = 573e3"""


    N_radial = 15
    PR_guess = 1.502

    ESM_vector = [1, 0.149, 0.121, 100]

    all_machines = optimise_machine_design(
                    mass_flow, design_T_in, design_p_in, scroll_p_out, PR_guess, N_radial, ESM_vector, condense, verbose=True)

    total_weight = 0
    for m in all_machines:
        if type(m) == component.RadialStage:
            print("chi_1", np.rad2deg(np.arctan(m.Vx_in/(m.R_mean_inlet*m.speed_rad))))
            print("chi_2", m.backsweep_angle)
            print("R_mean", m.R_mean_imp_exit - m.R_hub_inlet)
            print("span_1", m.bladeheight_in)
            print("span_2", m.bladeheight_imp_exit)
            print("spinner", m.R_hub_inlet)
            print("stator radius", m.R_stator_exit)
            print("stator angle", m.rotor_exit_flow_angle)
            print("RPM", m.speed)
            print("P in", m.gasflow_in.pressure)
            print("T in", m.gasflow_in.temperature)
            print("p out", m.gasflow_in.pressure * m.pressure_ratio_stage)
            print("")
        #print(m)
        #total_weight += m.weight
        
        pass

    #print("#####")
    #print(total_weight)
    # run_machines_day(all_machines, mass_flow, hours, design_T_in, design_p_in, verbose=False)
