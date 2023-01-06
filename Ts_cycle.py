import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as Axes
#from CoolProp.CoolProp import PropsSI
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib.figure import Figure



from co2_load import *

def plot_contours(ax, Trange, srange, hgrid, pgrid, dgrid):
    hcont = ax.contour(srange, Trange, hgrid, levels=np.linspace(300e3, 600e3, 31), linewidths=0.3, colors="orange")
    pcont = ax.contour(srange, Trange, pgrid, levels=[  1e2, 2.5e2, 5e2, 7.5e2, 
                                                        1e3, 2.5e3, 5e3, 7.5e3, 
                                                        1e4, 2.5e4, 5e4, 7.5e4, 
                                                        1e5, 2.5e5, 5e5, 7.5e5,
                                                        1e6, 2.5e6, 5e6, 7.5e6,
                                                        1e7, 2.5e7, 5e7, 7.5e7], linewidths=0.3, colors="b")
    
    
    
    ax.clabel(pcont, pcont.levels[::4], inline=False, fontsize=10)
    ax.clabel(hcont, hcont.levels[::5], inline=False, fontsize=10)
    #ax.contour(srange, Trange, dgrid, levels=[  1e-3, 2.5e-3, 5e-3, 7.5e-3, 
    #                                                        1e-2, 2.5e-2, 5e-2, 7.5e-2, 
    #                                                        1e-1, 2.5e-1, 5e-1, 7.5e-1, 
    ##                                                        1e0, 2.5e0, 5e0, 7.5e0,
    #                                                        1e1, 2.5e1, 5e1, 7.5e1], linewidths=0.3, colors="r")
    
    ax.contour(srange, Trange, pgrid, levels=[6e5], linewidths=1.5, colors="r")

def plot_sublimation(ax, offset=0, label=True):

    T_sub, s_sub = load_sublimation()
    T_triple = T_sub[-1]
    s_triple = s_sub[-1]

    if label:
        line_lab = "Sublimation line"
        triple_lab = "Triple point"
    else:
        line_lab = None
        triple_lab = None

    if offset == 0:
        ax.scatter(s_sub, T_sub, label=line_lab,s=10)
        ax.plot(s_sub, T_sub)
        ax.scatter([s_triple],[T_triple],label=triple_lab)
    else:
        T_sub = [T+offset for T in T_sub]
        ax.scatter(s_sub, T_sub, label=line_lab,s=10)
        ax.plot(s_sub, T_sub)

def plot_starting_conditions(ax, Trange, srange, pgrid):
    """bounding_pairs_Tp = [
    [254.51, 608.329],
    [169.73, 752.82],
    [185.341, 960.31],
    [242.753, 694.88]
    ] 

    bounding_pairs_Ts = []

    for pair in bounding_pairs_Tp:
        bounding_pairs_Ts.append(get_Tp_pair(Trange, srange, pgrid, pair))

    bounding_polygon = Polygon(bounding_pairs_Ts, closed=True, alpha=0.5)

    collection = PatchCollection([bounding_polygon], color="orange")
    ax.add_collection(collection)"""
    start_point = [185, 960]
    start_point_Ts = get_Tp_pair(Trange, srange, pgrid, start_point)
    plt.scatter(start_point_Ts[0], start_point_Ts[1], color="orange", label="Start point")

def plot_whole_graph(ax, Trange, srange, hgrid, pgrid, dgrid):
    plot_contours(ax, Trange, srange, hgrid, pgrid, dgrid)
    plot_sublimation(ax)
    #plot_sublimation(ax, 40, False)
    plot_starting_conditions(ax, Trange, srange, pgrid)

    ax.set_xlabel("Entropy, kJ/kg/K")
    ax.set_ylabel("Temperature, K")

    ax.set_ylim(140, 450)
    ax.set_xlim(1.5, 4.5)

    ax.set_title("Carbon dioxide Ts diagram at low temperatures")

def find_nearest_nd(arr, val):
    "Element in nd array `a` closest to the scalar value `a0`"
    arr = np.asarray(arr)
    idx = np.nanargmin(np.abs(arr - val))
    closeval = arr.flat[idx]
    return (np.where(arr==closeval))

def find_nearest_2_arrays(arr1, arr2, val1, val2, xgt=0, ygt=0, xlt=1e10, ylt=1e10):
    "Element in nd array `a` closest to the scalar value `a0`"
    arr1 = np.array(arr1)
    arr2 = np.array(arr2)

    for x in range(len(arr1)):
        for y in range(len(arr1[x])):
            if x < xgt:
                arr1[x, y] = 1e20
                arr2[x, y] = 1e20

    idx = np.nanargmin(np.sqrt(np.abs(arr1 - val1)**2 + np.abs(arr2 - val2)**2))
    closeval1 = arr1.flat[idx]
    return (np.where(arr1==closeval1))

def get_Tp_pair(Trange, srange, pgrid, Tp_pair, index=False):
    Ts_pairs = []

    Ts_pair = [Tp_pair[0]]

    T_index = find_nearest_nd(Trange, Tp_pair[0])[0]

    s_index = find_nearest_nd(pgrid[T_index,:][0], Tp_pair[1])[0]

    if index:
        return([s_index[0], T_index[0]])
    else:
        return([srange[s_index][0], Trange[T_index][0]])

def plot_ideal_process(ax, start_Tp, pressure_ratios, isen_efficiencies, intercool_temps, intercool_pressure_drops, 
Trange, srange, hgrid, pgrid, dgrid, label=None, linecolour="k", return_machine_params=False):

    in_Tp = start_Tp
    in_sT = get_Tp_pair(Trange, srange, pgrid, in_Tp)

    total_work = 0
    total_heat = 0

    machine_params = []
    

    for unit in range(len(pressure_ratios)):
        # Calculate for each unit (compressor + intercooler)
        # Plots for three stations:
        #   in (before compressor)
        #   out (after compressor, before intercooler)
        #   cool (after intercooler)
        
        pressure_ratio = pressure_ratios[unit]
        isen_efficiency = isen_efficiencies[unit]
        # Get the values of s and T that correspond to the starting T and p
        in_s_index = int(find_nearest_nd(srange, in_sT[0])[0])
        in_T_index = int(find_nearest_nd(Trange, in_sT[1])[0])
    
        in_sT_index = [in_s_index, in_T_index]

        out_p = in_Tp[1] * pressure_ratio

        

        # Find the temperature that closest corresponds to the value of end_p
        out_T_isen_index_2 = find_nearest_nd(pgrid[:,in_s_index], out_p)[0][0]
        out_T_isen = Trange[out_T_isen_index_2]


        # Get enthalpies for an efficiency calculation
        in_h = hgrid[in_T_index, in_s_index]
        out_h_isen = hgrid[out_T_isen_index_2][in_sT_index[0]]
        out_h_real = ((out_h_isen - in_h) / isen_efficiency) + in_h

        # Find the s and T for the output
        out_sT_index = find_nearest_2_arrays(hgrid, pgrid, out_h_real, out_p, xgt=out_T_isen*0)[::-1]
        out_sT_index = [int(i) for i in out_sT_index]
        out_sT = [srange[out_sT_index[0]], Trange[out_sT_index[1]]]

        work = out_h_real - in_h
        total_work += work

        #print(unit, in_sT[1], in_Tp[1]/1e5, dgrid[in_sT_index[1],in_sT_index[0]], dgrid[out_sT_index[1],out_sT_index[0]])

        # Draw the lines for the compressor
        if unit != 0:
            label=None
        ax.plot([in_sT[0], out_sT[0]], [in_sT[1], out_sT[1]], color=linecolour, label=label)
        ax.plot([in_sT[0], in_sT[0]], [in_sT[1], out_T_isen], color=linecolour, linewidth=0.5, linestyle="dashed")

        intercool_temp = intercool_temps[unit]

        machine_params.append(
            {
                "type":"c",
                "p in":in_Tp[1],
                "T in":in_Tp[0],
                "delta h": work,
                "pressure ratio": pressure_ratio
            }
        )

        # Find the p and T for the intercooler (if required)
        if out_sT[1] > (intercool_temp) and unit != len(pressure_ratios)-1:
            cool_p = out_p * (1 - intercool_pressure_drops[unit] - 0.01)
            cool_T = intercool_temp

            cool_sT = get_Tp_pair(Trange, srange, pgrid, [cool_T, cool_p])
            cool_sT_index = get_Tp_pair(Trange, srange, pgrid, [cool_T, cool_p], index=True)

            cool_h = in_h = hgrid[cool_sT_index[1], cool_sT_index[0]]
            
            heat = out_h_real - cool_h
            total_heat += heat

            ax.plot([out_sT[0], cool_sT[0]], [out_sT[1], cool_sT[1]], color=linecolour, label=label)

            in_Tp = [cool_T, cool_p]
            in_sT = cool_sT

            machine_params.append(
                {
                    "type":"h",
                    "p in":out_p,
                    "T in":out_sT[1],
                    "delta h":heat,
                }
            )
        
        else:
            in_Tp = [out_sT[1], out_p]
            in_sT = out_sT
        
        
        
    
    print("#####")
    print("Mean stage efficiency: ", np.average(isen_efficiencies))
    print("Final output pressure: ", in_Tp[1]/1e5, "bar")
    print("Total work done: ", total_work/1e3, "kJ/kg")
    print("Total heat extracted: ", total_heat/1e3, "kJ/kg")
    print("Output enthalpy", out_h_real)
    if out_p > 5.2e5:
        print("Ready to liquify")

    machine_params.insert(0, {
        "type": "cycle",
        "p out": in_Tp[1],
        "T out": in_Tp[0],
        "total work": total_work,
        "total cooling": total_heat
    })
    
    if return_machine_params:
        return(machine_params)


def main_load():

    Trange, srange, hgrid, pgrid, dgrid = get_grids()

    #fig = Figure(figsize=(5, 4), dpi=100)
    fig, ax = plt.subplots()

    start = [185.341, 900]

    print(fig, ax)

    plot_whole_graph(ax, Trange, srange, hgrid, pgrid, dgrid)


    # Icemelt utilsation cycle - T-XX-8
    #plot_ideal_process(ax, start, [2.1, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4], 
    #                            [0.89, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83], 
    #                            [210, 210, 210, 210, 220, 230, 240, 250], 
    #                           0.05, Trange, srange, hgrid, pgrid, dgrid, label="Real machines")

    # Radiator-driven cycle - T-XX-7
    #plot_ideal_process(ax, start, [2.1, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4], 
    #                            [0.89, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83], 
    #                            [250, 250, 250, 250, 250, 250, 250, 250], 
    #                           0.05, Trange, srange, hgrid, pgrid, dgrid, label="Real machines")

    # Ideal component cycle
    #plot_ideal_process(ax, start, [2.1, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4], 
    #                               [1, 1, 1, 1, 1, 1, 1, 1], 
    #                               [250, 250, 250, 250, 250, 250, 250, 250], 
    #                               0.0, Trange, srange, hgrid, pgrid, label="Ideal machines", linecolour="g")

    plt.legend()
    return(fig, ax, Trange, srange, hgrid, pgrid, dgrid)

def main_draw(Trange, srange, hgrid, pgrid, dgrid ):

    #fig = Figure(figsize=(5, 4), dpi=100)
    #ax = fig.add_subplot()
    fig, ax = plt.subplots()
    
    start = [185.341, 900]

    plot_whole_graph(ax, Trange, srange, hgrid, pgrid, dgrid)


    # Icemelt utilsation cycle - T-XX-8
    #plot_ideal_process(ax, start, [2.1, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4], 
    #                            [0.89, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83], 
    #                            [210, 210, 210, 210, 220, 230, 240, 250], 
    #                           0.05, Trange, srange, hgrid, pgrid, dgrid, label="Real machines")

    # Radiator-driven cycle - T-XX-7
    #plot_ideal_process(ax, start, [2.1, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4], 
    #                            [0.89, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83], 
    #                            [250, 250, 250, 250, 250, 250, 250, 250], 
    #                           0.05, Trange, srange, hgrid, pgrid, dgrid, label="Real machines")

    # Ideal component cycle
    #plot_ideal_process(ax, start, [2.1, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4], 
    #                               [1, 1, 1, 1, 1, 1, 1, 1], 
    #                               [250, 250, 250, 250, 250, 250, 250, 250], 
    #                               0.0, Trange, srange, hgrid, pgrid, label="Ideal machines", linecolour="g")

    plt.legend()
    return(fig, ax)

def generate_radiator_cycle(n_ICs, T_in, P_in, P_out, IC_out_temp, IC_pressure_drop_per_K):
    efficiency_guess = 0.7
    PR_total = P_out / P_in
    # Use a loop to set the desired compression ratio, taking into account the intercooler pressure drop
    for i in range(5):
        T_start = T_in
        P_start = P_in
        PR_block = np.power(PR_total, 1/(n_ICs+1))

        for block in range(n_ICs+1):
            P_end = PR_block * P_start
            delta_T = T_start * (PR_block**(0.31/1.31) - 1) / efficiency_guess
            T_end = T_start + delta_T

            if T_end > IC_out_temp:
                P_drop = IC_pressure_drop_per_K * abs(T_end - IC_out_temp)
                if P_drop > 0.99:
                    P_drop = 0.99
                P_start = P_end * (1 - P_drop - 0.01)
                T_start = IC_out_temp

            else:
                P_start = P_end
                T_start = T_end
        
        print(P_start, P_out)
        PR_total *= (P_out / P_start)
    
    print(PR_total, P_out / P_in, np.power(PR_total, 1/(n_ICs+1)))

    pressure_ratios = []
    efficiencies = []
    IC_temps = []
    IC_pressure_drops = []

    T_start = T_in
    P_start = P_in
    PR_block = np.power(PR_total, 1/(n_ICs+1))
    for block in range(n_ICs+1):
        pressure_ratios.append(PR_block)
        efficiencies.append(efficiency_guess)
        P_end = PR_block * P_start
        delta_T = T_start * (PR_block**(0.31/1.31) - 1) / efficiency_guess
        T_end = T_start + delta_T
        if T_end > IC_out_temp:
            P_drop = IC_pressure_drop_per_K * abs(T_end - IC_out_temp)
            if P_drop > 0.99:
                P_drop = 0.99
            P_start = P_end * (1 - P_drop)
            T_start = IC_out_temp

            IC_temps.append(IC_out_temp)
            IC_pressure_drops.append(P_drop)
        else:
            P_start = P_end
            T_start = T_end
            IC_temps.append(T_end+50)
            IC_pressure_drops.append(0)
    
    print(len(pressure_ratios), pressure_ratios[0]**len(pressure_ratios), PR_total)
    
    return(pressure_ratios, efficiencies, IC_temps, IC_pressure_drops)

def generate_icemelt_cycle(n_ICs, T_in, P_in, P_out, sublimation_line_margin, IC_pressure_drop_per_K):
    T_sub, p_sub = load_sublimation(Tp=True)

    interpolate_T_sub = interpolate.interp1d(p_sub, T_sub)

    efficiency_guess = 0.7
    PR_total = P_out / P_in
    # Use a loop to set the desired compression ratio, taking into account the intercooler pressure drop
    for i in range(5):
        T_start = T_in
        P_start = P_in
        PR_block = np.power(PR_total, 1/(n_ICs+1))
        for block in range(n_ICs+1):
            P_end = PR_block * P_start
            delta_T = T_start * (PR_block**(0.31/1.31) - 1) / efficiency_guess
            T_end = T_start + delta_T

            if P_end > max(p_sub):
                IC_out_temp = max(T_sub) + sublimation_line_margin
            else:
                IC_out_temp = interpolate_T_sub(P_end) + sublimation_line_margin

            if T_end > IC_out_temp:
                P_drop = IC_pressure_drop_per_K * (T_end - IC_out_temp)
                if P_drop > 0.99:
                    P_drop = 0.99
                P_start = P_end * (1 - P_drop)
                T_start = IC_out_temp
            else:
                P_start = P_end
                T_start = T_end
        
        print(P_start, P_out)
        PR_total *= (P_out / P_start)
    
    print(PR_total, P_out / P_in, np.power(PR_total, 1/(n_ICs+1)))

    pressure_ratios = []
    efficiencies = []
    IC_temps = []
    IC_pressure_drops = []

    T_start = T_in
    P_start = P_in
    PR_block = np.power(PR_total, 1/(n_ICs+1))
    for block in range(n_ICs+1):
        pressure_ratios.append(PR_block)
        efficiencies.append(efficiency_guess)
        P_end = PR_block * P_start
        delta_T = T_start * (PR_block**(0.31/1.31) - 1) / efficiency_guess
        T_end = T_start + delta_T

        if P_end > max(p_sub):
            IC_out_temp = max(T_sub) + sublimation_line_margin
        else:
            IC_out_temp = interpolate_T_sub(P_end) + sublimation_line_margin

        if T_end > IC_out_temp:
            P_drop = IC_pressure_drop_per_K * (T_end - IC_out_temp)
            if P_drop > 0.99:
                P_drop = 0.99
            P_start = P_end * (1 - P_drop)
            T_start = IC_out_temp

            IC_temps.append(IC_out_temp)
            IC_pressure_drops.append(P_drop)
        else:
            P_start = P_end
            T_start = T_end
            IC_temps.append(T_end+50)
            IC_pressure_drops.append(0)
    
    return(pressure_ratios, efficiencies, IC_temps, IC_pressure_drops)

#pressure_ratios, efficiencies, IC_temps, IC_pressure_drops = generate_radiator_cycle(10, 190, 900, 600e3, 273, 5e-4)
#pressure_ratios_i, efficiencies_i, IC_temps_i, IC_pressure_drops_i = generate_icemelt_cycle(5, 190, 900, 600e3, 40, 5e-4)
#pressure_ratios_i, efficiencies_i, IC_temps_i, IC_pressure_drops_i = generate_radiator_cycle(10, 190, 900, 600e3, 250, 10e-4, 1)
#pressure_ratios_j, efficiencies_j, IC_temps_j, IC_pressure_drops_j = generate_radiator_cycle(3, 190, 900, 600e3, 250, 10e-4, 1)

#fig, ax, Trange, srange, hgrid, pgrid, dgrid = main_load()

#plot_ideal_process(ax, [190, 900], pressure_ratios, efficiencies, IC_temps, IC_pressure_drops, Trange, srange, hgrid, pgrid, dgrid, label="1")
#plot_ideal_process(ax, [190, 900], pressure_ratios_i, efficiencies_i, IC_temps_i, IC_pressure_drops_i, Trange, srange, hgrid, pgrid, dgrid, label="1", linecolour="b")
#plot_ideal_process(ax, [190, 900], pressure_ratios_j, efficiencies_j, IC_temps_j, IC_pressure_drops_j, Trange, srange, hgrid, pgrid, dgrid, label="1", linecolour="r")

#plt.show()
"""
fig, ax, Trange, srange, hgrid, pgrid, dgrid = main_load()
plt.show()"""