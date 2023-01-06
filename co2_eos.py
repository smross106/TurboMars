"""
Combines CoolProp with Martynov et al to produce a full EOS for CO2, including below the triple point
https://onlinelibrary.wiley.com/doi/full/10.1002/ghg.1322
"""

import numpy as np
import CoolProp.CoolProp as CP
import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
import matplotlib.pyplot as plt
import matplotlib.axes as Axes
from scipy import ndimage
from scipy import interpolate
import os
import csv




def p(T, v):
    """
    Returns CO2 pressure

    p (float):  Pressure (Pa)
    v (float):  Specific volume (m3/kg)

    NOT WORKING - RETURNS BAD VALUES. NOT SURE WHY
    TRIED USING MOLAR DENSITIES, NO LUCK
    """
    R = 8.3144621
    
    a_0 = 1143e-3
    T_a = 343.55

    a_sv = a_0 * (1 - (T / T_a))

    b_0 = 24.45e-6
    c_b = 3.5e-4
    T_b = 39

    b_sv = b_0 * (1 - c_b * np.exp(T / T_b))

    p = (R * T / (v - b_sv)) - (a_sv / ((v * (v + b_sv)) + (b_sv * (v - b_sv))))

    return(p)

"""plot = PropertyPlot('CarbonDioxide', 'Ts', unit_system="SI", tp_limits="ACHP")
plot.calc_isolines(CoolProp.iHmass)
plot.calc_isolines(CoolProp.iP)
plot.calc_isolines(CoolProp.iQ)
plot.show()"""

#print(PropsSI("D", "T", 216.6, "P", 750, "CO2"))
"""gas_densities = []
enthalpy = []
pressures = np.logspace(1.6, 5.5, 100)

print(PropsSI("co2", "ptriple"))
ptriple = PropsSI("co2", "ptriple")
print(PropsSI("co2", "Ttriple"))
Ttriple = PropsSI("co2", "Ttriple")
print(ptriple, Ttriple)
htriple_cp = PropsSI("H", "T|gas", Ttriple, "P", ptriple, "CO2")
h_adjust = 729e3 - htriple_cp

for pres in pressures:
    gas_densities.append(PropsSI("D", "T|gas", 158.15, "P", pres, "CO2"))
    enthalpy.append(PropsSI("H", "T|gas", 158.15, "P", pres, "CO2"))



ax2 = ax1.twinx()
ax1.plot(pressures, gas_densities, 'g-')
ax1.scatter([0.02514e5],[0.08416], color="g")
ax2.plot(pressures, enthalpy, 'b-')
ax2.scatter([0.02514e5],[702.3e3-h_adjust], color="b")
ax2.spines['left'].set_color('g')
plt.show()
"""

def direct_isolines(ax):
    prange = []
    plabel = []
    for exponent in range(2, 7):
        prange.append(1 * (10**exponent))
        prange.append(2.5 * (10**exponent))
        prange.append(5 * (10**exponent))
        prange.append(7.5 * (10**exponent))
        if exponent == 2:
            factor = 1
            text = "Pa"
        elif 3 <= exponent and exponent < 6:
            factor = 1e3
            text = "kPa"
        
        elif 6 <= exponent and exponent < 9:
            factor = 1e6
            text = "MPa"
        
        plabel.append("{:.0f}".format(1 * (10**exponent)/factor)+text)
        plabel.append("{:.0f}".format(2.5 * (10**exponent)/factor)+text)
        plabel.append("{:.0f}".format(5 * (10**exponent)/factor)+text)
        plabel.append("{:.0f}".format(7.5 * (10**exponent)/factor)+text)

    hrange = np.linspace(300e3, 600e3, 61)

    Tgrid = np.zeros((len(prange), len(hrange)))
    sgrid = np.zeros((len(prange), len(hrange)))

    for p in range(len(prange)):
        for h in range(len(hrange)):
            try:
                Tgrid[p][h] = PropsSI("T", "P", prange[p], "Hmolar", hrange[h]/1e3*44, "CO2")
            except:
                Tgrid[p][h] = np.nan
            
            try:
                sgrid[p][h] = PropsSI("Smass", "P", prange[p], "Hmolar", hrange[h]/1e3*44, "CO2")/1e3
            except:
                sgrid[p][h] = np.nan

    for p in range(len(prange)):
        plt.plot(sgrid[p], Tgrid[p], color="b", linewidth=0.3)

        angle = 70#np.rad2deg(np.arctan((Tgrid[p][-2] - Tgrid[p][-1]) / (sgrid[p][-2] - sgrid[p][-1])))
        ax.annotate(plabel[p], xy=(sgrid[p][-1],Tgrid[p][-1]), xytext=(sgrid[p][-1],Tgrid[p][-1]),arrowprops=None, va="bottom", rotation=angle)

    for h in range(len(hrange)):
        ax.plot(sgrid[:,h], Tgrid[:,h], color="g", linewidth=0.3)
        angle=0
        ax.annotate("{:.0f} kJ/kg".format(hrange[h]/1e3), xy=(sgrid[:,h][0],Tgrid[:,h][0]), xytext=(sgrid[:,h][0],Tgrid[:,h][0]),arrowprops=None, va="center", rotation=angle)

def indirect_isolines(ax):
    Trange = np.linspace(140, 450, 621)
    srange = np.linspace(1.5, 4.5, 451)   
    #Trange = np.linspace(140, 450, 32)
    #srange = np.linspace(1.5, 4.5, 31)

    hgrid = np.zeros((len(Trange), len(srange)))
    pgrid = np.zeros((len(Trange), len(srange)))
    Dgrid = np.zeros((len(Trange), len(srange)))

    Ttriple = PropsSI("co2", "Ttriple")

    for T_index in range(len(Trange)):
        if T_index%10==0:print(Trange[T_index])
        for s_index in range(len(srange)):
            hgrid[T_index][s_index] = PropsSI("Hmass", "T", Trange[T_index], "Smass", srange[s_index]*1e3, "CO2")
            pgrid[T_index][s_index] = PropsSI("P", "T", Trange[T_index], "Smass", srange[s_index]*1e3, "CO2")
            Dgrid[T_index][s_index] = PropsSI("Dmass", "T", Trange[T_index], "Smass", srange[s_index]*1e3, "CO2")

            """if CP.PhaseSI("T", Trange[T_index], "Smass", srange[s_index]*1e3, "CO2") in ["liquid", "twophase"] and Trange[T_index] > Ttriple:
                Qgrid[T_index][s_index] = PropsSI("Q", "T", Trange[T_index], "Smass", srange[s_index]*1e3, "CO2")"""

    print("phase generated")

    #hgrid = scipy.ndimage.zoom(hgrid, 3)
    #pgrid = scipy.ndimage.zoom(pgrid, 3)
    #Qgrid = scipy.ndimage.zoom(Qgrid, 3)

    ax.contour(srange, Trange, hgrid, levels=np.linspace(300e3, 600e3, 61), colors="g", linewidths=0.3)
    ax.contour(srange, Trange, pgrid, levels=[  1e2, 2.5e2, 5e2, 7.5e2, 
                                                1e3, 2.5e3, 5e3, 7.5e3, 
                                                1e4, 2.5e4, 5e4, 7.5e4, 
                                                1e5, 2.5e5, 5e5, 7.5e5,
                                                1e6, 2.5e6, 5e6, 7.5e6], colors="b", linewidths=0.3)
    #ax.contour(srange, Trange, Qgrid, levels=[0, 0.2, 0.4, 0.6, 0.8, 1.0], colors="k", linewidths=0.3)

    return(srange, Trange, hgrid, pgrid, Dgrid)

def ps_interpolate(Trange, prange, sgrid, hgrid, Dgrid, s_low, s_high):
    srange = np.linspace(s_low, s_high, 50)

    pgrid = np.zeros((len(Trange),len(srange)))
    hgrid2 = np.zeros((len(Trange),len(srange)))
    Dgrid2 = np.zeros((len(Trange),len(srange)))

    ptriple = PropsSI("co2", "ptriple")
    Ttriple = PropsSI("co2", "Ttriple")
    striple_cp = PropsSI("Smass", "T|gas", Ttriple, "P", ptriple, "CO2")/1e3

    for T_index, T in enumerate(Trange):
        s_T = sgrid[T_index]
        func_p = interpolate.interp1d(np.flip(s_T), np.flip(prange), kind="cubic", fill_value="extrapolate")
        for s_index, s in enumerate(srange):
            if s >= min(s_T) and s <= max(s_T):
                pgrid[T_index][s_index] = func_p(s)
            else:
                pgrid[T_index][s_index] = np.nan

    for T_index, T in enumerate(Trange):
        s_T = sgrid[T_index]
        h_T = hgrid[T_index]
        D_T = Dgrid[T_index]
        func_h = interpolate.interp1d(np.flip(s_T), np.flip(h_T), kind="cubic", fill_value="extrapolate")
        func_D = interpolate.interp1d(np.flip(s_T), np.flip(D_T), kind="cubic", fill_value="extrapolate")

        for s_index, s in enumerate(srange):
            if s >= min(s_T) and s <= max(s_T):
                hgrid2[T_index][s_index] = func_h(s)
                Dgrid2[T_index][s_index] = func_D(s)

            else:
                hgrid2[T_index][s_index] = np.nan
                Dgrid2[T_index][s_index] = np.nan
        
    return(srange, pgrid, hgrid2, Dgrid2)

def indirect_isolines_Tp(ax):
    Trange = np.linspace(140, 300, 160)
    prange = np.logspace(1.5, 5.5, 128)  
    #Trange=np.linspace(140, 300, 16)
    #prange = np.logspace(1.5, 5.5, 16)

    hgrid = np.zeros((len(Trange), len(prange)))
    sgrid = np.zeros((len(Trange), len(prange)))
    Qgrid = np.zeros((len(Trange), len(prange)))
    Dgrid = np.zeros((len(Trange), len(prange)))

    Ttriple = PropsSI("co2", "Ttriple")

    for T_index in range(len(Trange)):
        for p_index in range(len(prange)):
            #hgrid[T_index][s_index] = PropsSI("Hmass", "T", Trange[T_index], "Smass", srange[s_index]*1e3, "CO2")
            try:
                if Trange[T_index] >= Ttriple:
                    sgrid[T_index][p_index] = PropsSI("Smass", "T", Trange[T_index], "P", prange[p_index], "CO2")/1e3
                    hgrid[T_index][p_index] = PropsSI("Hmass", "T", Trange[T_index], "P", prange[p_index], "CO2")/1e3
                    Dgrid[T_index][p_index] = PropsSI("Dmass", "T", Trange[T_index], "P", prange[p_index], "CO2")
                else:
                    sgrid[T_index][p_index] = PropsSI("Smass", "T|gas", Trange[T_index], "P", prange[p_index], "CO2")/1e3
                    hgrid[T_index][p_index] = PropsSI("Hmass", "T|gas", Trange[T_index], "P", prange[p_index], "CO2")/1e3
                    Dgrid[T_index][p_index] = PropsSI("Dmass", "T|gas", Trange[T_index], "P", prange[p_index], "CO2")
                
                if CP.PhaseSI( "T|gas", Trange[T_index], "P", prange[p_index], "CO2") in ["liquid", "twophase"] and Trange[T_index] > Ttriple:
                    Qgrid[T_index][p_index] = PropsSI("Q", "T", Trange[T_index], "Smass", prange[p_index]*1e3, "CO2")
                else:
                    Qgrid[T_index][p_index] = np.nan
            except:
                pass#print( Trange[T_index],  prange[p_index])

    print("phase generated")
    #print(Dgrid)

    srange, pgrid, hgrid2, Dgrid2 = ps_interpolate(Trange, prange, sgrid, hgrid, Dgrid, 1.5, 4.5)

    ax.contour(srange, Trange, pgrid, levels=[  1e2, 2.5e2, 5e2, 7.5e2, 
                                                1e3, 2.5e3, 5e3, 7.5e3, 
                                                1e4, 2.5e4, 5e4, 7.5e4, 
                                                1e5, 2.5e5, 5e5, 7.5e5,
                                                1e6, 2.5e6, 5e6, 7.5e6], colors="b", linewidths=0.3, linestyles="dashed")
    ax.contour(srange, Trange, hgrid2, levels=np.linspace(300, 500, 41), colors="g", linewidths=0.3, linestyles="dashed")

    #print(Dgrid2)
    Dcont = ax.contour(srange, Trange, Dgrid2, levels=[ 1e-2, 5e-2, 1e-1, 5e-1, 1e0, 5e0, 1e1, 5e1], colors="y", linewidths=1, linestyles="dashed")
    ax.clabel(Dcont, Dcont.levels, inline=True)

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
fig, ax = plt.subplots()

#indirect_isolines_Tp(ax)
"""srange, Trange, hgrid, pgrid, Dgrid = indirect_isolines(ax)

#print(hgrid)
#direct_isolines(ax)




with open(os.path.join(__location__, 'Ts_diagram.csv'), 'w',newline='') as f:
    writer = csv.writer(f)

    writer.writerow(["Carbon dioxide property chart"])
    writer.writerow(["Evaluated using CoolProp"])
    writer.writerow(["Rows of T[K], columns of entropy[kJ/kg/K]"])
    writer.writerow(["###"])
    
    writer.writerow(["Enthalpy [J/kg]"])
    writer.writerow(np.concatenate((["---"], srange)))
    for row_index, row in enumerate(hgrid):
        writer.writerow(np.concatenate(([Trange[row_index]], row)))
    writer.writerow(["###"])

    writer.writerow(["Pressure [Pa]"])
    writer.writerow(np.concatenate((["---"], srange)))
    for row_index, row in enumerate(pgrid):
        writer.writerow(np.concatenate(([Trange[row_index]], row)))
    writer.writerow(["###"])

    writer.writerow(["Density [kg/m3]"])
    writer.writerow(np.concatenate((["---"], srange)))
    for row_index, row in enumerate(Dgrid):
        writer.writerow(np.concatenate(([Trange[row_index]], row)))
    writer.writerow(["###"])

print("")


"""

Trange_out = []
srange_out = []
hgrid_out = []
pgrid_out = []
dgrid_out = []

with open(os.path.join(__location__, 'Ts_diagram.csv')) as csvfile:
    reader = csv.reader(csvfile)

    csv_rows = []
    for row in reader:
        csv_rows.append(row)
    
    block_breaks = [i for i in range(len(csv_rows)) if csv_rows[i][0] == '###']
    print(csv_rows[3][0]=="###", block_breaks)


    # Load in hgrid, the first one, along with temperature
    h_block = csv_rows[block_breaks[0]+2:block_breaks[1]]
    srange_out = [float(i) for i in h_block[0][1:]]
    Trange_out = [float(i[0]) for i in h_block[1:]]
    hgrid_out = [[float(j) for j in i[1:]] for i in h_block[1:]]

    # Load in pgrid and dgrid
    p_block = csv_rows[block_breaks[1]+2:block_breaks[2]]
    pgrid_out = [[float(j) for j in i[1:]] for i in p_block[1:]]

    d_block = csv_rows[block_breaks[2]+2:block_breaks[3]]
    dgrid_out = [[float(j) for j in i[1:]] for i in d_block[1:]]


ax.contour(srange_out, Trange_out, hgrid_out, levels=np.linspace(300e3, 600e3, 31))
#ax.contour(srange_out, Trange_out, pgrid_out, levels=[  1e2, 2.5e2, 5e2, 7.5e2, 
#                                                        1e3, 2.5e3, 5e3, 7.5e3, 
#                                                        1e4, 2.5e4, 5e4, 7.5e4, 
#                                                        1e5, 2.5e5, 5e5, 7.5e5,
#                                                        1e6, 2.5e6, 5e6, 7.5e6])

ax.contour(srange_out, Trange_out, dgrid_out, levels=[  1e-3, 2.5e-3, 5e-3, 7.5e-3, 
                                                        1e-2, 2.5e-2, 5e-2, 7.5e-2, 
                                                        1e-1, 2.5e-1, 5e-1, 7.5e-1, 
                                                        1e0, 2.5e0, 5e0, 7.5e0,
                                                        1e1, 2.5e1, 5e1, 7.5e1])


plt.show()


ptriple = PropsSI("co2", "ptriple")
Ttriple = PropsSI("co2", "Ttriple")
striple_cp = PropsSI("Smass", "T|gas", Ttriple, "P", ptriple, "CO2")/1e3
s_adjust = 4.25 - striple_cp

# Vukalovich data
T_vuk = [143.15, 148.15, 153.15, 158.15, 163.15, 168.15, 173.15, 178.15, 183.15, 188.15, 193.15, 198.15, 203.15, 208.15, 213.15, 216.55]
s_vuk = [5.414, 5.293, 5.18, 5.078, 4.981, 4.892, 4.809, 4.73, 4.658, 4.59, 4.525, 4.462, 4.401, 4.341, 4.281, 4.25]
s_vuk = [i - s_adjust for i in s_vuk]

"""ax.scatter(s_vuk, T_vuk, label="Sublimation line",s=10)
ax.plot(s_vuk, T_vuk)
ax.scatter([striple_cp],[Ttriple],label="Triple point")

ax.set_xlabel("Entropy, kJ/kg/K")
ax.set_ylabel("Temperature, K")

ax.set_ylim(140, 450)
ax.set_xlim(1.5, 4.5)

ax.set_title("Carbon dioxide Ts diagram at low temperatures")

plt.legend()
plt.show()"""

