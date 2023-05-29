import numpy as np
import os
import csv
#from CoolProp.CoolProp import PropsSI
from scipy import interpolate


def load_data(filename):
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))

    Trange_out = []
    srange_out = []
    hgrid_out = []
    pgrid_out = []
    dgrid_out = []

    with open(os.path.join(__location__, filename)) as csvfile:
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
    
    Trange_out = np.asarray(Trange_out)
    srange_out = np.asarray(srange_out)
    hgrid_out = np.asarray(hgrid_out)
    pgrid_out = np.asarray(pgrid_out)
    dgrid_out = np.asarray(dgrid_out)

    return(Trange_out, srange_out, hgrid_out, pgrid_out, dgrid_out)


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
    ax.contour(srange, Trange, dgrid, levels=[  1e-3, 2.5e-3, 5e-3, 7.5e-3, 
                                                            1e-2, 2.5e-2, 5e-2, 7.5e-2, 
                                                            1e-1, 2.5e-1, 5e-1, 7.5e-1, 
                                                            1e0, 2.5e0, 5e0, 7.5e0,
                                                            1e1, 2.5e1, 5e1, 7.5e1], linewidths=0.3, colors="r")

def load_sublimation(Tp=False):
    """ptriple = PropsSI("co2", "ptriple")
    Ttriple = PropsSI("co2", "Ttriple")
    striple_cp = PropsSI("Smass", "T", Ttriple, "P", ptriple, "CO2")/1e3
    s_adjust = 4.25 - striple_cp

    # Vukalovich data
    T_vuk = [143.15, 148.15, 153.15, 158.15, 163.15, 168.15, 173.15, 178.15, 183.15, 188.15, 193.15, 198.15, 203.15, 208.15, 213.15, 216.55]
    s_vuk = [5.414, 5.293, 5.18, 5.078, 4.981, 4.892, 4.809, 4.73, 4.658, 4.59, 4.525, 4.462, 4.401, 4.341, 4.281, 4.25]
    s_vuk = [i - s_adjust for i in s_vuk]

    p_vuk = []
    for i in range(len(T_vuk)):
        p_vuk.append(PropsSI("P", "T|gas", T_vuk[i], "Smass", s_vuk[i]*1e3, "CO2"))"""

    p_co2 = [2695, 4727, 7934, 12805, 19948, 30106, 44156, 63111, 88119, 81317098, 0, 0, 25069347, 71847802, 128805935, 168665455]
    s_co2 = [3.283, 3.165, 3.056, 2.955, 2.861, 2.774, 2.693, 2.616, 2.545, 2.477, 2.413, 2.352, 2.292, 2.234, 2.177, 2.138]
    T_co2 = [143.15, 148.15, 153.15, 158.15, 163.15, 168.15, 173.15, 178.15, 183.15, 188.15, 193.15, 198.15, 203.15, 208.15, 213.15, 216.55]

    if Tp:
        return(T_co2, p_co2)
    else:
        return(T_co2, s_co2)

def cleanup_grids_sublimation(Trange, srange, hgrid, pgrid, dgrid):
    # Take in a grid (dgrid, hgrid, pgrid etc) and remove all points below the sublimation curve
    T_sub, s_sub = load_sublimation()
    T_triple = T_sub[-1]
    s_triple = s_sub[-1]
    T_sub_curve = interpolate.interp1d(s_sub, T_sub)
    
    for T_index, T in enumerate(Trange):
        for s_index, s in enumerate(srange):
            if s < min(s_sub):
                # Point is to the left of the triple point
                if T < T_triple:
                    # Point is below the triple point
                    hgrid[T_index, s_index] = np.nan
                    pgrid[T_index, s_index] = np.nan
                    dgrid[T_index, s_index] = np.nan
            
            else:
                if s > min(s_sub) and s < max(s_sub):
                    # Point lies under the sublimation curve
                    if T < T_sub_curve(s):
                        hgrid[T_index, s_index] = np.nan
                        pgrid[T_index, s_index] = np.nan
                        dgrid[T_index, s_index] = np.nan
    
    return(hgrid, pgrid, dgrid)


def get_grids():
    Trange, srange, hgrid, pgrid, dgrid = load_data("Ts_diagram_repairs.csv")

    hgrid, pgrid, dgrid = cleanup_grids_sublimation(Trange, srange, hgrid, pgrid, dgrid)
    
    return(Trange, srange, hgrid, pgrid, dgrid)


# Temporary
"""def plot_contours(ax, Trange, srange, hgrid, pgrid, dgrid):
    hcont = ax.contour(srange, Trange, hgrid, levels=np.linspace(300e3, 600e3, 31), linewidths=0.3, colors="orange")
    pcont = ax.contour(srange, Trange, pgrid, levels=[  1e2, 2.5e2, 5e2, 7.5e2, 
                                                        1e3, 2.5e3, 5e3, 7.5e3, 
                                                        1e4, 2.5e4, 5e4, 7.5e4, 
                                                        1e5, 2.5e5, 5e5, 7.5e5,
                                                        1e6, 2.5e6, 5e6, 7.5e6,
                                                        1e7, 2.5e7, 5e7, 7.5e7], linewidths=0.3, colors="b")
    
    ax.clabel(pcont, pcont.levels[::4], inline=False, fontsize=10)
    ax.clabel(hcont, hcont.levels[::5], inline=False, fontsize=10)
    ax.contour(srange, Trange, dgrid, levels=[  1e-3, 2.5e-3, 5e-3, 7.5e-3, 
                                                            1e-2, 2.5e-2, 5e-2, 7.5e-2, 
                                                            1e-1, 2.5e-1, 5e-1, 7.5e-1, 
                                                            1e0, 2.5e0, 5e0, 7.5e0,
                                                            1e1, 2.5e1, 5e1, 7.5e1], linewidths=0.3, colors="r")"""