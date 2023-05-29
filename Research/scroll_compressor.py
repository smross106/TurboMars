import matplotlib.pyplot as plt
import xlwings as xw
from scipy import interpolate
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import os



#dataframe1 = pd.read_excel("Scroll compressor performance.xlsx")
ws = xw.Book(os.path.join(os.path.dirname(os.path.realpath(__file__)), "Scroll compressor performance.xlsx"))

AMBIENT_PRESSURE = 0.835

def interpolate_variables_vacuum(datadict, low_PR, high_PR):
    # Interpolate data over pressure range
    Q_func = interpolate.interp1d(datadict["PR for Q"], datadict["Q"])
    W_func = interpolate.interp1d(datadict["PR for W"], datadict["W"])
    PR_range = np.linspace(low_PR+1, high_PR+1, int(high_PR-low_PR))
    p_range = AMBIENT_PRESSURE / PR_range
    datadict["PR interp"] = PR_range
    datadict["P interp"] = p_range
    datadict["Q interp"] = Q_func(PR_range)
    datadict["W interp"] = W_func(PR_range)

def interpolate_variables_compressor(datadict, low_PR, high_PR):
    # Interpolate data over pressure range
    Q_func = interpolate.interp1d(datadict["Outlet P for Q"], datadict["Q"])
    W_func = interpolate.interp1d(datadict["Outlet P for W"], datadict["W"])
    PR_range = np.linspace(low_PR+1, high_PR+1, int(50*(high_PR-low_PR)))
    p_range = PR_range 
    datadict["PR interp"] = PR_range
    datadict["P interp"] = p_range
    datadict["Q interp"] = Q_func(p_range)
    datadict["W interp"] = W_func(p_range) / AMBIENT_PRESSURE

def massflow_vacuum(datadict):
    T = 293.15
    R = 287.058
    datadict["rho interp"] = datadict["P interp"] * 101.325e3 / (T * R)
    datadict["mass flow interp"] = datadict["rho interp"] * datadict["Q interp"] / (1000 * 60)

def massflow_compressor(datadict):
    T = 293.15
    R = 287.058
    rho_ambient = 101.325e3 / (T * R)
    datadict["rho interp"] = np.full(len(datadict["Q interp"]), rho_ambient)
    datadict["mass flow interp"] = rho_ambient * datadict["Q interp"] / (1000 * 60)

def power_estimate_adiabatic(datadict):
    T = 293.15
    Cp = 1005
    gamma = 1.4
    datadict["W adia interp"] = (np.power(datadict["PR interp"], (gamma-1)/gamma) - 1) * T * Cp * datadict["mass flow interp"]

def power_estimate_isothermal(datadict):
    T = 293.15
    R = 287.058
    datadict["W isot interp"] = np.log(datadict["PR interp"]) * T * R * datadict["mass flow interp"]

def polytropic_efficiency(datadict, adiabatic=True, gamma=1.4):
    if adiabatic:
        eta_isen = datadict["W adia interp"] / datadict["W interp"]
    else:
        eta_isen = datadict["W isot interp"] / datadict["W interp"]
    eta_poly = np.log(datadict["PR interp"]) / np.log(
        ((np.power(datadict["PR interp"], (gamma-1)/gamma) - 1)/eta_isen)+1
    )
    eta_poly *= (gamma-1) / gamma
    datadict["eta poly interp"] = eta_poly

def polytropic_efficiency_list(pressure_range, eta_isen, gamma=1.4):

    eta_poly = np.log(np.array(pressure_range)) / np.log(
        ((np.power(np.array(pressure_range), (gamma-1)/gamma) - 1)/np.array(eta_isen))+1
    )
    eta_poly *= (gamma-1) / gamma
    return(eta_poly)

def load_vacuum(sheet_name, low_PR, high_PR):
    vacsheet = ws.sheets[sheet_name]
    data_3000rpm = {
        "RPM": 3000,
        "Inlet P for Q": [x for x in vacsheet.range("A4:A100").value if x is not None],
        "PR for Q": AMBIENT_PRESSURE / np.array([x for x in vacsheet.range("A4:A100").value if x is not None],),
        "Q": [x for x in vacsheet.range("C4:C100").value if x is not None],
        "Inlet P for W": [x for x in vacsheet.range("D4:D100").value if x is not None],
        "PR for W": AMBIENT_PRESSURE / np.array([x for x in vacsheet.range("D4:D100").value if x is not None],),
        "W": [x for x in vacsheet.range("F4:F100").value if x is not None],
    }
    data_2000rpm = {
        "RPM": 2000,
        "Inlet P for Q": [x for x in vacsheet.range("G4:G100").value if x is not None],
        "PR for Q": AMBIENT_PRESSURE / np.array([x for x in vacsheet.range("G4:G100").value if x is not None]),
        "Q": [x for x in vacsheet.range("I4:I100").value if x is not None],
        "Inlet P for W": [x for x in vacsheet.range("J4:J100").value if x is not None],
        "PR for W": AMBIENT_PRESSURE / np.array([x for x in vacsheet.range("J4:J100").value if x is not None]),
        "W": [x for x in vacsheet.range("L4:L100").value if x is not None]
    }
    
    data = [data_3000rpm, data_2000rpm]
    for dic in data:
        interpolate_variables_vacuum(dic, low_PR, high_PR)
        massflow_vacuum(dic)
        power_estimate_adiabatic(dic)
        power_estimate_isothermal(dic)
        polytropic_efficiency(dic, False)

    return([data_3000rpm, data_2000rpm])

def load_compressor(sheet_name, low_PR, high_PR, rpm_3500=True, skip_3000=False):
    vacsheet = ws.sheets[sheet_name]
    if rpm_3500:
        data_3500rpm = {
            "RPM": 3500,
            "Outlet P for Q": [x+1 for x in vacsheet.range("A4:A100").value if x is not None],
            "Q": [x for x in vacsheet.range("B4:B100").value if x is not None],
            "Outlet P for W": [x+1 for x in vacsheet.range("C4:C100").value if x is not None],
            "W": [x/AMBIENT_PRESSURE for x in vacsheet.range("D4:D100").value if x is not None],
        }
        data_3000rpm = {
            "RPM": 3000,
            "Outlet P for Q": [x+1 for x in vacsheet.range("E4:E100").value if x is not None],
            "Q": [x for x in vacsheet.range("F4:F100").value if x is not None],
            "Outlet P for W": [x+1 for x in vacsheet.range("G4:G100").value if x is not None],
            "W": [x/AMBIENT_PRESSURE for x in vacsheet.range("H4:H100").value if x is not None],
        }
        data_2000rpm = {
            "RPM": 2000,
            "Outlet P for Q": [x+1 for x in vacsheet.range("I4:I100").value if x is not None],
            "Q": [x for x in vacsheet.range("J4:J100").value if x is not None],
            "Outlet P for W": [x+1 for x in vacsheet.range("K4:K100").value if x is not None],
            "W": [x/AMBIENT_PRESSURE for x in vacsheet.range("L4:L100").value if x is not None],
        }

        data = [data_3500rpm, data_3000rpm, data_2000rpm]
    else:
        if not skip_3000:
            data_3000rpm = {
                "RPM": 3000,
                "Outlet P for Q": [x+1 for x in vacsheet.range("A4:A100").value if x is not None],
                "Q": [x for x in vacsheet.range("B4:B100").value if x is not None],
                "Outlet P for W": [x+1 for x in vacsheet.range("C4:C100").value if x is not None],
                "W": [x/AMBIENT_PRESSURE for x in vacsheet.range("D4:D100").value if x is not None],
            }
        
        data_2000rpm = {
            "RPM": 2000,
            "Outlet P for Q": [x+1 for x in vacsheet.range("E4:E100").value if x is not None],
            "Q": [x for x in vacsheet.range("F4:F100").value if x is not None],
            "Outlet P for W": [x+1 for x in vacsheet.range("G4:G100").value if x is not None],
            "W": [x/AMBIENT_PRESSURE for x in vacsheet.range("H4:H100").value if x is not None],
        }

        if not skip_3000:
            data = [data_3000rpm, data_2000rpm]
        else:
            data = [data_2000rpm]

    for dic in data:
        interpolate_variables_compressor(dic, low_PR, high_PR)
        massflow_compressor(dic)
        power_estimate_adiabatic(dic)
        power_estimate_isothermal(dic)
        polytropic_efficiency(dic)
    return(data)

vacuum_sheets = ["ASq V16H030A"]
vacuum_dicts = []

compr_dicts = []

vacuum_dicts = vacuum_dicts + load_vacuum("ASq V16H030A", 5, 350)
vacuum_dicts = vacuum_dicts + load_vacuum("ASq V09H017A", 5, 140)
compr_dicts = compr_dicts + load_compressor("ASq P05H012A", 0.7, 2, False)
compr_dicts = compr_dicts + load_compressor("ASq P09H017A", 1.2, 3.4, True)
compr_dicts = compr_dicts + load_compressor("ASq P14H022A", 2.3, 6.4, True)
compr_dicts = compr_dicts + load_compressor("ASq P17H043B", 4.8, 12.7, True)
compr_dicts = compr_dicts + load_compressor("ASq P34H080A", 0.62, 1.5, False, True)


# Data sources from https://hal.science/hal-00926697/PDF/Scroll_compressor_modelling_for_hydrocarbons.pdf
MOXIE_pr = [108.57]
MOXIE_eta_i = 0.096061

sun_pr = [3.0, 3.5, 4.0, 4.5, 5.0]
sun_eta_i = [0.75781, 0.77990, 0.79362, 0.79721, 0.79470]
liu_2009_pr = [2.844, 2.423, 3.443, 3.826, 2.458]
liu_2009_eta_i = [0.59, 0.64, 0.57, 0.56, 0.64]
liu_2010_pr = [3.434, 3.030, 3.305, 3.946, 4.759, 3.434, 3.434, 3.434]
liu_2010_eta_i = [0.754, 0.717, 0.723, 0.706, 0.672, 0.649, 0.672, 0.761]
cho_pr = [3.457]
cho_eta_i = [0.7253]


fig, axs = plt.subplots(1,2, figsize=(15, 6))

for dic in vacuum_dicts:
    axs[0].plot(dic["PR interp"], dic["W adia interp"] / dic["W interp"], c="tab:orange")
    #plt.plot(dic["PR interp"], dic["W isot interp"] / dic["W interp"], c="tab:orange")

for dic in compr_dicts:
    axs[0].plot(dic["PR interp"], dic["W adia interp"] / dic["W interp"], c="tab:blue")
    #plt.plot(dic["PR interp"], dic["W isot interp"] / dic["W interp"], c="tab:orange")



axs[0].scatter(MOXIE_pr, MOXIE_eta_i, c="tab:green", label="MOXIE scroll", marker="+", s=40)
axs[0].scatter(liu_2009_pr, liu_2009_eta_i, c="tab:red", label="Liu et al 2009, CO2", marker="+", s=40)
axs[0].scatter(cho_pr, cho_eta_i, c="tab:gray", label="Cho et al 1996, air", marker="+", s=40)
axs[0].scatter(sun_pr, sun_eta_i, c="tab:pink", label="Sun et al, R22", marker="1", s=20)
axs[0].scatter(liu_2010_pr, liu_2010_eta_i, c="tab:brown", label="Liu et al 2010, R22", marker="1", s=20)


handles, labels = axs[0].get_legend_handles_labels()
patch = mpatches.Patch(color='tab:orange', label='AirSquared vacuum pumps')
handles.insert(0, patch)
patch = mpatches.Patch(color='tab:blue', label='AirSquared compressors')
handles.insert(0, patch) 


axs[0].set_xscale("log")
axs[0].set_xlim(1, 500)
axs[0].legend(handles=handles)
axs[0].set_ylim(0, 0.85)
ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*100))
axs[0].yaxis.set_major_formatter(ticks)

axs[0].set_title("Scroll compressor isentropic efficiency")
axs[0].set_xlabel("Pressure ratio")
axs[0].set_ylabel("Isentropic efficiency (%)")


for dic in vacuum_dicts:
    axs[1].plot(dic["PR interp"], dic["eta poly interp"], c="tab:orange")

for dic in compr_dicts:
    axs[1].plot(dic["PR interp"], dic["eta poly interp"], c="tab:blue")

axs[1].scatter(MOXIE_pr, polytropic_efficiency_list(MOXIE_pr, MOXIE_eta_i, 1.31), c="tab:green", label="MOXIE scroll", marker="+", s=40)
axs[1].scatter(liu_2009_pr, polytropic_efficiency_list(liu_2009_pr, liu_2009_eta_i, 1.31), c="tab:red", label="Liu et al 2009, CO2", marker="+", s=40)
axs[1].scatter(cho_pr, polytropic_efficiency_list(cho_pr, cho_eta_i), c="tab:gray", label="Cho et al 1996, air", marker="+", s=40)

handles, labels = axs[1].get_legend_handles_labels()
patch = mpatches.Patch(color='tab:orange', label='AirSquared vacuum pumps')
handles.insert(0, patch)
patch = mpatches.Patch(color='tab:blue', label='AirSquared compressors')
handles.insert(0, patch) 

axs[1].legend(handles=handles)

axs[1].set_xscale("log")
axs[1].set_xlim(1, 500)
axs[1].set_ylim(0, 0.8)

axs[1].set_title("Scroll compressor polytropic efficiency")
axs[1].set_xlabel("Pressure ratio")
axs[1].set_ylabel("Polytropic efficiency (%)")

ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*100))
axs[1].yaxis.set_major_formatter(ticks)





fig, ax1 = plt.subplots(figsize=(4.5, 2.83))
fig.subplots_adjust(left=0.125, right=0.95, bottom=0.175, top=0.925)

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":10})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

for dic in vacuum_dicts:
    ax1.plot(dic["PR interp"], dic["eta poly interp"], c="tab:orange")

for dic in compr_dicts:
    ax1.plot(dic["PR interp"], dic["eta poly interp"], c="tab:blue")

ax1.scatter(MOXIE_pr, polytropic_efficiency_list(MOXIE_pr, MOXIE_eta_i, 1.31), c="tab:green", label="MOXIE scroll", marker="+", s=40)
ax1.scatter(liu_2009_pr, polytropic_efficiency_list(liu_2009_pr, liu_2009_eta_i, 1.31), c="tab:red", label="Liu et al 2009, CO2", marker="+", s=40)
ax1.scatter(cho_pr, polytropic_efficiency_list(cho_pr, cho_eta_i), c="tab:gray", label="Cho et al 1996, air", marker="+", s=40)

handles, labels = ax1.get_legend_handles_labels()
patch = mpatches.Patch(color='tab:orange', label='AirSquared vacuum pumps')
handles.insert(0, patch)
patch = mpatches.Patch(color='tab:blue', label='AirSquared compressors')
handles.insert(0, patch) 

ax1.legend(handles=handles, handlelength=1)

ax1.set_xscale("log")
ax1.set_xlim(1, 500)
ax1.set_ylim(0.1, 0.8)

ax1.set_title("Scroll compressor polytropic efficiency", usetex=True)
ax1.set_xlabel("Pressure ratio", usetex=True)
ax1.set_ylabel("Polytropic efficiency (\%)", usetex=True)

ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*100))
ax1.yaxis.set_major_formatter(ticks)

plt.show()