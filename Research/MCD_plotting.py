import matplotlib.pyplot as plt
import xlwings as xw
from scipy import interpolate
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import os


ws = xw.Book(os.path.join(os.path.dirname(os.path.realpath(__file__)), "Mars Climate Database outputs.xlsx"))


def load_property(sheet_name):
    sheet = ws.sheets[sheet_name]

    start_i = 2
    end_i = 26
    data_out = []
    for i in range(start_i, end_i+1):
        data_out.append([float(x) for x in sheet.range("B{}:Z{}".format(i,i)).value])
    
    return(data_out)

def get_envelope(data):
    min_envelope = []
    max_envelope = []
    mean_envelope = []

    for i in range(len(data[0])):
        data_i = [j[i] for j in data]
        min_envelope.append(min(data_i))
        max_envelope.append(max(data_i))
        mean_envelope.append(np.mean(data_i))
    
    return(min_envelope, mean_envelope, max_envelope)

def plot_envelope(data, ax, title, col, labelx=True):
    if "blue" in col:
        desatcol = "#a8bed6"
    elif "orange" in col:
        desatcol = "#f8bf87"
    elif "green" in col:
        desatcol = "#afd5a9"    
    x = list(range(0, 25))
    data_min, data_avg, data_max = get_envelope(data)
    ax.fill_between(x, data_min, data_max, color=desatcol)
    ax.plot(x, data_avg, color=col)
    ax.set_title(title, usetex=True)
    if labelx:
        ax.set_xlabel("Solar time (hr)", usetex=True, fontsize=12)
    ax.set_xlim(0, 24)
    ax.set_xticks([0,6,12,18,24])
    ax.set_xticklabels([0,6,12,18,24], usetex=True)
    ax.grid(True)

fig, axs = plt.subplots(nrows=2, ncols=3, sharey="row", sharex="col", figsize=(5.76,3.8))

plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular'
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

fig.subplots_adjust(wspace=0.2)
fig.subplots_adjust(hspace=0.225)


plot_envelope(load_property("Erebus T"), axs[0,0], "Erebus temperature", "tab:blue",labelx=False)
plot_envelope(load_property("Erebus P"), axs[1,0], "Erebus pressure", "tab:blue")

plot_envelope(load_property("Arcadia T"), axs[0,1], "Arcadia temperature", "tab:orange",labelx=False)
plot_envelope(load_property("Arcadia P"), axs[1,1], "Arcadia pressure", "tab:orange")

plot_envelope(load_property("Phlegra T"), axs[0,2], "Phlegra temperature", "tab:green",labelx=False)
plot_envelope(load_property("Phlegra P"), axs[1,2], "Phlegra pressure", "tab:green")


axs[0,0].set_ylabel("Temperature (K)", usetex=True, fontsize=12)
axs[0,0].set_yticklabels([160, 180, 200, 220, 240, 260], usetex=True, fontsize=10)
axs[1,0].set_ylabel("Pressure (Pa)", usetex=True, fontsize=12)
axs[1,0].set_yticklabels([500, 600, 700, 800, 900], usetex=True,fontsize=10)

plt.show()