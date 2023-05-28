from design2 import *
from component import *


"""mass_flow = 0.0225
design_T_in = 185.46
design_p_in = 775.79
hours = [True, True, True, False, False, False, False, True]
condense = True
scroll_p_out = 573e3


N_radial = 15
PR_guess = 1.502

ESM_vector = [1, 0.149, 0.121, 100]

all_machines = optimise_machine_design(
                mass_flow, design_T_in, design_p_in, scroll_p_out, PR_guess, N_radial, ESM_vector, condense, verbose=False)

total_ESM = 0
total_weight = 0
total_power = 0
compr_weight = 0

prntcnt = 0
for m in all_machines:
    if type(m) == RadialStage:
        compr_weight += m.weight
    total_ESM += ESM(m.weight, m.power, m.cooling, ESM_vector)


print("#####")
print(total_ESM)
print(compr_weight)
#run_machines_day(all_machines, mass_flow, hours, design_T_in, design_p_in, verbose=False)"""


"""
import xlwings as xw
import matplotlib.pyplot as plt
import os


ws = xw.Book(os.path.join(os.path.dirname(os.path.realpath(__file__)), "Test case comparison.xlsx"))

sheetname = "Combined"

axstates = {}
centstates = {}
hxstates = {}
allstates = {}

for i in ws.sheets[sheetname].range("A3:C6").value:
    if i[2] != None:
        axstates[i[0]] = i[1:]
        allstates[i[0]] = i[1:]
for i in ws.sheets[sheetname].range("A8:C13").value:
    if i[2] != None:
        centstates[i[0]] = i[1:]
        allstates[i[0]] = i[1:]
for i in ws.sheets[sheetname].range("A15:C17").value:
    if i[2] != None:
        hxstates[i[0]] = i[1:]
        allstates[i[0]] = i[1:]


y_pos_ax = np.arange(0, len(axstates.keys()), 1)
y_pos_cent = np.arange(0, len(centstates.keys()), 1) + len(axstates.keys())
y_pos_hx = np.arange(0, len(hxstates.keys()), 1) + len(axstates.keys()) + len(centstates.keys())


fig, ax = plt.subplots(figsize=(112/25.4, 84/25.4))
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "errorbar.capsize":3,
    "lines.linewidth":0.75
})

plt.rc('text.latex')
plt.rc('font', family='sans-serif')

fig.tight_layout()
#fig.subplots_adjust(left=0.1, bottom=0.1)
fig.subplots_adjust(top=0.9, bottom=0.15)
ax.xaxis.grid(True)

ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('top')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.barh(y_pos_ax, [axstates[i][0]*100 for i in axstates.keys()], xerr = [axstates[i][1]*100 for i in axstates.keys()])
ax.barh(y_pos_cent, [centstates[i][0]*100 for i in centstates.keys()], xerr = [centstates[i][1]*100 for i in centstates.keys()])
ax.barh(y_pos_hx, [hxstates[i][0]*100 for i in hxstates.keys()], xerr = [hxstates[i][1]*100 for i in hxstates.keys()])
ax.set_yticks(np.arange(0, 12, 1))
ax.set_yticklabels([i for i in allstates.keys()], usetex=True, fontsize=10)

ax.set_xticks([-60,-40,-20,0,20,40])
ax.set_xticklabels(["-60","-40","-20","0","20","40"], usetex=True, fontsize=10)

ax.set_title(r"Model errors vs experimental results", usetex=True, fontsize=12)
ax.set_xlabel(r"Error \%", usetex=True, fontsize=10)


ax.invert_yaxis()

plt.show()

"""




"""

import xlwings as xw
import matplotlib.pyplot as plt
import os


ws = xw.Book(os.path.join(os.path.dirname(os.path.realpath(__file__)), "Optimiser sensitivity.xlsx"))

sheetname = "Sheet1"

states = {}
for i in ws.sheets[sheetname].range("A2:H24").value:
    if i[2] != None:
        states[i[0]] = i[1:]


boundary_modifiers = []
boundary_modifiers_small = []

model_modifiers = []
model_modifiers_small = []

for key in states.keys():
    if key != "Baseline":
        mod = 100 * (states[key][1] / states["Baseline"][1] - 1)
        states[key].append(mod)
    
        if "Output" in key or "Input" in key:
            boundary_modifiers.append([abs(mod), mod,key])
            if abs(mod) < 2:
                boundary_modifiers_small.append([abs(mod), mod,key])
        else:
            model_modifiers.append([abs(mod), mod,key])
            if abs(mod) < 2:
                model_modifiers_small.append([abs(mod), mod,key])

y_pos_b = np.arange(len(boundary_modifiers))
y_pos_m = np.arange(len(model_modifiers)) + len(boundary_modifiers)

boundary_modifiers.sort()
model_modifiers.sort()



fig, ax = plt.subplots(figsize=(2.96, 2.22))
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":10
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

fig.tight_layout()
fig.subplots_adjust(left=0.1, bottom=0.175)
ax.xaxis.grid(True)

ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('top')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.barh(y_pos_b, [i[1] for i in boundary_modifiers], align="center")
ax.barh(y_pos_m, [i[1] for i in model_modifiers], align="center")

ax.set_yticks(np.concatenate((y_pos_b,y_pos_m)))
ax.set_yticklabels([i[2] for i in boundary_modifiers+model_modifiers], usetex=True, fontsize=8, fontdict={"horizontalalignment":"center"})
ax.set_xticks([-4,-2,0,2,4,6,8])
ax.set_xticklabels(["-4","-2","0","2","4","6","8"], usetex=True, fontsize=10)
ax.set_title(r"ESM (kg-eq) sensitivity analysis", usetex=True, fontsize=12)
ax.set_xlabel(r"ESM change \%", usetex=True, fontsize=10)



fig, ax = plt.subplots(figsize=(2.96, 2.22))
plt.rcParams.update({
    "text.usetex": True,
    'mathtext.default': 'regular',
    "font.size":10
})
plt.rc('text.latex')
plt.rc('font', family='sans-serif')

fig.tight_layout()
fig.subplots_adjust(left=0.1, bottom=0.175)
ax.xaxis.grid(True)

y_pos_b = np.arange(len(boundary_modifiers_small))
y_pos_m = np.arange(len(model_modifiers_small)) + len(boundary_modifiers_small)

boundary_modifiers_small.sort()
model_modifiers_small.sort()

ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('top')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.barh(y_pos_b, [i[1] for i in boundary_modifiers_small], align="center")
ax.barh(y_pos_m, [i[1] for i in model_modifiers_small], align="center")

ax.set_yticks(np.concatenate((y_pos_b,y_pos_m)))
ax.set_yticklabels([i[2] for i in boundary_modifiers_small+model_modifiers_small], usetex=True, fontsize=8, fontdict={"horizontalalignment":"center"})

ax.set_title(r"ESM (kg-eq) sensitivity analysis", usetex=True, fontsize=12)
ax.set_xlabel(r"ESM change \%", usetex=True, fontsize=10)

plt.show()

"""