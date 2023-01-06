import component

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import time
import design

flow_cold = component.GasFlow(0.05, 250, 274787)
flow_hot = component.GasFlow(0.05, 300.5, 549574)


# ax1 = design.optimise_axial(flow, flow.delta_h_PR(1.1), 0.85)

rad1 = design.optimise_radial(flow_cold, flow_cold.delta_h_PR(2), 0.85)
print(rad1)
print(rad1.A_stator_exit)


#hx1 = design.optimise_hx(flow_hot, flow_hot.delta_h_delta_T(abs(250-flow_hot.temperature)), 3.448283283164784e-05, flow_hot.pressure*0.03)
hx_cond = design.optimise_hx(flow_hot, 72.4e3, 0.00016509001192732648, flow_hot.pressure*0.03)
print(hx_cond)
"""print(hx1.duct_length, hx1.diffuser_length, hx1.nozzle_length)
print(hx1.pitch_diameter_ratio*hx1.tube_diameter*1000)
print(hx1.bulk_area, hx1.coolant_velocity, hx1.pitch_diameter_ratio)
print(hx1.required_area)
print(hx1.n_tubes_wide)"""

