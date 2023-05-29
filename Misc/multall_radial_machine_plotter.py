# Plots machines based on circular arcs

import numpy as np
import matplotlib.pyplot as plt

n_points_machine = 7
n_points_entry = 2
n_points_exit = 2

inlet_radius = 0.15
curvature_radius = 0.35
angle_turned = 90

velocity_scaling = "linear"
velocity_out = 0.6

axial = []
radial = []
velocity = []

theta = 0

for i in range(n_points_machine):
    theta = np.deg2rad((angle_turned / (n_points_machine-1))*i)
    x = np.sin(theta) * curvature_radius
    r = inlet_radius + curvature_radius * (1 - np.cos(theta))

    if velocity_scaling == "linear":
        vel = 1 - (1 - velocity_out) * (i / (n_points_machine + n_points_exit))
    
    axial.append(x)
    radial.append(r)
    velocity.append(vel)


entry_dx = axial[1] - axial[0]
exit_x = axial[-1]
exit_r = radial[-1]
exit_dx = np.cos(theta) * curvature_radius
exit_dr = radial[-1] - radial[-2]

for i in range(n_points_exit):
    x = exit_x + (i+1) * exit_dx
    r = exit_r + (i+1) * exit_dr

    if velocity_scaling == "linear":
        vel = 1 - (1 - velocity_out) * ((i + n_points_machine+1) / (n_points_machine + n_points_exit))
    
    axial.append(x)
    radial.append(r)
    velocity.append(vel)

for i in range(n_points_entry):
    axial.insert(0, axial[0] - entry_dx)
    radial.insert(0, inlet_radius)
    velocity.insert(0, 1.0)

print(axial)
print(radial)

axial_str = ""
radial_str = ""
vel_str = ""

for i in axial:
    axial_str += "{:.5f}".format(i)
    axial_str += "  "

for i in radial:
    radial_str += "{:.5f}".format(i)
    radial_str += "  "

for i in velocity:
    vel_str += "{:.5f}".format(i)
    vel_str += "  "

print("Axial coords")
print(axial_str)
print("Radial coords")
print(radial_str)
print("Velocities")
print(vel_str)