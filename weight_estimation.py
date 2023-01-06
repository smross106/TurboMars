import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors


def axial_compressor_weight_1(mean_radius, tip_radius, speed_rpm, N_stages):
    # Sagerser scaling
    # https://ntrs.nasa.gov/api/citations/19720005136/downloads/19720005136.pdf

    diameter = mean_radius * 2
    U_tip = speed_rpm * np.pi / 30 * tip_radius
    
    L_over_D = 0.2 + (0.234  - 0.218 * ((2 * mean_radius - tip_radius) / tip_radius)) * N_stages

    weight = (24.2 * np.power(diameter, 2.2) * np.power(N_stages, 1.2) *
                np.power(U_tip / 335, 0.5) * (1 + (L_over_D / (0.2 + 0.081 * N_stages))) )

    # Factor for additional structure
    weight *= 1.18
    
    return(weight)

def axial_compressor_weight_2(mean_radius, tip_radius, solidity, speed_rpm, N_stages, pressure, yield_stress, density):
    
    # From Sagerser scaling
    stage_length = (0.2 + (0.234  - 0.218 * ((2 * mean_radius - tip_radius) / tip_radius)) * N_stages) * mean_radius * 2

    # Casing weights
    hub_radius = 2*mean_radius - tip_radius
    # Assume both inner and outer casings are designed around SF=3 from hoop stress or 3mm, whichever greater
    min_casing_thickness = pressure * tip_radius / (yield_stress)
    casing_thickness = max(min_casing_thickness, 3e-3)
    casing_weight = 2 * np.pi * stage_length * casing_thickness * (tip_radius + hub_radius) * density

    # Blade weights
    # Assume 100% of stage length has a blade in it -> sets chord (balances inter-stage spacings with blade stagger)
    # Blade number per stage set from solidity
    # Assume uniform thickness, equal to 5% chord
    chord = stage_length / (N_stages*2)
    n_blades = 2 * N_stages * 2 * np.pi * mean_radius / chord
    volume_blades = n_blades * (chord * 0.05*chord * (tip_radius - hub_radius))
    blades_weight = volume_blades * density

    print(casing_weight, blades_weight, n_blades)

def radial_compressor_weight(radius_inlet_hub, radius_inlet_tip, radius_outlet_tip, bladeheight_outlet, yield_stress, density, pressure):
    # Best-in-class is NASA WATE code, restricted to US citizens only
    # Instead we shall approximate

    # Assume both casing and volute are quarter-circle revolutions
    revolved_radius_tip = radius_outlet_tip - radius_inlet_tip
    revolved_radius_hub = revolved_radius_tip + bladeheight_outlet
    revolved_tip_centroid = (1 - 4/(3*np.pi)) * revolved_radius_tip + radius_inlet_tip
    #revolved_hub_centroid = (1 - 4/(3*np.pi)) * revolved_radius_hub + radius_inlet_hub

    # Uses Pappus's Centroid Theorem to get surface area  tip
    surface_area_tip = 2 * np.pi * revolved_tip_centroid * (np.pi/2) * revolved_radius_tip

    # Find minimum thickness to support a hoop stress based on pressure and the 
    # arithmetic mean of inlet and outlet tip radii

    thickness_minimum = pressure * (radius_inlet_hub + radius_inlet_tip) / (2*yield_stress)

    # Casing thickness is same but with safety factor 2, 2mm maximum thickness
    casing_thickness = max(2e-3, thickness_minimum*2)

    casing_weight = ((surface_area_tip * casing_thickness)) *  density

    # Also by Pappus centroid theorem
    hub_section_area = (revolved_radius_hub * (revolved_radius_hub+radius_inlet_hub)) - (0.25 * np.pi * revolved_radius_hub**2)

    revolved_hub_centroid = ((revolved_radius_hub+radius_inlet_hub) * revolved_radius_hub * 0.5 * (revolved_radius_hub+radius_inlet_hub)) 

    revolved_hub_centroid -= (radius_inlet_hub + (revolved_radius_hub *(1 - 4/(3*np.pi)))) * (0.25 * np.pi * revolved_radius_hub**2)
    revolved_hub_centroid /= hub_section_area

    hub_volume = (2 * np.pi * revolved_hub_centroid) * hub_section_area
    hub_weight = hub_volume*density

    # Add factor for diffuser plane
    # Assume radius increases by 1.5
    radius_diffuser_out = radius_outlet_tip * 1.5
    diffuser_weight = np.pi * (radius_diffuser_out**2 - radius_outlet_tip**2) * (2*casing_thickness) * density

    # Add fixings approximation
    compressor_weight = (casing_weight + diffuser_weight + hub_weight) * 1.18

    return(compressor_weight)

def cryocooler_weight(cooling_power, optimistic=False):
    if optimistic:
        mass_per_power = 0.1
        figure_of_merit = 0.25
    else:
        mass_per_power = 0.2
        figure_of_merit = 0.125
    
    electrical_power = cooling_power / figure_of_merit
    cryocooler_weight = cooling_power * mass_per_power

    return(cryocooler_weight, electrical_power)

def air_heat_exchanger_weight(cooling_power, air_speed, air_density, REL_values=False):
    # Using simplified methodology from https://intech-gmbh.com/heat_exchangers_calc_and_select/
    # Heat exchanger tubes assumed to be 5mm tubes with 0.5mm wall thickness, made of copper
    # Good design - delta T assumed to be 10 degrees throughout
    # More advanced model may deploy https://arc.aiaa.org/doi/pdf/10.2514/6.2021-3711
    # Assume perfect heat exchange on liquid side, Nusselt number on air side scales as Re^0.8

    if REL_values:
        diameter = 2e-3
        wall_thickness = 4e-5
    else:
        diameter = 5e-3
        wall_thickness = 0.5e-3

    Re = air_speed * air_density * diameter / 14.4e-6
    Nu_air = 0.029 * np.power(Re, 0.8) * 0.9244
    h = Nu_air * 0.011 /  diameter

    delta_T = 10
    area = cooling_power / (h * delta_T)

    hx_weight = area * wall_thickness * 8960

    print("a", area)

    return(hx_weight)

def cryocooler_structure_weight(cooling_power, cycle_length_hours, yield_stress, density):
    # Use Lukas Schrenk scaling law http://lukas-schrenk.com/files/MasterThesis_Schrenk.pdf
    
    # Assume 80% of cooling power goes to CO2 collection, and that CO2 forms in a sphere
    m_co2 = cooling_power * (cycle_length_hours * 3600) * 0.8 / (598000 + 49000)
    R_co2 = np.power((m_co2 / 1500) * 4/(3*np.pi), 1/3)

    # Based on assumed scalings. Temperature set for SF=2 at 20 bar
    R_chamber = 1.5 * R_co2
    H_chamber = 3 * R_co2
    thickness = max(3e-3, 2* 20e5 * R_chamber / (yield_stress))

    chamber_weight = 2*np.pi*thickness * (R_chamber+thickness)**2 + np.pi*H_chamber*(2*np.pi*thickness + thickness**2)
    chamber_weight*=density
    return(chamber_weight)

def motor_weight(mechanical_power, hts=False):
    # Models from SUAVE software used
    # https://www.researchgate.net/publication/314090333_SUAVE_An_Open-Source_Environment_for_Conceptual_Vehicle_Design_and_Optimization
    if hts == False:
        # Non superconducting motor
        # Downstream sources: https://ntrs.nasa.gov/api/citations/20150000747/downloads/20150000747.pdf used
        #                    https://apps.dtic.mil/sti/pdfs/AD1043171.pdf
        motor_weight = 0.888 * np.power(mechanical_power/1000, 0.8997)
    else:
        # Motor with high-temperature superconducting coils
        # Technically outperforms regular motor above 1.9kW but additional mass from radiators may outweigh this saving
        motor_weight = 1.034 * np.power(mechanical_power/1000, 0.6616)
    
    return(motor_weight)

def piston_compressor_weight():
    pass
"""
ax_comp_w = axial_compressor_weight_1(0.1, 0.15, 3000, 8)
print(radial_compressor_weight(0, 0.15, 0.3, 0.1, 300e6, 2700, 1e4))

print(air_heat_exchanger_weight(2500, 50, 0.6))
#print(ax_comp_w)
#print(rad_comp_w)

cc_cooling_power = (598000 + 49000) * 1000 / (8*3600)
cc_w, cc_electrical_power = cryocooler_weight(cc_cooling_power)
cc_s_w = 2*cryocooler_structure_weight(cc_cooling_power, 4, 300e6, 2700)
print(cc_w, cc_s_w)
print(cc_electrical_power)"""

mass_flows = np.linspace(0.005, 0.1, 20)

mass_flow = 0.05

# Solve for compression system
stage_densities_in = [0.0260, 0.0444, 0.0904, 0.211, 0.474, 1.10, 2.46, 5.83]
stage_densities_out = [0.0444, 0.0839, 0.1698, 0.4007, 0.911, 2.081, 4.693, 10.997]

stage_pressures = [900, 1890, 4310, 9825, 22400, 51070, 116500, 265000, 6370000]
stage_efficiencies = [0.85, 0.73, 0.73, 0.73, 0.71, 0.69, 0.67, 0.65]

#0.05kg/s
stage_speeds = [7500, 15000, 15000, 40000, 40000, 40000, 60000, 60000]
stage_radius_out = [0, 250e-3, 175e-3, 60e-3, 60e-3, 60e-3, 40e-3, 37.5e-3]
flow_coeff = [1.0, 0.5, 0.5, 0.5, 0.2, 0.1, 0.07, 0.04]
stage_type = ["A", "C", "C", "C", "C", "C", "C", "C"]

#stage_speeds = [50000, 100000, 100000, 100000, 100000, 100000, 100000, 100000]

# Stage 1 is axial, 8 stages
volume_flow = mass_flow/stage_densities_in[0]
U_tip = np.sqrt(180 * 1.31 * 189) * 0.7
Vx = U_tip * 0.5
area_1 = volume_flow / Vx
# Set blade speed at 7500rpm, tip speed at 0.7 Mach
R_tip = U_tip / (stage_speeds[0] * np.pi/30)
R_mean = 0.5 * (R_tip + np.sqrt(R_tip**2 - area_1/np.pi))
mass_stage_1 = axial_compressor_weight_1(R_mean, R_tip, stage_speeds[0], 8)
work_stage_1 = mass_flow * (stage_pressures[1] - stage_pressures[0]) / stage_densities_in[0]
work_stage_1 /= stage_efficiencies[0]
print("Stage 1 mean radius: {:.1f}mm, bladeheight {:.1f}mm".format(R_mean/1e-3, 2*(R_tip-R_mean)/1e-3))
print("Reynolds number 1st stage: ", int(R_mean * (stage_speeds[0] * np.pi/30) * R_tip**2 / 14.4e-6))

stage_masses = [mass_stage_1]
stage_works = [work_stage_1]
HX_masses = []

for stage in range(1, len(stage_densities_in)):

    if stage_type[stage] == "C":
        # Centrifugal stage
        R_tip_out = stage_radius_out[stage]
        U_tip = R_tip_out * (stage_speeds[stage] * np.pi/30)

        delta_h = U_tip**2 * 0.85

        Vx = U_tip * flow_coeff[stage]

        area_in = mass_flow / (stage_densities_in[stage] * Vx)
        area_out = mass_flow / (stage_densities_out[stage] * Vx)
        # Use specific speed to set R_mean
        #delta_h = U_tip**2 * 0.85
        #R_tip_out = U_tip / (stage_speeds[stage] * np.pi/30)
        

        R_mean_in = 0.5 * 4.5 * mass_flow / (np.power(delta_h, 0.25) * stage_densities_in[stage])
        if area_in / (4 * np.pi * R_mean_in) > R_mean_in:
            R_tip_in = np.sqrt(area_in/np.pi)
            R_mean_in = R_tip_in/2
            R_hub_in = 0
            bladeheight_in = R_tip_in
        else:
            y = area_in / (4 * np.pi * R_mean_in)
            bladeheight_in = y*2
            R_tip_in = R_mean_in + y
            R_hub_in = R_mean_in - y

        
        bladeheight_out = area_out / (2 * np.pi * R_tip_out)

        stage_mass = radial_compressor_weight(R_hub_in, R_tip_in, R_tip_out, bladeheight_out, 500e6, 2700, stage_pressures[stage])
        stage_power = delta_h * mass_flow / stage_efficiencies[stage]

        stage_masses.append(stage_mass)
        stage_works.append(stage_power/0.8)

        print("Stage {} inlet radius: {:.1f}mm, outlet radius {:.1f}mm, outlet bladeheight {:.2f}mm, delta h {:.1f}kJ/kg".format(
            stage+1, R_tip_in/1e-3, R_tip_out/1e-3, bladeheight_out/1e-3, delta_h/1000))
        #print("\t delta h {:.1f}kJ/kg".format(delta_h/1000))

    hx_mass = air_heat_exchanger_weight(delta_h * mass_flow, Vx, stage_densities_out[stage])
    print(delta_h*mass_flow)
    HX_masses.append(hx_mass)

motor_weights = [motor_weight(stage_works[0]),
    motor_weight(stage_works[1]),
    motor_weight(stage_works[2]),
    motor_weight(stage_works[3]+stage_works[4]+stage_works[5]),
    motor_weight(stage_works[6]+stage_works[7]),]

print(stage_masses)
print("Turbo mass: ", int(sum(stage_masses) + sum(HX_masses) + sum(motor_weights)))
#print(stage_works)
print("Turbo power: ", int(sum(stage_works)))

cryocooler_cooling = (598000 + 49000) * mass_flow

cc_w, cc_electrical_power = cryocooler_weight(cryocooler_cooling)
cc_s_w = 2*cryocooler_structure_weight(cryocooler_cooling, 4, 300e6, 2700)
print("Cryo weight: ", int(cc_w + cc_s_w))
print("Cryo power: ", int(cc_electrical_power))



cmap1 = plt.cm.cool # Axial stages
cmap2 = plt.cm.Purples # Radial stages
cmap3 = plt.cm.Greens # Heat exchangers
cmap4 = plt.cm.Greys # Low speed motors
cmap5 = plt.cm.Reds # High speed motors


inner_colours = [cmap1(.8), cmap2(.8), cmap3(.8), cmap4(.5), cmap5(.8)]

outer_colours = [
    cmap1(.8),
    *cmap2(np.linspace(.8, .2, 7)),
    *cmap3(np.linspace(.8, .2, 7)),
    *cmap4(np.linspace(.5, .2, 3)),
    *cmap5(np.linspace(.8, .6, 2))
]
#outer_colors = cmap(np.arange(7)*6)
#inner_colors = [cmap(np.linspace(.6, .1, 9)),
#                cmap(np.linspace(.6, .2, 4)),
#                cmap(np.linspace(.6, .2, 2))]

crude_labels = ["Axial stages", "Radial stages", "Heat exchangers", "Low-speed motors", "High-speed motors"]
detailed_labels = ["Axial stages", 
    "Radial stage 1", "Radial stage 2", "Radial stage 3", "Radial stage 4", "Radial stage 5", "Radial stage 6", "Radial stage 7",
    "HX 1", "HX 2", "HX 3", "HX 4", "HX 5", "HX 6", "HX 7",
    "Axial 7500RPM motor", "Radial 7500RPM motor", "Radial 15kRPM motor", "Radial 40kRPM motor", "Radial 50kRPM motor"]
plt.pie(stage_masses+HX_masses+motor_weights, radius=1, colors=outer_colours, labels=detailed_labels, rotatelabels=True)
plt.pie([stage_masses[0], sum(stage_masses[1:]), sum(HX_masses), sum(motor_weights[0:3]), sum(motor_weights[3:])], radius=0.5, autopct='%1.1f%%',colors=inner_colours)

plt.show()

