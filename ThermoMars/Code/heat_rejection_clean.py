import numpy as np
import matplotlib.pyplot as plt

"""
References:
 - [1] - A solar azimuth formula that renders circumstantial treatment unnecessary without compromising mathematical rigor: 
        Mathematical setup, application and extension of a formula based on the subsolar point and atan2 function, Zhang et al 2021
 - [2] - Über die Extinktion des Lichtes in der Erdatmosphäre, Schoenberg 1929
 - [3] - Mars Fact Sheet, https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
 - [4] - Thermal control of MSL Rover "Curiosity" using an Active Fluid Loop, Birur 2013 
 - [5] - Advanced Active Thermal Control Systems Architecture Study, Hanford and Ewert
 - [6] - Electric feed systems for liquid propellant rockets, Rachov et al
 - [7] - Mass Optimization and Experimental Validation of a Forced-Convection Heat Exchanger for Mars Surface Waste 
        Heat Rejection, Colgan 2023

"""

SOL_HRS = 24.65
MARS_RADIUS = 3396.2e3
MARS_ATMOSPHERE_OPTICAL_HEIGHT = 8624
MARS_ECC = 0.093377

pump_efficiency = 0.68
pump_motor_power_density = 1/3240 # kg/W
# Values from rocket turbopumps, [6]

heat_exchanger_power_density = 2.72/1000 # kg/W
heat_pump_efficiency = 0.50
# From [5]

REL_heat_exchanger_power_density = heat_exchanger_power_density/100
# From some friends at Reaction Engines who "can beat your numbers by a factor of 100"

def generate_solar_data(num_points, times, latitude, areocentric_longitude, atmosphere_tau, baseline_solar_heat=590, hours_per_sol=24.936):
        """
        Generate lists of solar intensities, altitudes and azimuths for a single day.

        Solar position implements method from Reference [1] with necessary modifications for Martian orbital parameters.

        Martian atmospheric optimal height, for use in calculation from Source [2]. Calculation is based on constant-density
        hydrostatic argument, with data taken from [3]
        """
        
        axis_obliquity = np.deg2rad(hours_per_sol)
        subsolar_latitude = np.arcsin(np.sin(axis_obliquity) * np.sin(np.deg2rad(areocentric_longitude)))

        solar_intensity_space = baseline_solar_heat * np.power(
        (1 + MARS_ECC * np.cos(np.deg2rad(areocentric_longitude - 248))) / (1 - MARS_ECC), 2)

        solar_intensities = np.zeros((num_points))
        solar_altitudes = np.zeros((num_points))
        solar_azimuths = np.zeros((num_points))
        solar_vectors = np.empty((num_points, 3))


        for point in range(num_points):
            subsolar_longitude = np.deg2rad(-15 * (times[point] - (SOL_HRS/2)))

            Sx = np.cos(subsolar_latitude) * np.sin(subsolar_longitude)
            Sy = (np.cos(latitude) * np.sin(subsolar_latitude)) - (
                np.sin(latitude) * np.cos(subsolar_latitude) * np.cos(subsolar_longitude))
            Sz = (np.sin(np.deg2rad(latitude)) * np.sin(subsolar_latitude)) + (
                np.cos(np.deg2rad(latitude)) * np.cos(subsolar_latitude) * np.cos(subsolar_longitude))

            solar_altitude = np.rad2deg(np.arcsin(Sz))
            solar_azimuth_south_clockwise = np.rad2deg(np.arctan2(-Sx, Sy))
            solar_azimuth = solar_azimuth_south_clockwise


            # Calculate relative air mass for scaling the optical depth of the atmosphere
            relative_air_mass = np.sqrt(np.power(MARS_RADIUS + MARS_ATMOSPHERE_OPTICAL_HEIGHT, 2) - 
                np.power(MARS_RADIUS * np.cos(np.deg2rad(solar_altitude)), 2)) - MARS_RADIUS * np.sin(np.deg2rad(solar_altitude))
            relative_air_mass /= MARS_ATMOSPHERE_OPTICAL_HEIGHT

            optical_depth = np.exp(- atmosphere_tau * relative_air_mass)

            solar_intensity = solar_intensity_space * optical_depth

            if Sz <= 0:
                solar_intensity = 0
                solar_altitude = 0

            solar_intensities[point] = solar_intensity
            solar_altitudes[point] = solar_altitude
            solar_azimuths[point] = solar_azimuth

            # x is east, y north, z up
            solar_vectors[point] = [Sx, Sy, Sz]
            solar_vectors[point] = solar_vectors[point]/np.linalg.norm(solar_vectors[point])

        return(solar_intensities, solar_vectors)

def get_insolation(panel_vector, time, latitude, areocentric_longitude=270, atmosphere_tau=0.1):
    solar_intensities, solar_vectors = generate_solar_data(1, [time], latitude, areocentric_longitude, atmosphere_tau)

    solar_intensity = solar_intensities[0]
    solar_vector = solar_vectors[0]
    
    panel_vector_normal = panel_vector/np.linalg.norm(panel_vector)
    insolation_factor = np.dot(panel_vector_normal, solar_vector)
    if insolation_factor<0:
        insolation_factor=0
    return(insolation_factor*solar_intensity)

def radiator_loss_evaluation(heat_load_w, ground_temp, air_temp, radiator_temp, rad_orientation, latitude, time,
                         rad_infrared_emit=0.85, rad_solar_absorb=0.15, 
                         ground_reflect=0.3, ground_emit=0.94, air_emit=0.225, radiator_area_mass=10):
    
    SB_const = 5.67E-8
    
    infrared_out = (radiator_temp ** 4) * SB_const * rad_infrared_emit
    if rad_orientation == "h":
        infrared_in = SB_const * rad_infrared_emit *((air_temp**4) * air_emit)
        indirect_in = 0
    else:
        infrared_in = SB_const * rad_infrared_emit * (
            ((air_temp**4) * air_emit * 0.5) + 
            ((ground_temp**4) * ground_emit * 0.5))
        indirect_in = get_insolation([0,0,1], time, latitude) * ground_reflect * rad_solar_absorb * 0.5
    
    if rad_orientation == "h":
        solar_in = get_insolation([0,0,1], time, latitude) * rad_solar_absorb # horizontal
    
    elif rad_orientation == "ns":
        solar_in = (get_insolation([0,1,0], time, latitude) + get_insolation([0,-1,0], time, latitude)) * rad_solar_absorb # N-S
    elif rad_orientation == "ew":
        solar_in = (get_insolation([1,0,0], time, latitude) + get_insolation([-1,0,0], time, latitude)) * rad_solar_absorb # E-W
    

    radiator_loss_per_m2 = infrared_out - infrared_in - solar_in - indirect_in

    radiator_mass = radiator_area_mass * heat_load_w / radiator_loss_per_m2

    total_radiator_pressure_drop_Pa = 371.8 * heat_load_w / radiator_loss_per_m2
    # Hanford and Ewert [5] quote this number as the ISS requirement, Pa/m2 of radiator area
    total_mass_flow = heat_load_w / (10 * 4750)
    # Approximate properties for ammonia at 70 bar, same as ISS radiators
    # Assume delta T in radiators is 10 degrees K
    pumping_power = total_radiator_pressure_drop_Pa * (total_mass_flow/630)
    pumping_mass = pumping_power * pump_motor_power_density

    return({
        "Total mass": radiator_mass+pumping_mass, 
        "ESM":(radiator_mass+pumping_mass) + 0.1*pumping_power,
        "Loss per m2": radiator_loss_per_m2,
        "Radiator mass": radiator_mass,
        "Power": pumping_power,
        "Aux mass": pumping_mass})
    

def forced_convection_hx(average_coolant_temp, heat_load_w, air_temp, air_pressure):
    # All data digitised and curvefit from [7]
    # Baseline values
    # Steel HX, ambient temperature 260K, 600Pa, average coolant temperature 412.5K
    #if heat_load_w < 1000:
    #    heat_load_w = 1000
        # Model does not cover low power systems, assume 1kWth is the minimum size
    steel_hx_mass = 4.1715E-4 * np.power(heat_load_w, 1.0204)
    hx_fan_power_w = heat_load_w / 24
    hx_pumping_power_w = 1.4348E-9 * np.power(heat_load_w, 2.0914)

    delta_T = average_coolant_temp - air_temp
    air_density_scaled = air_pressure / air_temp
    #ambient_scale_factor = np.exp(
    #    0.0429 - 0.516*np.log(delta_T/152.5) - 0.66*np.log(air_density_scaled / (600/260))
    #)
    ambient_scale_factor = 1.636 - 0.562*(air_density_scaled / (600/260))
    ambient_scale_factor /= (delta_T/152.5)

    
    if average_coolant_temp > 600:
        material_scale_factor = 1.0
    else:
        material_scale_factor = 0.8
    
    hx_mass = steel_hx_mass * ambient_scale_factor * material_scale_factor
    power = hx_fan_power_w+hx_pumping_power_w
    pump_mass = hx_pumping_power_w * pump_motor_power_density

    #return(hx_mass + 0.15*(hx_fan_power_w+hx_pumping_power_w), hx_fan_power_w+hx_pumping_power_w)
    return({
        "Total mass": hx_mass+pump_mass,
        "ESM": (hx_mass+pump_mass) + 0.1*power,
        "Power": power,
        "Aux mass": pump_mass
    })


def compressor_mass(power_kw):
    mass = 31.83 * np.power(power_kw, 0.476)

    return(mass)


def ESM_radiator_curvefit(T_reject):
    # kg/W
    delta_T = T_reject - 235
    if delta_T < 0:
        return(np.nan)
    one_over_ESM = (5.858E-6 * np.power(delta_T, 2)) - (1.004E-4 * delta_T) + (1.498E-2)
    ESM = 1/one_over_ESM

    return(ESM/1000)

def ESM_FCHX_curvefit(T_reject):
    # kg/W
    delta_T = T_reject - 235
    if delta_T < 0:
        return(np.nan)
    one_over_ESM = 2.853E-2 * np.log(delta_T) + 7.468E-2
    ESM = 1/one_over_ESM

    return(ESM/1000)

def increase_T_heatpump(T_input, T_output, heat_in, power_ESM = 0.1, use_REL=False, eta=0.5):
    carnot = T_output / (T_output-T_input)
    
    CoP = carnot * eta

    work_in = heat_in / CoP
    heat_out = heat_in + work_in

    compressor_weight = compressor_mass(work_in/1000)
    if not use_REL:
        hx_in_weight = heat_in * heat_exchanger_power_density
        hx_out_weight = heat_out * heat_exchanger_power_density
    else:
        hx_in_weight = heat_in * REL_heat_exchanger_power_density
        hx_out_weight = heat_out * REL_heat_exchanger_power_density
    
    return(compressor_weight+hx_in_weight+hx_out_weight+work_in*power_ESM, heat_out, work_in, CoP)