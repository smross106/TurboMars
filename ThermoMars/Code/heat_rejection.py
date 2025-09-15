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

def generate_solar_data(num_points, times, latitude, areocentric_longitude, atmosphere_tau):
        """
        Generate lists of solar intensities, altitudes and azimuths for a single day.

        Solar position implements method from Reference [1] with necessary modifications for Martian orbital parameters.

        Martian atmospheric optimal height, for use in calculation from Source [2]. Calculation is based on constant-density
        hydrostatic argument, with data taken from [3]
        """
        
        axis_obliquity = np.deg2rad(24.936)
        subsolar_latitude = np.arcsin(np.sin(axis_obliquity) * np.sin(np.deg2rad(areocentric_longitude)))

        solar_intensity_space = 590 * np.power(
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
    

#times = [6.16, 7.39, 8.63, 9.86, 11.09, 12.32, 13.56, 14.79, 16.02, 17.26, 18.49]

def zigzag_temp_function(time):
    avg_ground_temp = 225
    ground_temp_swing = 70
    avg_air_temp = 205
    air_temp_swing = 6.05
    max_temp_time = 13.35
    min_temp_time = 6.16
    day_length = 24.66

    if time < min_temp_time:
        time_since_max = (day_length - max_temp_time) + time
        ratio = 1 - (time_since_max / (day_length - max_temp_time + min_temp_time))
    elif time >= min_temp_time and time < max_temp_time:
        ratio = (time - min_temp_time) / (max_temp_time - min_temp_time)
    else:
        ratio = 1 - ((time - max_temp_time) / (day_length - max_temp_time + min_temp_time))
    
    ground_temp = avg_ground_temp - (ratio - 0.5) * ground_temp_swing
    air_temp = avg_air_temp - (ratio - 0.5) * air_temp_swing

    return(ground_temp, air_temp)


def sin_linear_temp_function(time, T_peak, T_min, time_peak, latitude, areocentric_longitude=270):
        """
        Generate the temperatures at each point in the day - can be called with either ground or ambient

        Model used is a sinusoidal fit during the day, rising from a minimum at sunset to a maximum
            either at local noon, or time_peak if set
        Between sunset and sunrise, a linear fit is used to bridge the temperatures.
        This approximately matches the data from [4]
        """
        
        axis_obliquity = np.deg2rad(24.936)
        solar_declination = np.arcsin(np.sin(axis_obliquity) * np.sin(np.deg2rad(areocentric_longitude)))

        time_sunrise = (SOL_HRS / (2 * np.pi)) * np.arccos(-np.tan(np.deg2rad(latitude)) * np.tan(solar_declination))
        time_sunset = SOL_HRS - time_sunrise

        # Setting up amplitude for the cos wave
        T_sine_peak = T_peak 
        T_sine_trough = (T_min - T_sine_peak * np.cos((time_sunrise - time_peak) * 2 * np.pi / SOL_HRS)) / (
            1 - np.cos((time_sunrise - time_peak) * 2 * np.pi / SOL_HRS))
        
        def T_sinusoid(time):
            T = (T_sine_peak - T_sine_trough) * (np.cos((time - time_peak) * 2 * np.pi / SOL_HRS)) + T_sine_trough
            return(T)

        # Set up parameters for the linear fit
        T_sunset = T_sinusoid(time_sunset)
        T_sunrise = T_min
        T_midnight = (T_sunset + T_sunrise) / 2

        def T_linear(time):
            if time > time_sunset:
                T = (T_sunset - T_midnight) * ((SOL_HRS - time) / (SOL_HRS - time_sunset)) + T_midnight
            else:
                T = (T_sunrise - T_midnight) * (time / time_sunrise) + T_midnight
            return(T)
        


        if time > time_sunrise and time < time_sunset:
            # Time is during the day
            T = T_sinusoid(time)
        else:
            T = T_linear(time)
            
        
        return(T)


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


def day_evaluate(evaluates, latitude, 
                 max_ground_temp, min_ground_temp, max_air_temp, min_air_temp, max_temp_time, system_temps, 
                 target_heat=1e3, radiator_area_mass=6.93):
    evaluated = []
    valid_evaluates = ["hx", "rh", "rew", "rns"]
    for eval in evaluates:
        if eval in valid_evaluates:
            evaluated.append([])
        else:
            raise ValueError
    
    times = np.linspace(0, 24.6, 50)
    for t in times:
        ground_temp = sin_linear_temp_function(t, max_ground_temp, min_ground_temp, max_temp_time, latitude)
        air_temp = sin_linear_temp_function(t, max_air_temp, min_air_temp, max_temp_time, latitude)
        for eval_i, eval in enumerate(evaluates):
            if eval in valid_evaluates:
                if eval == "hx":
                    hx_data = forced_convection_hx(system_temps[eval_i], target_heat, air_temp, 800)
                    evaluated[eval_i].append(hx_data)
                elif eval[0] == "r":
                    radiator_orientation = eval[1:]
                    radiator_data = radiator_loss_evaluation(
                        target_heat, ground_temp, air_temp, system_temps[eval_i], 
                        radiator_orientation, latitude, t, radiator_area_mass=radiator_area_mass)
                    
                    

                    evaluated[eval_i].append(radiator_data)

    return(times, evaluated)


#print(forced_convection_hx((402+334)/2, 10E3, 260, 600))

"""times, data = day_evaluate(["hx", "rew"], 39, 260, 260-70, 260, 260-40, 13.35, 
                                     [300, 300], target_heat=100E3)
print(data[0][10])
print(data[1][10])"""


peak_temps = np.linspace(200, 300, 50)
radiator_temps = np.linspace(250, 500, 70)

radiator_data_max = np.zeros((peak_temps.size, radiator_temps.size))
radiator_data_avg = np.zeros((peak_temps.size, radiator_temps.size), dtype=object)
radiator_ESM_avg = np.zeros((peak_temps.size, radiator_temps.size))
hx_data_max = np.zeros((peak_temps.size, radiator_temps.size))
hx_data_avg = np.zeros((peak_temps.size, radiator_temps.size), dtype=object)
hx_ESM_avg = np.zeros((peak_temps.size, radiator_temps.size))

for pt_i, peak_temp in enumerate(peak_temps):
    for rt_i, radiator_temp in enumerate(radiator_temps):
        times, output = day_evaluate(["hx", "rew"], 39, peak_temp, peak_temp-70, peak_temp, peak_temp-40, 13.35, 
                                     [radiator_temp, radiator_temp], target_heat=100E3)
        all_hx_ESM = [i["ESM"] for i in output[0] if i["ESM"] > 0]
        all_radiators_ESM = [i["ESM"] for i in output[1] if i["ESM"] > 0]

        if peak_temp < radiator_temp:
            hx_data_max[pt_i, rt_i] = max(all_hx_ESM)
            radiator_data_max[pt_i, rt_i] = max(all_radiators_ESM)
        else:
            hx_data_max[pt_i, rt_i] = np.nan
            radiator_data_max[pt_i, rt_i] = np.nan

        def closest_ESM(data, ESM_target):
            closest = data[0]
            for i in data:
                if abs(i["ESM"] - ESM_target) < closest["ESM"]:
                    closest = i
            return(closest)

        average_hx_ESM = sum(all_hx_ESM)/len(all_hx_ESM)
        hx_closest = closest_ESM(output[0], average_hx_ESM)
        hx_data_avg[pt_i, rt_i] = hx_closest
        hx_ESM_avg[pt_i, rt_i] = hx_closest["ESM"]

        average_radiator_ESM = sum(all_radiators_ESM)/len(all_radiators_ESM)
        radiator_closest = closest_ESM(output[1], average_radiator_ESM)
        radiator_data_avg[pt_i, rt_i] = radiator_closest
        radiator_ESM_avg[pt_i, rt_i] = radiator_closest["ESM"]
    
    print(int(pt_i/peak_temps.size * 100), "%")

print(np.nanmin(hx_data_max), np.nanmax(hx_data_max))
print(np.nanmin(radiator_data_max), np.nanmax(radiator_data_max))
print(np.nanmin(hx_ESM_avg), np.nanmax(hx_ESM_avg))
print(np.nanmin(radiator_ESM_avg), np.nanmax(radiator_ESM_avg))


plt.plot(radiator_temps, [i/100 for i in hx_data_max[17]], label="FCHX")
plt.plot(radiator_temps, [i/100 for i in radiator_data_max[17]], label="Radiator")
plt.ylabel("ESM of heat rejection, kg/kWth")
plt.xlabel("Temperature of heat rejection system")
plt.title("Equivalent System Mass of different heat rejection systems \n 100kWth systems, maximum environment temperature 236K")
plt.legend()
plt.show()



lvls_rad = [200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000]
lvls_hx = [200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000]
lvls_delta = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

fig, ax = plt.subplots()

plt.contourf(peak_temps, radiator_temps, np.transpose(hx_data_max), levels=lvls_hx)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, np.transpose(hx_data_max), levels=lvls_hx)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Forced-convection heat exchanger \n equivalent system mass, 100kWth capacity")

fig, ax = plt.subplots()

plt.contourf(peak_temps, radiator_temps, np.transpose(radiator_data_max), levels=lvls_rad)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, np.transpose(radiator_data_max), levels=lvls_rad)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Lightweight east-west radiator \n equivalent system mass, 100kWth capacity")

fig, ax = plt.subplots()

plt.contourf(peak_temps, radiator_temps, np.transpose(hx_ESM_avg), levels=lvls_hx)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, np.transpose(hx_ESM_avg), levels=lvls_hx)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Forced-convection heat exchanger w/ thermal storage \n equivalent system mass, 100kWth capacity")

fig, ax = plt.subplots()

plt.contourf(peak_temps, radiator_temps, np.transpose(radiator_ESM_avg), levels=lvls_rad)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, np.transpose(radiator_ESM_avg), levels=lvls_rad)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Lightweight east-west radiator w/ thermal storage \n equivalent system mass, 100kWth capacity")

fig, ax = plt.subplots()

radiator_improvement_percent = 100 * (np.transpose(radiator_data_max) - np.transpose(radiator_ESM_avg)) / np.transpose(radiator_data_max)

plt.contourf(peak_temps, radiator_temps, radiator_improvement_percent, levels=lvls_delta)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, radiator_improvement_percent, levels=lvls_delta)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Improvement in radiator ESM (% reduction) \n w/ thermal storage")

fig, ax = plt.subplots()

hx_improvement_percent = 100 * (np.transpose(hx_data_max) - np.transpose(hx_ESM_avg)) / np.transpose(hx_data_max)


plt.contourf(peak_temps, radiator_temps, hx_improvement_percent, levels=lvls_delta)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, hx_improvement_percent, levels=lvls_delta)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Improvement in heat exchanger ESM (% reduction) \n w/ thermal storage")



fig, ax = plt.subplots()

improvement_percent = 100 * (np.transpose(radiator_ESM_avg) - np.transpose(hx_ESM_avg)) / np.transpose(radiator_ESM_avg)
lvls_delta2 = [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

plt.contourf(peak_temps, radiator_temps, improvement_percent, levels=lvls_delta2)
plt.plot([250, 300], [250, 300], c="k")
plt.colorbar()

CS = ax.contour(peak_temps, radiator_temps, improvement_percent, levels=lvls_delta2)
ax.clabel(CS, CS.levels, inline=True, fontsize=10, colors="k")

plt.xlabel("Environment maximum temperature (K)")
plt.ylabel("Heat rejection system temperature (K)")
plt.title("Improvement in heat rejection ESM (% reduction) \n w/ thermal storage, vs radiator")


plt.show()

"""times, data = day_evaluate(["hx", "rew", "rns", "rh"], 25, 260, 190, 222, 180, 13.35, np.linspace(298, 298, 4), 10E3)

plt.plot(times, [i["ESM"] for i in data[0]], label="Heat exchanger")
plt.plot(times, [i["ESM"] for i in data[1]], label="Radiator east-west")
plt.plot(times, [i["ESM"] for i in data[2]], label="Radiator north-south")
plt.plot(times, [i["ESM"] for i in data[3]], label="Radiator horizontal")

print("Average HX mass", sum([i["ESM"] for i in data[0]])/len(times))
print("Average E-W radiator mass", sum([i["ESM"] for i in data[1]])/len(times))
print("Average N-S radiator mass", sum([i["ESM"] for i in data[2]])/len(times))
print("Ratio", min(sum([i["ESM"] for i in data[1]]),sum([i["ESM"] for i in data[2]]))/sum([i["ESM"] for i in data[0]]))

print(data[1][0])
print(data[0][0])

plt.legend()
plt.show()"""