import numpy as np
from scipy import optimize


class GasFlow(object):
    def __init__(self, mass_flow, temperature, pressure, carbon_dioxide=True):
        self.mass_flow = mass_flow
        self.temperature = temperature
        self.pressure = pressure

        self.carbon_dioxide = carbon_dioxide
        # If carbon_dioxide=False, gas is air
        if self.carbon_dioxide:
            self.Rgas = 189
            self.cp = 830
            self.gamma = 1.31
            self.viscosity = 14e-6
            self.Pr = 0.79
            self.k = 0.011

        else:
            self.Rgas = 287
            self.cp = 1005
            self.gamma = 1.40
            self.viscosity = 17e-6
            self.Pr = 0.72
            self.k = 0.016

        self.gm1 = self.gamma - 1

    def density(self):
        return(self.pressure / (self.Rgas * self.temperature))
    
    def density_deltah(self, deltah):
        temperature_new = self.temperature + (deltah / self.cp)
        density_new = np.power(temperature_new/self.temperature, 1/self.gm1) * self.density()

        return(density_new)
    
    def pressure_ratio(self, delta_T, efficiency):
        isentropic_T_ratio = 1 + delta_T / (self.temperature * efficiency)
        pressure_ratio = np.power(isentropic_T_ratio, self.gamma / self.gm1)
        return(pressure_ratio)
    
    def speed_sound(self):
        return(np.sqrt(self.gamma * self.Rgas * self.temperature))
    
    def density_static(self, velocity):
        mach = velocity / self.speed_sound()
        rho_stag = self.density()
        rho_static = rho_stag * np.power(1 + self.gm1/2 * mach**2, 1/self.gm1)

        return(rho_static)
    
    def delta_h_PR(self, PR):
        TR = np.power(PR, self.gm1 / self.gamma)
        delta_T = (TR - 1) * self.temperature
        delta_h = self.cp * delta_T
        return(delta_h)
    
    def delta_h_delta_T(self, delta_T):
        delta_h = delta_T * self.cp
        return(delta_h)

class Component(object):
    def __init__(self):
        pass

class AxialStage(Component):
    def __init__(self, gasflow_in, pressure_ratio_stage, speed, work_coeff, flow_coeff, efficiency_guess):
        self.gasflow_in = gasflow_in

        self.pressure_ratio_stage = pressure_ratio_stage

        self.work_stage = self.gasflow_in.delta_h_PR(self.pressure_ratio_stage) / efficiency_guess

        self.speed = speed
        self.speed_rad = speed * np.pi / 30
        self.work_coeff = work_coeff
        self.flow_coeff = flow_coeff

        self.efficiency = efficiency_guess

        self.R_mean_in = 0
        self.Vx_in = 0
        self.A_in = 0
        self.bladeheight_in = 0
        self.R_tip_in = 0
        self.weight = 0

        self.calculate_geometry_in()

    def calculate_geometry_in(self):

        self.R_mean_in = np.sqrt(self.work_stage / (self.work_coeff * self.speed_rad**2))

        self.Vx_in = self.speed_rad * self.R_mean_in * self.flow_coeff
        self.A_in = self.gasflow_in.mass_flow / (self.Vx_in * self.gasflow_in.density())
        self.bladeheight_in = self.A_in / (2 * np.pi * self.R_mean_in)
        self.R_tip_in = self.R_mean_in + (self.bladeheight_in/2)
        self.R_hub_in = self.R_mean_in - (self.bladeheight_in/2)

        self.A_stator_exit = self.A_in

    def weight_estimate(self):
        diameter = self.R_mean_in * 2
        U_tip = self.R_tip_in * self.speed_rad

        N_stages = 1

        L_over_D = 0.2 + (0.234  - 0.218 * ((2 * self.R_mean_in - self.R_tip_in) / self.R_tip_in)) * N_stages

        weight = (24.2 * np.power(diameter, 2.2) * np.power(N_stages, 1.2) *
                np.power(U_tip / 335, 0.5) * (1 + (L_over_D / (0.2 + 0.081 * N_stages))) )
        
        self.weight = weight

        return(weight)

    def ESM_estimate(self):
        self.power = self.work_stage * self.gasflow_in.mass_flow / self.efficiency
        power_weight = self.power * .121

        self.ESM = self.weight + power_weight

        return(self.ESM)

    def efficiency_cordier(self, verbose=False):
        # Estimate efficiency based on the Cordier plot
        # Assume peak efficiency of 85% at specific speed=2.0
        # Drop-off of 10% by a change of 1.5 in either direction
        specific_speed = (np.sqrt(self.gasflow_in.mass_flow / self.gasflow_in.density_deltah(self.work_stage)) * 
                self.speed * np.power(self.work_stage, -0.75))

        log_ss_diff = np.log10(specific_speed/2.0)
        efficiency = 0.85 - ((log_ss_diff**2) * 1.614)

        return(efficiency)
    
    def efficiency_smyth(self, verbose=False):
        peak_flow_param = 2.5
        peak_work_param = 0.55

        peak_flow_param *= (self.A_in / (self.R_mean_in**2)) / 3.2

        flow_param = self.gasflow_in.mass_flow / (
            self.gasflow_in.density() * self.R_mean_in**3 * self.speed_rad)
        work_param = self.work_stage / np.power(self.R_mean_in * self.speed_rad, 2)

        # Axis 1 is approximately horizontal, with slight downward slope
        # Axis two is approximately vertical, with slight forward slope
        axis_1_distance = (flow_param - peak_flow_param)*np.cos(np.pi/6) - (work_param - peak_work_param)*np.sin(np.pi/6)
        axis_2_distance = (flow_param - peak_flow_param)*np.sin(np.pi/6) + (work_param - peak_work_param)*np.cos(np.pi/6)

        # Approximate distance between 2% efficiency contours in Smyth chart
        contour_distance = 0.15
        semimajor_axis_1 = 0.8
        semimajor_axis_2 = 1.8

        normalised_distance = np.sqrt(np.power(axis_1_distance/semimajor_axis_1, 2) 
            + np.power(axis_2_distance/semimajor_axis_2, 2))

        if normalised_distance <= 2:
            efficiency = 0.92
        else:
            efficiency = 0.92 - 0.02*((normalised_distance-1)/contour_distance)
        
        if verbose:
            print(efficiency, peak_flow_param, flow_param, peak_work_param, work_param)

        return(efficiency)

    def estimate_efficiency(self, verbose=False):
        
        impeller_efficiency_estimate = self.efficiency_smyth(verbose)

        if impeller_efficiency_estimate<0.01:
            impeller_efficiency_estimate = 0.01
        
        stator_efficiency_estimate = impeller_efficiency_estimate

        self.efficiency = impeller_efficiency_estimate * stator_efficiency_estimate
        
        return(self.efficiency)

    def gasflow_out(self):
        delta_T = self.work_stage / self.gasflow_in.cp
        pressure_ratio = self.gasflow_in.pressure_ratio(delta_T, self.estimate_efficiency())
        T_out = self.gasflow_in.temperature + delta_T
        p_out = self.gasflow_in.pressure * pressure_ratio

        gasflow = GasFlow(self.gasflow_in.mass_flow, T_out, p_out)
        return(gasflow)

    def __repr__(self):
        out_str = ""
        out_str += " \t An axial stage, with entry pressure of {:.0f}Pa".format(self.gasflow_in.pressure)
        out_str += "\n"
        out_str += "Inlet mean radius: {:.0f}mm with bladeheight of {:.2f}mm and speed of {:.0f}RPM".format(
            self.R_mean_in*1000, self.bladeheight_in*1000, self.speed)
        out_str += "\n"
        out_str += "Total installed weight of {:.2f}kg".format(self.weight)
        out_str += "\n"
        out_str += "Equivalent system mass of {:.1f}kg-eq".format(self.ESM)
        out_str += "\n"
        out_str += "Power consumption of {:.1f}W".format(self.power)
        out_str += "\n"
        out_str += "Estimated efficiency of {:.1f}%".format(self.estimate_efficiency()*100)
        return(out_str)

class RadialStage(Component):
    def __init__(self, gasflow_in, pressure_ratio_stage, speed, work_coeff, flow_coeff, efficiency_guess, radius_hub_inlet, diffusion_ratio=3):
        self.gasflow_in = gasflow_in

        self.pressure_ratio_stage = pressure_ratio_stage

        self.work_stage = self.gasflow_in.delta_h_PR(self.pressure_ratio_stage) / efficiency_guess

        self.speed = speed
        self.speed_rad = speed * np.pi / 30
        self.flow_coeff = flow_coeff
        self.work_coeff = work_coeff

        self.efficiency = efficiency_guess
        
        self.R_hub_inlet = radius_hub_inlet

        self.diffusion_ratio = diffusion_ratio

        self.calculate_geometry_in()
        self.calculate_geometry_rotor_out()
        self.calculate_geometry_stator()
    
    def calculate_geometry_in(self):

        self.R_tip_inlet = np.sqrt(
            (np.power(2, 1/3) * np.power(self.gasflow_in.mass_flow / (self.gasflow_in.density()*np.pi), 2/3) / np.power(self.speed_rad, 2/3)) + 
            np.power(self.R_hub_inlet, 2)
        )

        self.Vx_in = self.gasflow_in.mass_flow / (
            self.gasflow_in.density()*np.pi * ((self.R_tip_inlet**2) - (self.R_hub_inlet**2)))

        self.R_mean_inlet = (self.R_tip_inlet + self.R_hub_inlet) / 2
        self.bladeheight_in = self.R_tip_inlet - self.R_hub_inlet

        self.A_in = 2 * np.pi * self.R_mean_inlet * self.bladeheight_in
    
    def calculate_geometry_rotor_out(self):
        # Using backsweep angle to find blade speed for the desired enthalpy change in the rotor
        enthalpy_rotor = self.work_stage
        slip = 0.85

        #blade_speed_exit = np.sqrt(enthalpy_rotor / (slip * (1 + (self.flow_coeff * np.tan(self.backsweep_angle_rad)))))
        blade_speed_exit = np.sqrt(enthalpy_rotor / self.work_coeff)
        self.backsweep_angle_rad = np.arctan((self.work_coeff/slip - 1) / self.flow_coeff)
        self.backsweep_angle = np.rad2deg(self.backsweep_angle_rad)
        
        R_mean_impeller_exit = blade_speed_exit / self.speed_rad
        
        self.R_mean_imp_exit = R_mean_impeller_exit

        self.V_r_imp_exit = blade_speed_exit * self.flow_coeff
        self.V_theta_imp_exit = slip*(blade_speed_exit + self.V_r_imp_exit*np.tan(self.backsweep_angle_rad))

        self.rotor_exit_flow_angle_rad = np.arctan(self.V_theta_imp_exit/self.V_r_imp_exit)
        self.rotor_exit_flow_angle = np.rad2deg(self.rotor_exit_flow_angle_rad)

        self.V_imp_exit = np.sqrt(self.V_r_imp_exit**2 + self.V_theta_imp_exit**2)
        self.A_imp_exit = self.gasflow_in.mass_flow / (self.V_imp_exit * self.gasflow_in.density())
        self.bladeheight_imp_exit = self.A_imp_exit / (2 * np.pi * self.R_mean_imp_exit)

        # Weisner correlation
        self.N_blades_rotor = np.round(np.power(np.sqrt(np.cos(self.backsweep_angle_rad))/(1-slip), 1/0.7))

    def calculate_geometry_stator(self):
        # Choose a length and area ratio along the optimum efficiency line
        stator_flow_angle_rad = abs(self.rotor_exit_flow_angle_rad)

        self.R_stator_in = self.R_mean_imp_exit * 1.1
        self.A_stator_in = 2 * np.pi * self.R_stator_in * self.bladeheight_imp_exit * np.cos(stator_flow_angle_rad)
        
        self.V_r_stator_in = self.V_r_imp_exit / 1.1
        self.V_theta_stator_in = self.V_theta_imp_exit / 1.1
        self.V_stator_in = self.V_imp_exit / 1.1

        self.A_stator_exit = self.A_stator_in * self.diffusion_ratio
        self.R_stator_exit = self.A_stator_exit / (2 * np.pi * np.cos(stator_flow_angle_rad) * self.bladeheight_imp_exit)
        self.V_stator_exit = self.V_stator_in / self.diffusion_ratio

        # Using a fit from the diffuser performance map provided

        L_W_diffuser = 0.6643 * np.power(self.diffusion_ratio, 2.7023) 
    
        # Using properties for a logarithmic spiral
        b = np.cos(stator_flow_angle_rad) / np.sin(stator_flow_angle_rad)
        a = self.R_stator_in
        L_diffuser = a * np.sqrt(1 + np.power(b,2)) / b
        self.diffuser_blade_length = L_diffuser
        W_diffuser = L_diffuser / L_W_diffuser
        self.N_blades_stator = np.round((2 * np.pi * self.R_stator_in) / W_diffuser)

    def weight_estimate(self, yield_stress, density):
        revolved_radius_tip = self.R_mean_imp_exit - self.R_tip_inlet
        revolved_radius_hub = revolved_radius_tip + self.bladeheight_imp_exit
        revolved_tip_centroid = (1 - 4/(3*np.pi)) * revolved_radius_tip + self.R_tip_inlet
        #revolved_hub_centroid = (1 - 4/(3*np.pi)) * revolved_radius_hub + radius_inlet_hub

        # Uses Pappus's Centroid Theorem to get surface area tip
        surface_area_tip = 2 * np.pi * revolved_tip_centroid * (np.pi/2) * revolved_radius_tip

        # Find minimum thickness to support a hoop stress based on pressure and the 
        # arithmetic mean of inlet and outlet tip radii

        thickness_minimum = self.gasflow_in.pressure * (self.R_hub_inlet + self.R_tip_inlet) / (2*yield_stress)

        # Casing thickness is same but with safety factor 2, 5mm maximum thickness
        casing_thickness = max(5e-3, thickness_minimum*2)

        casing_weight = ((surface_area_tip * casing_thickness)) *  density

        # Also by Pappus centroid theorem
        hub_section_area = (revolved_radius_hub * (revolved_radius_tip+self.R_tip_inlet)) - (0.25 * np.pi * revolved_radius_hub**2)

        revolved_hub_centroid = ((revolved_radius_tip+self.R_tip_inlet) * revolved_radius_hub * 0.5 * (revolved_radius_tip+self.R_tip_inlet)) 

        revolved_hub_centroid -= (self.R_hub_inlet + (revolved_radius_hub *(1 - 4/(3*np.pi)))) * (0.25 * np.pi * revolved_radius_hub**2)
        revolved_hub_centroid /= hub_section_area

        hub_volume = (2 * np.pi * revolved_hub_centroid) * hub_section_area
        hub_weight = hub_volume*density

        # Add factor for diffuser plane
        # Assume radius increases by 2
        diffuser_plate_area = (np.power(self.R_stator_exit,2) - np.power(self.R_stator_in, 2)) * np.pi
        diffuser_blade_area = self.N_blades_stator * self.diffuser_blade_length * self.bladeheight_imp_exit

        diffuser_weight = (diffuser_plate_area * 2 * casing_thickness) + (1e-3 * diffuser_blade_area)
        diffuser_weight *= density

        # Add fixings approximation
        compressor_weight = (casing_weight + diffuser_weight + hub_weight) * 1.18

        if self.R_mean_imp_exit < self.R_mean_inlet:
            compressor_weight = 1e5
        
        self.rotor_weight = (casing_weight + hub_weight) * 1.18
        self.stator_weight = diffuser_weight

        self.weight = compressor_weight
        return(self.weight)
    
    def ESM_estimate(self):
        self.power = self.work_stage * self.gasflow_in.mass_flow / self.efficiency
        power_weight = self.power * .121

        self.ESM = self.weight + power_weight

        return(self.ESM)

    def efficiency_cordier(self):
        # Estimate efficiency based on the Cordier plot
        # Assume peak efficiency of 85% at specific speed=2.0
        # Drop-off of 10% by a change of 1.5 in either direction
        specific_speed = (np.sqrt(self.gasflow_in.mass_flow / self.gasflow_in.density_deltah(self.work_stage)) * 
                self.speed * np.power(self.work_stage, -0.75))

        log_ss_diff = np.log10(specific_speed/0.5)
        efficiency = 0.78 - ((log_ss_diff**2) * 1.614)

        return(efficiency)
    
    def efficiency_smyth(self, verbose):
        peak_flow_param = 2.5
        peak_work_param = 5.5

        # Smyth constraint one - assume A/r_1^2 = 5 for radials
        peak_flow_param *= (self.A_in / (self.R_mean_inlet**2)) / 5

        # Smyth constraint two - assume A2/A1 = 0.85 for radials
        peak_work_param *= (self.A_imp_exit / self.A_in) / 0.85

        # Smyth constraint three - assume r2/r1 = 2.2 for radials
        peak_flow_param *= (self.R_mean_imp_exit/self.R_mean_inlet) / 2.2


        flow_param = self.gasflow_in.mass_flow / (
            self.gasflow_in.density() * self.R_mean_inlet**3 * self.speed_rad)
        work_param = self.work_stage / np.power(self.R_mean_inlet * self.speed_rad, 2)

        if verbose:
            print("flow", flow_param, peak_flow_param)
            print("work", work_param, peak_work_param)

        # Axis 1 is approximately horizontal, with slight downward slope
        # Axis two is approximately vertical, with slight forward slope
        axis_1_distance = (flow_param - peak_flow_param)*np.cos(np.pi/6) - (work_param - peak_work_param)*np.sin(np.pi/6)
        axis_2_distance = (flow_param - peak_flow_param)*np.sin(np.pi/6) + (work_param - peak_work_param)*np.cos(np.pi/6)

        # Approximate distance between 2% efficiency contours in Smyth chart
        contour_distance = 0.3
        semimajor_axis_1 = 3
        semimajor_axis_2 = 2

        normalised_distance = np.sqrt(np.power(axis_1_distance/semimajor_axis_1, 2) 
            + np.power(axis_2_distance/semimajor_axis_2, 2))

        if normalised_distance <= 2:
            efficiency = 0.94
        else:
            efficiency = 0.94 - 0.02*((normalised_distance-1)/contour_distance)

        return(efficiency)

    def estimate_efficiency(self, verbose=False):
        
        impeller_efficiency_estimate = self.efficiency_smyth(verbose)

        if impeller_efficiency_estimate<0.01:
            impeller_efficiency_estimate = 0.01
        if self.work_coeff <= 0:
            impeller_efficiency_estimate = 0.001
        if self.flow_coeff <= 0:
            impeller_efficiency_estimate = 0.001
        
        # Correlation from diffuser performance map
        stator_efficiency = 1 - (0.7328 * np.power(self.diffusion_ratio, -0.997))

        self.efficiency = impeller_efficiency_estimate * stator_efficiency
        self.efficiency_impeller = impeller_efficiency_estimate
        self.efficiency_stator = stator_efficiency
        
        return(self.efficiency)

    def gasflow_out(self):
        delta_T = self.work_stage / self.gasflow_in.cp
        pressure_ratio = self.gasflow_in.pressure_ratio(delta_T, self.estimate_efficiency())
        T_out = self.gasflow_in.temperature + delta_T
        p_out = self.gasflow_in.pressure * pressure_ratio

        gasflow = GasFlow(self.gasflow_in.mass_flow, T_out, p_out)
        return(gasflow)

    def __repr__(self):
        out_str = ""
        out_str += " \t A radial stage, with entry pressure of {:.0f}Pa".format(self.gasflow_in.pressure)
        out_str += "\n"
        out_str += "Inlet mean radius: {:.1f}mm, inlet hub radius of {:.1f}mm and bladeheight of {:.2f}mm".format(
            self.R_mean_inlet*1000, self.R_hub_inlet*1000, self.bladeheight_in*1000)
        out_str += "\n"
        out_str += "Impeller outlet radius: {:.1f}mm, blade height of {:.2f}mm and speed of {:.0f}RPM".format(
            self.R_mean_imp_exit*1000, self.bladeheight_imp_exit*1000, self.speed)
        out_str += "\n"
        out_str += "Diffuser outlet radius: {:.1f}mm".format(
            self.R_stator_exit*1000)
        out_str += "\n"
        out_str += "Total installed weight of {:.2f}kg".format(self.weight)
        out_str += "\n"
        out_str += "Equivalent system mass of {:.1f}kg-eq".format(self.ESM)
        out_str += "\n"
        out_str += "Power consumption of {:.1f}W".format(self.power)
        out_str += "\n"
        out_str += "Estimated efficiency of {:.1f}%".format(self.estimate_efficiency()*100)
        return(out_str)
    

class HeatExchanger(Component):
    def __init__(self, gasflow_in, delta_h, air_speed, pressure_drop):
        self.gasflow_in = gasflow_in
        self.delta_h = delta_h
        self.air_speed = air_speed
        self.pressure_drop = pressure_drop

        self.REL_values = False
    
    def weight_estimate(self):
        # Using simplified methodology from https://intech-gmbh.com/heat_exchangers_calc_and_select/
        # Heat exchanger tubes assumed to be 5mm tubes with 0.5mm wall thickness, made of copper
        # Good design - delta T assumed to be 10 degrees throughout
        # More advanced model may deploy https://arc.aiaa.org/doi/pdf/10.2514/6.2021-3711
        # Assume perfect heat exchange on liquid side, Nusselt number on air side scales as Re^0.8

        if self.REL_values:
            diameter = 2e-3
            wall_thickness = 4e-5
        else:
            diameter = 5e-3
            wall_thickness = 0.5e-3

        cooling_power = self.gasflow_in.mass_flow * self.delta_h
        Re = self.air_speed * self.gasflow_in.density() * diameter / 14.4e-6
        Nu_air = 0.029 * np.power(Re, 0.8) * 0.9244
        h = Nu_air * 0.011 /  diameter
        delta_T = 10
        area = cooling_power / (h * delta_T)

        hx_weight = area * wall_thickness * 8960

        self.weight = hx_weight
        return(hx_weight)

class HeatExchangerAdvanced(Component):
    def __init__(self, cooling_power, airflow, inlet_area, bulk_area_ratio, 
            tube_diameter=2e-3, pitch_diameter_ratio=1.25, 
            wall_thickness=40e-6, wall_thermal_conductivity=300, geometry="wall",
            coolant="VLT", coolant_velocity=2, neglect_peak_velocity=False):
        self.cooling_power = cooling_power
        self.airflow = airflow

        self.bulk_area_ratio = bulk_area_ratio
        self.inlet_area = inlet_area
        self.bulk_area = self.inlet_area * bulk_area_ratio
        self.tube_diameter = tube_diameter
        self.pitch_diameter_ratio = pitch_diameter_ratio
        self.geometry = geometry
        # geometry is "row", "staggered" or "wall"
        # row is distinct tubes with gaps in a grid
        # staggered is distinct tubes with gaps, staggered
        # wall is tubes sintered into a plate, with roughness height equal to 30% of tube diameter
        # Pitch diameter ratio is the crosswise/transverse pitch
        
        self.wall_thickness = wall_thickness
        self.wall_thermal_conductivity = wall_thermal_conductivity

        self.coolant_velocity = coolant_velocity
        if coolant == "glycol":
            # Coolant in fluid side is ethylene glycol
            # Data from https://wiki.anton-paar.com/en/automotive-antifreeze/
            # and https://www.researchgate.net/publication/230074410_Thermal_Conductivity_Density_Viscosity_and_Prandtl-Numbers_of_Ethylene_Glycol-Water_Mixtures
            # VW CG13 at 60:40 ratio with water, freezes at -50C
            # Properties evaluated at -30C, viscosity-weighted average
            self.coolant_inlet = 228
            self.coolant_viscosity = 0.107
            self.coolant_density = 1120
            self.coolant_dT = 20
            self.coolant_thermal_conductivity = 0.32
            self.coolant_Pr = 600
        elif coolant == "VLT":
            # Coolant in fluid side is Therminol VLT
            # Data from https://www.therminol.com/sites/therminol/files/documents/TF-21_Therminol_VLT.pdf
            # Properties evaluated at -45C
            self.coolant_inlet = 200
            self.coolant_viscosity = 2.16e-3
            self.coolant_density = 804
            self.coolant_dT = 50
            self.coolant_thermal_conductivity = 0.1169
            self.coolant_Pr = 31.04
        
        elif coolant == "supercrit He":
            # Coolant in fluid side is supercritical helium at 20 bar pressure
            # Data from NIST databook
            # With constant delta T of 50C (250K to 330K), evaluate properties at 240K
            self.coolant_inlet = 200
            self.coolant_viscosity = 1.7231e-5
            self.coolant_density = 3.9639
            self.coolant_dT = 50
            self.coolant_thermal_conductivity = 0.13517
            self.coolant_Pr = 0.6623
        
        elif coolant == "ammonia":
            # Coolant in fluid side is liquid ammonia at 20 bar pressure
            # Data from NIST databook
            # With constant delta T of 50C (250K to 330K), evaluate properties at 240K
            self.coolant_inlet = 200
            self.coolant_viscosity = 0.00025696
            self.coolant_density = 683.88
            self.coolant_dT = 50
            self.coolant_thermal_conductivity = 0.58261
            self.coolant_Pr = 1.9613
        
        elif coolant == "concept":
            # Coolant in fluid side is liquid ammonia at 20 bar pressure
            # Data from NIST databook
            # With constant delta T of 50C (250K to 330K), evaluate properties at 240K
            self.coolant_inlet = 200
            self.coolant_viscosity = 1.0e-4
            self.coolant_density = 667
            self.coolant_dT = 50
            self.coolant_thermal_conductivity = 0.1
            self.coolant_Pr = 1.5
        
        self.neglect_peak_velocity = neglect_peak_velocity
        
        self.make_geometry()
    
    def make_geometry(self):
        self.duct_size = np.sqrt(self.bulk_area)
        self.inlet_size = np.sqrt(self.inlet_area)
        # Geometry of inlet diffuser
        self.diffusion_ratio = self.bulk_area / self.inlet_area

        self.diffuser_length = 0.6643 * np.power(self.diffusion_ratio, 2.7023/2) * self.inlet_size

        # Geometry of outlet nozzle
        nozzle_angle = 30
        self.nozzle_length = (self.diffusion_ratio - 1) / (4 * np.sin(np.deg2rad(nozzle_angle))) * self.inlet_size

        # Geometry of actual HX duct
        self.heat_transfer_per_unit_area = self.heat_transfer_row()
        self.required_area = self.cooling_power / self.heat_transfer_per_unit_area

        
        self.n_tubes_wide = np.floor(self.duct_size / (self.pitch_diameter_ratio * self.tube_diameter)) - 1
        if self.n_tubes_wide < 1:
            self.n_tubes_wide = 1

        if self.geometry == "wall":
            self.duct_length = self.required_area / (2 * self.n_tubes_wide * self.duct_size)
            self.n_rows = np.ceil(self.duct_length / self.tube_diameter)
            self.n_tubes = self.n_rows * self.n_tubes_wide
        else:
            area_per_row = self.duct_size * np.pi * self.tube_diameter * self.n_tubes_wide
            self.n_rows = np.ceil(self.required_area / area_per_row)

            self.n_tubes = self.n_rows * self.n_tubes_wide
            self.duct_length = self.n_tubes * (self.pitch_diameter_ratio * self.tube_diameter)
    
    def weight_estimate(self):
        # Assume a tube density of 9000
        one_tube_weight = self.duct_size * self.wall_thickness * np.pi * self.tube_diameter * 9000
        tubes_weight = one_tube_weight * self.n_tubes

        self.tubes_weight = tubes_weight

        # Assume a duct material density of 2700 (aluminium)
        # and safety factor of 3 on hoop stress, or 1mm thick
        duct_wall_area = 4 * self.duct_length * self.duct_size
        duct_wall_thickness = max(self.airflow.pressure * (self.duct_size/2) / 133e6, 1e-3)
        duct_weight = duct_wall_area * duct_wall_thickness * 2700

        # Diffuser wall area
        diffuser_wall_area = 0.5 * (self.duct_size + self.inlet_size) * self.diffuser_length
        diffuser_wall_area *= 4 
        # 4 walls in a diffuser
        diffuser_weight = diffuser_wall_area * duct_wall_thickness * 2700

        # Nozzle wall area
        nozzle_wall_area = 0.5 * (self.duct_size + self.inlet_size) * self.nozzle_length
        nozzle_wall_area *= 4 
        # 4 walls in a nozzle
        nozzle_weight = nozzle_wall_area * duct_wall_thickness * 2700

        self.weight = duct_weight + tubes_weight + diffuser_weight + nozzle_weight

    def heat_transfer_row(self):
        # Heat transfer inside the tube
        D_internal = self.tube_diameter - (2 * self.wall_thickness)
        Re_coolant = self.coolant_velocity * self.coolant_density * D_internal / self.coolant_viscosity

        if Re_coolant < 2300:
            # Coolant flow is laminar
            Nu_coolant = 3.66
        else:
            # Coolant flow is turbulent
            # Unlikely for glycol
            # Dottis-Boelter correlation
            Nu_coolant = 0.023 * np.power(Re_coolant, 0.8) * np.power(self.coolant_Pr, 1/3)
        
        heat_transfer_coefficient_coolant = Nu_coolant * self.coolant_thermal_conductivity / D_internal
        coolant_thermal_resistance = 1 / (heat_transfer_coefficient_coolant * np.pi * D_internal)

        # Wall thermal resistance
        wall_thermal_resistance = np.log(self.tube_diameter / D_internal) / (2 * np.pi * self.wall_thermal_conductivity)

        
        # Gas thermal resistance
        bulk_velocity = self.airflow.mass_flow / (self.bulk_area * self.airflow.density())
        if self.neglect_peak_velocity:
            peak_velocity = bulk_velocity
        else:
            peak_velocity = bulk_velocity * self.pitch_diameter_ratio / (self.pitch_diameter_ratio - 1)

        rho_tubes = self.airflow.density_static(peak_velocity)
        if self.geometry == "row" or self.geometry == "staggered":
            # Implements https://thermopedia.com/cn/content/1212/
            
            Re_tubes = peak_velocity * rho_tubes * self.tube_diameter / self.airflow.viscosity

            c = 1
            m = 1
            if self.geometry == "staggered":
                if Re_tubes < 3e2:
                    c = 1.309
                    m = 0.360
                elif Re_tubes < 2e5:
                    c = 0.273
                    m = 0.635
                elif Re_tubes < 2e6:
                    c = 0.124
                    m = 0.700
            else:
                if Re_tubes < 3e2:
                    c = 0.742
                    m = 0.431
                elif Re_tubes < 2e5:
                    c = 0.211
                    m = 0.651
                elif Re_tubes < 2e6:
                    c = 0.116
                    m = 0.700
            
            Nu_gas = c * np.power(Re_tubes, m) * np.power(self.airflow.Pr, 1/3)
            heat_transfer_coefficient_gas = Nu_gas * self.airflow.k / self.tube_diameter
            gas_thermal_resistance = 1 / (heat_transfer_coefficient_gas * np.pi * self.tube_diameter)
        
        elif self.geometry == "wall":
            D_hydraulic = 2 * (self.pitch_diameter_ratio - 1) * self.tube_diameter
            Re_tubes = peak_velocity * rho_tubes * D_hydraulic / self.airflow.viscosity

            f = self.friction_factor_walls()
                
            if Re_tubes < 2300:
                # Cengel pg 437
                Nu_gas = 8.24
            
            elif Re_tubes < 3e3:
                # Gnielinski relationship, Cengel pg 441

                Nu_gas = (f/8) * (Re_tubes - 1000) * self.airflow.Pr / (
                    1 + 12.7 * np.power(f/8, 0.5) * (np.power(self.airflow.Pr, 2/3) - 1)
                )

            elif 3e3 <= Re_tubes and Re_tubes< 5e6:
                # Gnielinski relationship, Cengel pg 441
                Nu_gas = (f/8) * (Re_tubes - 1000) * self.airflow.Pr / (
                    1 + 12.7 * np.power(f/8, 0.5) * (np.power(self.airflow.Pr, 2/3) - 1)
                )
            else:
                # Chilton-Colburn analogy
                Nu_gas = 0.125 * f * Re_tubes * np.power(self.airflow.Pr, 1/3)

            heat_transfer_coefficient_gas = Nu_gas * self.airflow.k / D_hydraulic
            gas_thermal_resistance = 1 / (heat_transfer_coefficient_gas * np.pi * self.tube_diameter)

        thermal_resistance = coolant_thermal_resistance + wall_thermal_resistance + gas_thermal_resistance

        heat_flow_per_unit_area = self.coolant_dT / thermal_resistance

        return(heat_flow_per_unit_area)

    def gas_pressure_drop(self):
        bulk_velocity = self.airflow.mass_flow / (self.bulk_area * self.airflow.density())   
        if self.neglect_peak_velocity:
            peak_velocity = bulk_velocity
        else:
            peak_velocity = bulk_velocity * self.pitch_diameter_ratio / (self.pitch_diameter_ratio - 1)

        rho_tubes = self.airflow.density_static(peak_velocity)

        if self.geometry == "row" or self.geometry == "staggered":

            Re_tubes = peak_velocity * rho_tubes * self.tube_diameter / self.airflow.viscosity
            # Implements equations from https://thermopedia.com/cn/content/1212/
            
            # Euler number correlations
            c0 = 0
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0

            if self.geometry == "staggered":
                # Use staggered heat exchanger correlations
                if self.pitch_diameter_ratio <= 1.375:
                    # Use correlation for a=1.25
                    if Re_tubes < 1e3:
                        c0 =  0.795
                        c1 =  0.247e3
                        c2 =  0.335e3
                        c3 = -0.155e4
                        c4 =  0.241e4
                    elif Re_tubes >= 1e3 and Re_tubes < 2e6:
                        c0 =  0.245 
                        c1 =  0.339e4
                        c2 = -0.984e7
                        c3 =  0.132e11
                        c4 = -0.599e13
                elif self.pitch_diameter_ratio <= 1.75:
                    # Use correlation for a=1.5
                    if Re_tubes < 1e3:
                        c0 =  0.683
                        c1 =  0.111e2
                        c2 = -0.973e2
                        c3 = -0.426e3
                        c4 =  0.574e3
                    elif Re_tubes >= 1e-3 and Re_tubes < 2e6:
                        c0 =  0.203
                        c1 =  0.248e4
                        c2 = -0.758e7
                        c3 =  0.104e11
                        c4 = -0.482e13
                elif self.pitch_diameter_ratio <= 2.25:
                    # Use correlation for a=2
                    if Re_tubes < 1e2:
                        c0 =  0.713
                        c1 =  0.448e2
                        c2 = -0.126e3
                        c3 = -0.582e3
                        c4 =  0
                    elif Re_tubes >= 1e2 and Re_tubes < 1e4:
                        c0 =  0.343 
                        c1 =  0.303e3
                        c2 = -0.717e5
                        c3 =  0.880e7
                        c4 = -0.380e9
                    elif Re_tubes >= 1e4 and Re_tubes < 1e6:
                        c0 =  0.162 
                        c1 =  0.181e4
                        c2 =  0.792e8
                        c3 = -0.165e13
                        c4 =  0.872e16
                elif self.pitch_diameter_ratio <= 2.75:
                    # Use correlation for a=2.5
                    if Re_tubes < 5e3:
                        c0 =  0.330
                        c1 =  0.989e2
                        c2 = -0.148e5
                        c3 =  0.192e7
                        c4 = -0.862e8
                    elif Re_tubes >= 5e3 and Re_tubes < 2e6:
                        c0 =  0.119 
                        c1 =  0.498e4
                        c2 = -0.507e8
                        c3 =  0.251e12
                        c4 = -0.463e15
            else:
                # Use inline square bank correlations
                if self.pitch_diameter_ratio <= 1.375:
                    # Use correlation for a=1.25
                    if Re_tubes < 2e3:
                        c0 =  0.272
                        c1 =  0.207e3
                        c2 =  0.103e3
                        c3 = -0.286e3
                        c4 =  0
                    elif Re_tubes >= 2e-3 and Re_tubes < 2e6:
                        c0 =  0.267 
                        c1 =  0.249e3
                        c2 = -0.927e7
                        c3 =  0.100e11
                        c4 =  0
                elif self.pitch_diameter_ratio <= 1.75:
                    # Use correlation for a=1.5
                    if Re_tubes < 2e3:
                        c0 =  0.263
                        c1 =  0.867e2
                        c2 = -0.202e0
                        c3 =  0
                        c4 =  0
                    elif Re_tubes >= 2e-3 and Re_tubes < 2e6:
                        c0 =  0.235 
                        c1 =  0.197e4
                        c2 = -0.124e8
                        c3 =  0.312e11
                        c4 = -0.274e14
                elif self.pitch_diameter_ratio <= 2.25:
                    # Use correlation for a=2
                    if Re_tubes < 800:
                        c0 =  0.188
                        c1 =  0.566e2
                        c2 = -0.646e3
                        c3 =  0.601e4
                        c4 = -0.183e5
                    elif Re_tubes >= 800 and Re_tubes < 2e6:
                        c0 =  0.247 
                        c1 = -0.595
                        c2 =  0.150
                        c3 = -0.137
                        c4 =  0.396

            eul =   (c0 + 
                    c1 * np.power(Re_tubes, -1) + 
                    c2 * np.power(Re_tubes, -2) + 
                    c3 * np.power(Re_tubes, -3) +
                    c4 * np.power(Re_tubes, -4))


            pressure_drop_row = eul * (0.5 * rho_tubes * peak_velocity**2)

            self.pressure_drop_hx = pressure_drop_row * self.n_rows

        elif self.geometry == "wall":
            D_hydraulic = 2 * (self.pitch_diameter_ratio - 1) * self.tube_diameter

            f = self.friction_factor_walls()

            self.pressure_drop_hx = (self.duct_length / D_hydraulic) * f * (0.5 * rho_tubes * peak_velocity**2)
    
        # Pressure drop from diffusers

        # Based on performance map interpolation, loss coefficient is about 0.15 for all ratios
        # 5% loss in nozzle

        inlet_velocity = self.airflow.mass_flow / (self.inlet_area * self.airflow.density())
        self.pressure_drop_diffuser = (0.5 * self.airflow.density() * np.power(inlet_velocity, 2)) * (0.15 + 0.05)

        self.pressure_drop = self.pressure_drop_hx + self.pressure_drop_diffuser


        # TODO Pressure drop from bends
        # Implements https://www.thermopedia.com/content/577/

    def coolant_pumping_power(self):
        # Find the power required to pump the coolant through the tubes
        D_internal = self.tube_diameter - (2 * self.wall_thickness)
        Re_coolant = self.coolant_velocity * self.coolant_density * D_internal / self.coolant_viscosity

        if Re_coolant < 2300:
            # Coolant flow is laminar
            f = 64 / Re_coolant
        else:
            # Coolant flow is turbulent
            # Unlikely for glycol
            # Haaland correlation for f
            roughness = 0.002e-3
            a = -1.8 * np.log(6.9/Re_coolant + np.power(roughness/D_internal / 3.7, 1.11))
            f = 1 / np.power(a, 2)
        
        pressure_drop_one_tube = f * (self.duct_size / D_internal) * self.coolant_density * self.coolant_velocity**2 / 2
        pressure_drop_exchanger = pressure_drop_one_tube * self.n_rows

        pumping_power = pressure_drop_exchanger * self.coolant_velocity * np.pi * (D_internal**2 / 4) * self.n_tubes_wide

        self.pumping_power = pumping_power

    def friction_factor_walls(self):
        bulk_velocity = self.airflow.mass_flow / (self.bulk_area * self.airflow.density())
        peak_velocity = bulk_velocity * self.pitch_diameter_ratio / (self.pitch_diameter_ratio - 1)
        rho_tubes = self.airflow.density_static(peak_velocity)
        D_hydraulic = 2 * (self.pitch_diameter_ratio - 1) * self.tube_diameter
        Re_tubes = peak_velocity * rho_tubes * D_hydraulic / self.airflow.viscosity

        if Re_tubes < 2100:
            # Laminar flow
            f = 16 / Re_tubes
        
        elif Re_tubes < 1e4:
            # Transition region: use greater of the two
            f_lam = 16 / Re_tubes

            # Haaland approximation
            rough = 0.15 * self.tube_diameter / D_hydraulic
            a = -1.8 * np.log(6.9/Re_tubes + np.power(rough / 3.7, 1.11))
            f_turb = 1 / np.power(a, 2)

            f = max(f_lam, f_turb)
        
        elif 1e4 <= Re_tubes and Re_tubes < 4e8:
            rough = 0.3 * self.tube_diameter / D_hydraulic
            a = -1.8 * np.log(6.9/Re_tubes + np.power(rough / 3.7, 1.11))
            f = 1 / np.power(a, 2)
        else:
            #print(Re_tubes, bulk_velocity, self.bulk_area)
            rough = 0.3 * self.tube_diameter / D_hydraulic
            a = -1.8 * np.log(6.9/Re_tubes + np.power(rough / 3.7, 1.11))
            f = 1 / np.power(a, 2)
        
        return(f)

    def estimate_ESM(self):
        self.power_weight = self.pumping_power * .121

        self.ESM = self.weight + self.power_weight

    def __repr__(self):
        out_str = ""
        relative_pd = 100 * self.pressure_drop / self.airflow.pressure
        if self.geometry == "wall":
            out_str += " \t A sintered tube-wall heat exchanger, with entry pressure of {:.0f}Pa and a drop of {:.1f}%".format(
                self.airflow.pressure, relative_pd)
        elif self.geometry == "staggered":
            out_str += " \t A staggered tube heat exchanger, with entry pressure of {:.0f}Pa and a drop of {:.1f}%".format(
                self.airflow.pressure, relative_pd)
        elif self.geometry == "row":
            out_str += " \t A tube array heat exchanger, with entry pressure of {:.0f}Pa and a drop of {:.1f}%".format(
                self.airflow.pressure, relative_pd)
        out_str += "\n"
        out_str += "Cooling power: {:.0f}W with temperature drop of {:.1f}K ".format(
            self.cooling_power, self.cooling_power / (self.airflow.mass_flow * self.airflow.cp))
        out_str += "\n"
        out_str += "Total number of tube passes: {:.0f}".format(self.n_tubes)
        out_str += "\n"
        out_str += "Total installed length of {:.1f}mm, inlet width of {:.1f}mm and duct width of {:.1f}mm".format(
            (self.duct_length + self.diffuser_length + self.nozzle_length)*1000, self.inlet_size*1000, self.duct_size*1000)
        out_str += "\n"
        out_str += "Total installed weight of {:.2f}kg".format(self.weight)
        out_str += "\n"
        out_str += "Coolant pumping power of {:.2f}W".format(self.pumping_power)
        out_str += "\n"
        out_str += "Equivalent system mass of {:.1f}kg-eq".format(self.ESM)

        return(out_str)
        

"""bc_flow = GasFlow(0.05, 230, 1900)

radial_1 = RadialStage(bc_flow, 11e3, 6000, 1.0, 0.75, 0.88, 0.1, 60)

print(radial_1.weight_estimate(500e6, 2700))
print(radial_1.weight_estimate(2000e6, 1700))"""
"""suitflow = GasFlow(0.0329, 300, 101325, carbon_dioxide=False)
suit_hx = HeatExchangerAdvanced(200, suitflow, 0.01131, 0.01131, 
    tube_diameter=1.6e-3, wall_thickness=0.5e-3,
    geometry="wall", pitch_diameter_ratio=16.6, coolant_velocity=0.5)
suit_hx.gas_pressure_drop()
suit_hx.coolant_pumping_power()
suit_hx.weight_estimate()
print(suit_hx.pressure_drop)
print(suit_hx.duct_length*1000)
print(suit_hx.pumping_power)
print(suit_hx.n_tubes_wide)
print(suit_hx.n_tubes)
print(suit_hx.tubes_weight)
print(suit_hx.weight - suit_hx.tubes_weight)
print(suit_hx.weight)"""
#bc_flow = GasFlow(0.05, 300, 1900)
#wall_hx = HeatExchangerAdvanced(2000, bc_flow, 0.015, 0.045, geometry="wall", pitch_diameter_ratio=2)
#wall_hx.gas_pressure_drop()
#tube_hx = HeatExchangerAdvanced(2000, bc_flow, 0.015, 0.045, geometry="row", pitch_diameter_ratio=2)
#tube_hx.gas_pressure_drop()
#print(wall_hx.pressure_drop, tube_hx.pressure_drop)
#print(wall_hx.duct_length, tube_hx.duct_length)
#print(wall_hx.n_tubes, tube_hx.n_tubes)
#print(wall_hx.n_rows, tube_hx.n_rows)
#print(wall_hx.n_tubes_wide, tube_hx.n_tubes_wide)
