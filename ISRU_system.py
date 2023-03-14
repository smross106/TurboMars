import numpy as np
import matplotlib.pyplot as plt
from component import Component, CondensingHX, BufferTank

"""
This makes full use of Lukas Schrenk's MIT masters thesis

Available at https://web.archive.org/web/20191011015910/http://lukas-schrenk.com/files/MasterThesis_Schrenk.pdf
"""

class WaterElectrolyserBasic(Component):
    def __init__(self, mass_flow, fudge_mass=4*3, fudge_power=1):
        self.mass_flow = mass_flow
        self.fudge_mass = fudge_mass
        self.fudge_power = fudge_power
    
        self.cell_voltage = 1.7
        self.current_density = 25e3
        self.r_cell = 0.09/2
        self.A_cell = np.power(self.r_cell, 2) * np.pi

        self.rho_titanium = 4000
        self.rho_membrane = 0.5

        self.power_coefficient = 45700 * 3600
    
    def estimate_power(self):
        self.power = self.mass_flow * self.power_coefficient
    
    def estimate_mass(self):
        self.N_cells = self.power / (self.cell_voltage * self.current_density * self.A_cell)
        self.m_one_cell = self.A_cell * ((4e-3*self.rho_titanium) + (2.4e-3*self.rho_titanium*0.6)+self.rho_membrane)
        self.m_cells = self.m_one_cell * self.N_cells

        self.H_stack = self.N_cells * 6.7e-3
        self.t_case = 1e-3

        self.m_case = ((2*np.power(self.r_cell+self.t_case,2)*np.pi*self.t_case) + 
                       (2*self.t_case*(self.r_cell+self.t_case)*np.pi*self.H_stack)) * self.rho_titanium

        self.weight = self.m_cells + self.m_case
    
    def estimate_ESM(self):
        self.estimate_power()
        self.estimate_mass()

        self.power *= self.fudge_power
        self.weight *= self.fudge_mass

class WaterElectrolyserScaled(Component):
    def __init__(self, water_mass_flow, h24=True):
        # Basis is IMISPPS electrolyser (33W/kg)
        # https://ascelibrary.org/doi/epdf/10.1061/%28ASCE%29AS.1943-5525.0000201
        # Also used is Lukas Shcrenk electrolyser (75W/kg)
        # and https://www.sciencedirect.com/science/article/pii/B9780444527455002884 (50-500W/kg)
        self.mass_flow = water_mass_flow
        self.h24 = h24
        self.night = False

        self.specific_power = 50 #W/kg
        self.specific_energy = 45700 * 3600 * (2/18) # J/kg H2O electrolysed

        self.estimate_ESM()

    def estimate_weight(self):
        self.weight = self.power / self.specific_power
    
    def estimate_power(self):
        self.power = self.specific_energy * self.mass_flow
    
    def estimate_ESM(self):
        self.estimate_power()
        self.estimate_weight()
        self.cooling = 0
        
class ReactorSeparatorScaled(Component):
    def __init__(self, total_output_mass_flow, h24):
        # Basis is flight-like IMISPPS reactor, heat exchangers, separators, etc 
        # https://ascelibrary.org/doi/epdf/10.1061/%28ASCE%29AS.1943-5525.0000201
        self.mass_flow = total_output_mass_flow
        self.h24 = h24
        self.night = False

        self.basis_output_mass_flow = 0.7 / (24 * 3600)
        self.basis_power = 186
        self.basis_mass = 11
        self.scaling_index_mass = 0.75

        self.estimate_ESM()
    
    def estimate_weight(self):
        self.weight = self.basis_mass * np.power(self.mass_flow/self.basis_output_mass_flow, self.scaling_index_mass)
    
    def estimate_power(self):
        self.power = self.basis_power * (self.mass_flow/self.basis_output_mass_flow)
    
    def estimate_ESM(self):
        self.estimate_weight()
        self.estimate_power()
        self.cooling = 0

class CryocoolScaled(Component):
    def __init__(self, mass_flow, h24, gas_inlet_temperature=298, liquid_outlet_temperature=111):
        # System to cool gaseous oxygen/methane stream to liquifaction point using a cryocooler
        # Basis is flight-like IMISPPS cryocooler and distillation column
        # https://ascelibrary.org/doi/epdf/10.1061/%28ASCE%29AS.1943-5525.0000201
        self.mass_flow = mass_flow
        self.h24 = h24
        self.night = False
        self.gas_inlet_temperature = gas_inlet_temperature
        self.liquid_outlet_temperature = liquid_outlet_temperature
        

        self.basis_output_mass_flow = 0.7 / (24 * 3600)
        self.basis_power = 165
        self.basis_mass = 7
        self.scaling_index_mass = 0.8
        self.basis_gas_inlet_temperature = 298
        self.basis_liquid_outlet_temperature = 111

        self.estimate_ESM()
    
    def estimate_weight(self):
        self.weight = self.basis_mass * np.power(self.mass_flow/self.basis_output_mass_flow, self.scaling_index_mass)
    
    def estimate_power(self):
        self.power = self.basis_power * (self.mass_flow/self.basis_output_mass_flow) * (
            (self.gas_inlet_temperature - self.liquid_outlet_temperature) / (self.basis_gas_inlet_temperature-self.basis_liquid_outlet_temperature))
    
    def estimate_ESM(self):
        self.estimate_weight()
        self.estimate_power()
        self.cooling = self.power

class RodwellWaterScaled(Component):
    def __init__(self, mass_flow, h24):
        # Data from 
        # https://ntrs.nasa.gov/api/citations/20180007948/downloads/20180007948.pdf
        self.mass_flow = mass_flow
        self.h24 = h24
        self.night = False

        self.basis_system_mass = 1413.7 + 150.5
        self.basis_system_flow = 379 / (24 * 3600)

        self.power_coefficient = 2.28e6

        self.estimate_ESM()

    def estimate_power(self):
        self.power = self.power_coefficient * self.mass_flow

    def estimate_mass(self):
        # Mass of drilling system
        self.fixed_mass = 3207.0

        # Scale the baseline mass
        self.scaled_mass = np.ceil(self.mass_flow/self.basis_system_flow) * self.basis_system_mass

        # Assume 3km pipeline to base
        volume_flow = self.mass_flow / 1000
        pipe_radius = np.sqrt(volume_flow/np.pi)
        pipe_radius = max(pipe_radius, 2e-3)
        
        self.pipe_mass = 2 * np.pi * pipe_radius * 3000 * 5e-3 * 8000

        print(self.fixed_mass, self.scaled_mass, self.pipe_mass)
        
        self.weight = self.fixed_mass + self.scaled_mass + self.pipe_mass
    
    def estimate_ESM(self):
        self.estimate_power()
        self.estimate_mass()
        self.cooling = 0

class WaterExtractorSoilBasic(Component):
    def __init__(self, mass_flow, ice_concentration, batch_time, use_heat=False, fudge_mass=1.8, fudge_power=5):
        self.mass_flow = mass_flow
        self.ice_concentration = ice_concentration
        self.batch_time = batch_time
        self.use_heat = use_heat

        self.fudge_mass = fudge_mass
        self.fudge_power = fudge_power

        self.regolith_density = 1800
        self.volume_flow = 100*self.mass_flow / (self.regolith_density*self.ice_concentration)
        self.batch_volume = self.volume_flow * self.batch_time

        self.rho_aluminium = 2700

        self.estimate_ESM()

    def estimate_mass(self):
        t_material = 5e-3
        # Oven
        V_oven = 0.75 * self.mass_flow/(self.ice_concentration/100) * self.batch_time / self.regolith_density
        h_oven = np.power(2 * V_oven / (np.pi * np.power(np.tan(np.pi/12),2)), 1/3)
        r_oven = h_oven * np.tan(np.pi/12)
        m_oven = np.pi * r_oven * (r_oven + np.sqrt(np.power(h_oven,2)+np.power(r_oven,2)))
        m_oven *= t_material * self.rho_aluminium


        # Auger system
        self.d_auger1 = 0.0956 * np.log(self.volume_flow * 3600) + 0.0029
        d_shaft1 = self.d_auger1/5
        t_blade1 = self.d_auger1/15
        self.L_auger1 = h_oven*3
        m_auger1 = ((np.pi * np.power(d_shaft1/2,2) * self.L_auger1) + 
                    (np.sqrt(np.power(0.3*self.d_auger1, 2) + np.power(self.d_auger1/(2*np.pi), 2)) 
                     * 2*np.pi*self.L_auger1*t_blade1*self.L_auger1/self.d_auger1))
        m_auger1 *= self.rho_aluminium

        d_auger2 = r_oven
        d_shaft2 = d_auger2/5
        t_blade2 = d_auger2/120
        L_auger2 = h_oven
        m_auger2 = ((np.pi * np.power(d_shaft2/2,2) * L_auger2) + 
                    (np.sqrt(np.power(0.3*d_auger2, 2) + np.power(d_auger2/(16*np.pi), 2)) 
                     * 2*np.pi*L_auger2*t_blade2*L_auger2/(8*d_auger2)))
        m_auger2 *= self.rho_aluminium

        # Hopper and sieve
        V_hopper = 2 * self.batch_volume * 1.1
        h_hopper = 0.5
        a_hopper = np.sqrt(V_hopper/h_hopper)
        m_hopper = (a_hopper + 4*h_hopper)*a_hopper * t_material * self.rho_aluminium

        m_sieve = np.power(a_hopper,2) * t_material * self.rho_aluminium * 0.5

        # Sweep gas heater
        m_sweep_gas = 2.37 * self.mass_flow * 3600

        # Sweep gas heat exchanger
        l_hx = h_oven*6
        R_hx = 0.05
        m_hx = l_hx * np.pi * np.power(R_hx + t_material, 2) * t_material

        # Water filter
        m_filter = 200 * self.mass_flow * 3600 * 3.5e-4 + 1443 * 0.61 * 4 * np.power(a_hopper/100,2)

        m_pump = 0.54 * np.ceil(self.mass_flow*12)

        self.weight = m_oven + m_auger1 + m_auger2 + m_hopper + m_sieve + m_sweep_gas + m_hx + m_filter + m_pump

    def estimate_power(self):
        heating_energy = ((1000 * 100 / self.ice_concentration) * 100) + 1.5*335e3 + 4200*50 +1.5*2260e3
        P_heating = heating_energy * self.mass_flow / 0.9

        P_pump = 10 * np.ceil(self.mass_flow*12)

        # Power for augers
        rpm_auger1 = (-13.83*np.log(self.volume_flow*3600) + 138.36)/3
        S_auger1 = 2000 * np.power(self.d_auger1,2)
        P_auger1 = 10*746*(S_auger1*rpm_auger1 + 0.7*self.volume_flow*self.regolith_density) / 1e6

        print(heating_energy)
        print(P_auger1, P_heating)

        if self.use_heat:
            self.power = P_pump + P_auger1
        else:
            self.power = P_heating + P_pump + P_auger1

    def estimate_ESM(self):
        self.estimate_mass()
        self.estimate_power()

        self.power *= self.fudge_power
        self.weight *= self.fudge_mass

        self.cooling = 0


# IMMISPS
"""el = WaterElectrolyserScaled(10*0.55*1.11E-6*9)
el.estimate_ESM()

re = ReactorSeparatorScaled(10*0.7/(24*3600))
re.estimate_ESM()

methcool = CryocoolScaled(10*0.7/(24*3600))
methcool.estimate_ESM()"""

# Starship-sized
co2_flow = 1000 
water_flow = 800
ox_flow = 1440
methane_flow = 360

co2_24h = True
water_24h = True
react_24h = True
# Otherwise operate for 8 hours, store 24h worth of supply

if co2_24h:
    co2_density = 12
else:
    co2_density = 1101
water_density = 1000

h24 = 24*3600
h12 = 12*3600
h8 = 8*3600
h6 = 6*3600
h4 = 4*3600

components = []

if water_24h:
    rw = RodwellWaterScaled(water_flow/h24, water_24h)
else:
    rw = RodwellWaterScaled(water_flow/h8, water_24h)
    bfw = BufferTank(water_flow/water_density)
    components.append(bfw)
components.append(rw)



# Use CO2 evaporation to cool gases
if co2_24h:
    temperature = 298 
else:
    temperature = 298 - (298-210) * co2_flow/(methane_flow+ox_flow)

if react_24h:
    el1 = WaterElectrolyserScaled(water_flow/h24, react_24h)
    sab = ReactorSeparatorScaled(methane_flow/h24, react_24h)
    # About 12.5% of the water flow is recycled through the electrolyser
    el2 = WaterElectrolyserScaled(water_flow/(8*h24), react_24h)

    cc_meth = CryocoolScaled(methane_flow/h24, react_24h, gas_inlet_temperature=temperature)
    cc_ox = CryocoolScaled(ox_flow/h24, react_24h, gas_inlet_temperature=temperature, liquid_outlet_temperature=90)
else:
    el1 = WaterElectrolyserScaled(water_flow/h8, react_24h)
    sab = ReactorSeparatorScaled(methane_flow/h8, react_24h)
    # About 12.5% of the water flow is recycled through the electrolyser
    el2 = WaterElectrolyserScaled(water_flow/(8*h8), react_24h)

    cc_meth = CryocoolScaled(methane_flow/h8, react_24h, gas_inlet_temperature=temperature)
    cc_ox = CryocoolScaled(ox_flow/h8, react_24h, gas_inlet_temperature=temperature, liquid_outlet_temperature=90)

components.append(el1)
components.append(sab)
components.append(cc_meth)
components.append(cc_ox)

total_W = 0
total_P_day = 0
total_P_night = 0
total_C_day = 0
total_C_night = 0
for i in components:
    total_W += i.weight
    if i.h24:
        total_P_day += i.power
        total_P_night += i.power
        total_C_day += i.cooling
        total_C_night += i.cooling
    else:
        if i.night:
            total_P_night += i.power
            total_C_night += i.cooling
        else:
            total_P_day += i.power
            total_C_day += i.cooling

total_P_capacity = total_P_day + total_P_night
total_E = total_P_night/1000 * 12

total_C_capacity = max(total_C_day, total_C_night)

print("W", int(total_W), "P", int(total_P_capacity), "C", int(total_C_capacity), "E", int(total_E))
print("ESM", total_W+.149*total_P_capacity+.121*total_C_capacity+5*total_E)