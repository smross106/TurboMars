""""
Duplication of MEANGEN functionality for compressor design
"""
from numbers import Real
import numpy as np

class Gas(object):
    def __init__(self):
        pass

class RealGas(Gas):
    def __init__(self, R, gamma):
        self._R = R
        self._gamma = gamma

        self._cp = self.R * self.gamma / (self.gamma - 1)
        self._fga = (self.gamma - 1) / self.gamma
    
    def state_lowv(self, h_o, s, V, p_ref, h_ref, s_ref):
        """
        State for low-speed (M<0.3) 
        """
        gamma = self.gamma
        gg = self.gamma / (self.gamma - 1)
        h = h_o - (0.5 * V**2)
        P = p_ref * np.power(h/h_ref, gg) * np.exp((s - s_ref) / self.R)
        T = h/self.cp
        rho = P / (self.R * T)
        a = np.sqrt(gamma * self.R * T)
        Q = 0.0

        return(P, T, rho, Q, V, gamma, a)


    
    @property
    def R(self):
        return(self._R)
    
    @property
    def gamma(self):
        return(self._gamma)

    @property
    def cp(self):
        return(self._cp)
    
    @property
    def fga(self):
        return(self._fga)
    
    @property
    def h(self, T, p):
        return(self.cp * T)

class Machine(object):
    def __init__(self, gas, nstages, T_o_in, P_o_in, design_basis_radius):
        """"
        Initialise the Machine object, which stores multiple Stages as well as whole-machine paramaters
        
        Inputs:
            gas (Gas):                  Gas object for property functions
            nstages (int):              number of stages
            T_o_in (float):             Stagnation temperature at inlet [K]
            P_o_in (float):             Stagnation pressure at inlet [bar]
            design_basis_radius (str):  Radius reference for design: "H", "M" or "T"
        """

        # Read in the stage inputs
        self.turbo_type = "C"

        self.gas = gas

        # Correct P_o_in to Pa
        self.T_o_in = T_o_in
        self.P_o_in = P_o_in * 1e5

        self.nstages = nstages

        # Only accept design basis radius of "H", "M" or "T"
        assert design_basis_radius.upper() in ["H", "M", "T"], "Design basis radius invalid"
        self.design_basis_radius = design_basis_radius.upper

        # Initialisation calculations

        self.nrows = self.nstages * 2
        self.nstations = self.nrows * 9


        # Gas properties
        self.h_o_in = self.gas.h(self.T_o_in, self.P_o_in)
        

class CompressorStage(object):
    """
    Stores a CompressorStage and runs calculations based on it
    """
    def __init__(self, flo_typ, design_hub, mass_flow, rotation_speed, gas):
        #
        #       General property definitions
        #
        assert flo_typ.upper() in ["AXI", "MIX"], "flo_typ must be 'AXI' or 'MIX'"
        assert design_hub.upper() in ["H", "M", "T"], "design_hub must be 'H', 'M' or 'T'"

        self.flo_typ = flo_typ.upper()
        self.design_hub = design_hub.upper()

        #
        #       Top-level design parameters
        #
        assert type(gas) in [Gas, RealGas], "Gas must be a valid object"

        self.mass_flow = mass_flow
        self.rotation_speed = rotation_speed
        self.omega = self.rotation_speed * np.pi / 30
        self.gas = gas

        #
        #region  Get default parameters (MEANGEN suggested)
        #

        # Thickness of leading and trailing edge relative to axial chord
        self.Tk_LE = 0.02
        self.Tk_TE = 0.01

        # Maximum thickness of stator and rotor relative to axial chord
        self.Tk_max_S = 0.10
        self.Tk_max_R = 0.075

        # Fraction of axial chord at which maximum thickness occurs for stator/rotor
        self.x_Tk_max_S = 0.45
        self.x_Tk_max_R = 0.40

        # Fraction of axial chord over which the LE/TE is modified
        self.x_mod_LE = 0.02
        self.x_mod_TE = 0.01

        # Shape parameter for blade thickness
        self.Tk_type = 2

        # Zweifel coefficient
        self.zweifel = 0.5

        # Diffusion factor
        self.D_fac = 0.35

        # Exponent for transforming axial position
        self.expo = 1.35

        # Quasi-orthogonal angles at LE and TE of stator and rotor. Degrees
        self.Q_R_LE = 88
        self.Q_R_TE = 92
        self.Q_S_LE = 92
        self.Q_S_TE = 88

        # Axial chord of rotor, stator. Metres
        self.axchrd_R = 0.05
        self.axchrd_S = 0.04

        # Gap between rows and stages as a fraction of rotor axial chord
        self.rowgap = 0.25
        self.stagegap = 0.5

        # Deviation angles, incidence angles of rotor and stator. Degrees
        self.Devn_R = 5
        self.Devn_S = 5
        self.A_inc_R = -2
        self.A_in_S = -2

        # Isentropic efficiency estimate
        self.eta = 0.9

        # Blockage factors at rotor LE and stator TE
        self.fblock_R_LE = 0.0
        self.fblock_S_TE = 0.0

        # Smoothing number and factor
        self.Nsmooth = 5
        self.smooth_fac = 0.1



        #endregion

        # Initialise internal flow parameters
    
    def set_inlet_conditions(self, p_o_in, t_o_in):
        pass
    
