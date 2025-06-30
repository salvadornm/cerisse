import math

# Constants (if needed)
R_AIR = 287.0  # Specific gas constant for air [J/(kg·K)]

class nozzle_functions:
    @staticmethod
    def Pstag(P: float, M: float, gamma: float) -> float:
        """Compute stagnation pressure"""
        term = 1.0 + 0.5 * (gamma - 1.0) * M**2
        exponent = gamma / (gamma - 1.0)
        return P * math.pow(term, exponent)

    @staticmethod
    def Tstag(T: float, M: float, gamma: float) -> float:
        """Compute stagnation temperature"""
        return T * (1.0 + 0.5 * (gamma - 1.0) * M**2)

    @staticmethod
    def Pchok(P0: float, gamma: float) -> float:
        """Compute choked pressure (Mach 1)"""
        exponent = gamma / (gamma - 1.0)
        return P0 * math.pow(2.0 / (gamma + 1.0), exponent)

    @staticmethod
    def Tchok(T0: float, gamma: float) -> float:
        """Compute choked temperature (Mach 1)"""
        return T0 * (2.0 / (gamma + 1.0))

    @staticmethod
    def masschok(T0: float, P0: float, gamma: float) -> float:
        """Compute choked mass flow rate [kg/s]"""
        gamma_p1_o2 = 0.5 * (gamma + 1.0)
        term = math.pow(1.0 / gamma_p1_o2, gamma_p1_o2 / (gamma - 1.0))
        return P0 * math.sqrt(gamma / (R_AIR * T0)) * term

    @staticmethod
    def Pnozz(Mach: float, P0: float, gamma: float) -> float:
        """Compute pressure as a function of Mach number"""
        if Mach < 0.0:
            return 0.0
        factor = 1.0 + 0.5 * (gamma - 1.0) * Mach**2
        return P0 * math.pow(factor, -gamma / (gamma - 1.0))

