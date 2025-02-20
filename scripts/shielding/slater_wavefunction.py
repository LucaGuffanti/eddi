# slater_wavefunction.py
# This scripts computes the electronic slater wavefunction of an atomic species, starting from the electron configuration of the atom,
# considering the atomic shielding constants following the slater rules.

import numpy as np
import scipy.integrate as spi
import math as math

from constants import *
from shielding_constant_calculator import ShieldingConstantsCalculator

class SlaterWaveFunction():
    """Class that computes the Slater wavefunction for an atomic species, given its atomic number.
    Attributes:
        norm_epsilon (float): Tolerance for normalization constant error.
        number_electrons_epsilon (float): Tolerance for electron number error.
        atomic_number (int): Atomic number of the species.
        shielding_constants_calculator (ShieldingConstantsCalculator): Instance to compute shielding constants.
        electron_configurations (list): List of electron configurations.
        radial_wave_functions (list): List of radial wave functions.
        density_functions (list): List of density functions.
    Methods:
        __init__(atomic_number: int):
            Initializes the SlaterWaveFunction with the given atomic number.
        construct_wavefunction(electron_configuration: str) -> tuple:
            Constructs the Slater radial wavefunction and the total density function for a given electron configuration.
        check_global_density():
            Verifies that the density functions are overall correct by computing the integral throughout the domain.
        density(r: float) -> float:
            Computes the total electron density at a given distance from the nucleus.
        compute_density_in_interval(start: float, end: float, step: float) -> list:
            Computes the electron density in an interval of distances from the nucleus.
        electrons(r: float) -> float:
            Computes the number of electrons at a given distance from the nucleus.
        compute_electrons_in_interval(start: float, end: float, step: float) -> list:
            Computes the number of electrons at an incremental radius from the nucleus by integrating the density functions.
        compute_single_density_functions(start: float, end: float, step: float) -> list:
            Computes the densities for each wavefunction separately.
        compute_single_electrons(start: float, end: float, step: float) -> list:
            Computes the number of electrons for each wavefunction separately.
    """

    def __init__(self, atomic_number: int):
        self.norm_epsilon = 1e-9
        self.number_electrons_epsilon = 1e-9

        self.atomic_number = atomic_number
        self.shielding_constants_calculator = ShieldingConstantsCalculator()
        self.shielding_constants_calculator.compute(atomic_number, atomic_number)

        self.shielding_constants_calculator.print_data()

        self.electron_configurations = list(self.shielding_constants_calculator.SHIELDING_CONSTANTS_TABLE.items())
        print(self.electron_configurations)

        self.radial_wave_functions = list()
        self.density_functions = list()
        for conf in self.electron_configurations:
            # print(conf[1][0])
            (wf, density) = self.construct_wavefunction(conf)
            self.radial_wave_functions.append(wf)
            self.density_functions.append(density)

        self.check_global_density()

    def construct_wavefunction(self, electron_configuration: str) -> tuple:
        """Constructs the Slater radial wavefunction and the total density function for a given electron configuration, considering
        all the electrons that are present in the level.

        Args:
            electron_configuration (str): string describing the electron configuration of the atom

        Returns:
            functions (tuple): a tuple containing the wavefunction and the density function
        """
        
        print("=== Constructing wavefunction for :", electron_configuration[0][1],",",electron_configuration[0][2], "===")
        
        n_star = int(electron_configuration[0][2][0])
        zeta   = electron_configuration[1][3]
        n_e    = electron_configuration[1][0]

        print(" Orbital     n   = ", electron_configuration[0][2][0])
        print(" Orbital     l   = ", electron_configuration[0][2][1])
        print(" Slater orb  n_* = ", n_star)
        print(" Zeta        ζ   = ", zeta)
        print(" Electrons   n_e = ", n_e)

        wf = lambda r: r**(n_star - 1) * np.exp(-zeta * r)
        # Normalize the wavefunction in the interval [0, +inf]

        print("----Wave function data----")
        integral, _ = spi.quad(lambda r: wf(r) ** 2 * r**2, 0, np.inf)
        print(" ∫Ψ*(r)Ψ(r)r^2dr = ", integral)
        norm_constant_computed = 1 / np.sqrt(integral)
        norm_constant_theoretical = (2 * zeta) ** (n_star + 0.5) / np.sqrt(math.factorial(2 * n_star))
        print(" Th. N           = ", norm_constant_theoretical)
        print(" Compt. N        = ", norm_constant_computed)
        print(" Error      |e|  = ", abs(norm_constant_theoretical - norm_constant_computed))
        error_match = abs(norm_constant_theoretical - norm_constant_computed) < self.norm_epsilon
        print(" check: |e| < ε_n ?  = ", abs(norm_constant_theoretical - norm_constant_computed) < self.norm_epsilon)
        
        n_constant = 0
        if (error_match):
            print("    -> Will use computed normalization constant.")
            n_constant = norm_constant_computed
        else:
            print("    -> !!!Large error. Check calculations!!!")
            print("       Will use theoretical normalization constant.")
            n_constant = norm_constant_theoretical

        wf = lambda r: n_constant * r**(n_star - 1) * np.exp(-zeta * r)
        density = lambda r: n_e * (n_constant * r**(n_star - 1) * np.exp(-zeta * r))**2

        print("----Density Wave function check----")
        integral, _ = spi.quad(lambda r: density(r) * r**2, 0, np.inf)
        print(" ∫ρ(r)r^2dr  = ", integral)
        print(" n_e         = ", n_e)
        print(" check: |n_e - ∫ρ(r)r^2dr| < ε_r  = ", abs(n_e - integral) < self.number_electrons_epsilon)

        return (wf, density)

    def check_global_density(self):
        """Verifies that the density functions are overall correct by computing the integral throughout the domain
        """

        print("=== Global Density Check ===")

        computed_electrons = 0
        for density in self.density_functions:
            integral, _ = spi.quad(lambda r: density(r) * r**2, 0, np.inf)
            computed_electrons += integral

        print(" Σ∫ρ(r)r^2dr = ", computed_electrons)
        print(" n_e         = ", self.atomic_number)
        print(" check: |n_e - Σ∫ρ(r)r^2dr| < ε_r = ", abs(self.atomic_number - computed_electrons) < self.number_electrons_epsilon)

    def density(self, r: float) -> float:
        """Computes the total electron density at a given distance from the nucleus.

        Args:
            r (float): distance from the nucleus

        Returns:
            float: total electron density at the given distance
        """
        return sum([density(r) for density in self.density_functions])
    
    def electrons(self, r:float) -> float:
        """Computes the number of electrons at a given distance from the nucleus.
        
        Args:
            r (float): distance from the nucleus

        Returns:
            float: total electron number at the given distance
        """
        electrons = 0
        for density in self.density_functions:
            integral, _ = spi.quad(lambda r: density(r) * r**2, 0, r)
            electrons += integral
        return electrons

    def compute_density_in_interval(self, start: float, end: float, step: float) -> list:
        """Computes the electron density in an interval of distances from the nucleus.

        Args:
            start (float): starting distance from the nucleus
            end (float): ending distance from the nucleus
            step (float): step size

        Returns:
            tuple: list of distances and list of electron densities at the given distances
        """
        assert start <= end, "Start distance must be less than or equal to end distance."
        assert step > 0, "Step size must be positive."

        res_dist = []
        res_dens = []

        for r in np.arange(start, end, step):
            density = self.density(r)  
            res_dist.append(r)
            res_dens.append(density)        
        return (res_dist, res_dens)
    
    def compute_electrons_in_interval(self, start: float, end: float, step: float) -> list:
        """Computes the number of electrons at an incremental radius from the nucleus by integrating the density functions.

        Args:
            start (float): starting distance from the nucleus
            end (float): ending distance from the nucleus
            step (float): step size

        Returns:
            list: list of distances and list of computed electrons at the given distance
        """
        
        assert start <= end, "Start distance must be less than or equal to end distance."
        assert step > 0, "Step size must be positive."

        res_dist = []
        res_elec = []

        for r in np.arange(start, end, step):
            res_elec.append(self.electrons(r))
            res_dist.append(r)
        return (res_dist, res_elec)
    
    def compute_single_density_functions(self, start: float, end: float, step: float) -> list:
        """Computes the density functions for the individual electrons in the atom.

        Args:
            start (float): starting distance from the nucleus
            end (float): ending distance from the nucleus
            step (float): step size

        Returns:
            tuple: list of distances and matrix of electron densities at the given distances
        """
        assert start <= end, "Start distance must be less than or equal to end distance."
        assert step > 0, "Step size must be positive."

        res_dist = [r for r in np.arange(start, end, step)]
        res_matrix = []
        res_dens = []

        for density in self.density_functions:
            res_dens = []
            for r in res_dist:
                res_dens.append(density(r))
            res_matrix.append(res_dens)
        
        return (res_dist, res_matrix)
    
    def compute_single_electron_functions(self, start: float, end: float, step: float) -> list:
        """Computes the number of electrons for each wavefunction separately.

        Args:
            start (float): Starting distance from the nucleus
            end (float): Ending distance from the nucleus
            step (float): Step size

        Returns:
            list: List of distances and matrix of computed electrons at the given distances for each wavefunction
        """
        res_dist = [r for r in np.arange(start, end, step)]
        res_matrix = []
        res_dens = []

        for density in self.density_functions:
            res_dens = []
            for r in res_dist:
                integral, _ = spi.quad(lambda r: density(r) * r**2, 0, r)
                res_dens.append(integral)
            res_matrix.append(res_dens)
        return (res_dist, res_matrix)

if __name__ == '__main__':
    wf = SlaterWaveFunction(26)
    print(wf.density(1.0))
