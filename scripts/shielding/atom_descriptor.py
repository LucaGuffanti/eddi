import numpy as np
from slater_wavefunction import SlaterWaveFunction
from constants import *
from utils.download_data import load_vdw_radii
from config import *

data = load_vdw_radii(vdw_radii_path=VDW_RADII_PATH, periodic_table_path=PERIODIC_TABLE_PATH, verbose=True)

class AtomDescriptor:
    """
    A class to describe an atom using its atomic number and position.
    It uses the Slater wave function to calculate the electron density.
    """

    def __init__(self, atomic_number: int, position: tuple, float= 0,charge: float=0):
        """
        Initialize the AtomDescriptor with an atomic number and position.

        Parameters:
        atomic_number (int): The atomic number of the atom.
        position (tuple): The (x, y, z) coordinates of the atom.
        """
        self.atomic_number = atomic_number
        self.position = np.array(position)
        self.slater_wf = SlaterWaveFunction(self.atomic_number)
        self.nuclear_radius = 0
        self.charge = 0

        
        self.vdw_radius = (data[self.atomic_number] / 100) * 1.8897259886
        

    def radial_coordinate_density(self, radius):
        """
        Calculate the electron density at a given radial distance from the nucleus.

        Parameters:
        radius (float): The radial distance from the nucleus.

        Returns:
        float: The electron density at the given radius.
        """
        if radius < self.nuclear_radius:
            return 0
        return self.slater_wf.density(radius)

    def cartesian_coordinate_density(self, position):
        """
        Calculate the electron density at a given Cartesian coordinate.

        Parameters:
        position (tuple): The (x, y, z) coordinates where the density is calculated.

        Returns:
        float: The electron density at the given position.
        """
        position = np.array(position)
        radius = np.linalg.norm(position - self.position)
        if radius < self.nuclear_radius:
            return 0
        return self.slater_wf.density(radius)