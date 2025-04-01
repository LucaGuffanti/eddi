# shieldingconstant_calculator.py
# 
# This scripts uses Slater rules to compute the shielding constants of electrons in an atom, for all
# atoms in a range chosen by the user, defaulted to [1, 6].
#
# Reference:
# (Slater, J. C. Phys. Rev. 36 (1): 57–64. Bibcode:1930PhRv...36...57S. doi:10.1103/PhysRev.36.57).

import matplotlib.pyplot as plt
import os
import csv

from shielding_constants import compute_shielding_constant, compute_clementi_constant
from constants import *


ORBITALS_GROWTH_ORDERING = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p']
SINGLE_ORBITALS = { orbital: [f"{orbital}{i}" for i in range(1, MAX_ELECTRONS_PER_LEVEL[orbital[1]] + 1)] for orbital in ORBITALS_GROWTH_ORDERING }


class ShieldingConstantsCalculator:
    def __init__(self, verbose=False):
        self.ORBITALS_GROWTH_ORDERING  = ORBITALS_GROWTH_ORDERING
        self.SINGLE_ORBITALS           = SINGLE_ORBITALS
        self.SHIELDING_CONSTANTS_TABLE = None
        self.ATOMIC_NUMBER_TO_CONFIGURATION = None
        self.verbose = verbose
        if self.verbose:
            print("=====INITIALIZING ELECTRON CONFIGURATIONS=====")

    def generate_configuration(self, atomic_number: int) -> str:
        """Generates the electron configuration of an atom given its atomic number.
        
        Args:
            atomic_number (int): atomic number of the atom
        
        Returns:
            str: electron configuration of the atom
        """
        configuration = ""
        remaining_electrons = atomic_number

        for orbital in self.ORBITALS_GROWTH_ORDERING:
            if remaining_electrons == 0:
                break

            if remaining_electrons < MAX_ELECTRONS_PER_LEVEL[orbital[1]]:
                configuration += f"{orbital}{remaining_electrons}"
                break

            configuration += f"{orbital}{MAX_ELECTRONS_PER_LEVEL[orbital[1]]} "
            remaining_electrons -= MAX_ELECTRONS_PER_LEVEL[orbital[1]]

        return configuration

    def compute(self, lower_z: int = 1, upper_z: int = 6):
        assert lower_z >= 1 and lower_z <= upper_z, "Invalid range of atomic numbers"
        self.LOWER_BOUND_Z = lower_z
        self.UPPER_BOUND_Z = upper_z

        if self.ATOMIC_NUMBER_TO_CONFIGURATION is None:
            self.ATOMIC_NUMBER_TO_CONFIGURATION = {
                i : self.generate_configuration(i).strip() for i in range(self.LOWER_BOUND_Z, self.UPPER_BOUND_Z + 1)
            }

        if self.verbose:
            print("=====COMPUTING SHIELDING CONSTANTS=====")

        self.SHIELDING_CONSTANTS_TABLE = {
            (self.LOWER_BOUND_Z + z, configuration, target[0]+target[1]) : (compute_shielding_constant(configuration, target[0]+target[1], verbose=self.verbose))
                for z, configuration in enumerate(self.ATOMIC_NUMBER_TO_CONFIGURATION.values())
                for target in configuration.split(" ")  
                if target[1] in ['s', 'd', 'f'] 
        }


    def print_data(self):
        """Prints the data in the shielding constants table"""
       
        print("=====SHIELDING CONSTANTS=====")
        for key, value in self.SHIELDING_CONSTANTS_TABLE.items():
            print(f"AtomicNumber: {key[0]}, Configuration: {key[1]}, TargetElectron: {key[2]}")
            print(f"ElectronsInLevel: {value[0]}, ShieldingConstant: {value[1]:7f}, EffectiveNuclearCharge: {value[2]:7f}, ZetaCoefficient: {value[3]:7f}")
            print("_______________________________________")

    def print_data_to_csv(self, csv_file: str = "shielding_constants.csv"):
        """Prints the data in the shielding constants table to a CSV file"""
        if self.verbose:
            print("=====SAVING DATA TO CSV=====")
        try:
            with open(csv_file, mode='w') as file:
                writer = csv.writer(file)
                writer.writerow(["AtomicNumber", "Configuration", "TargetElectron", "ElectronsInLevel" ,"ShieldingConstant", "EffectiveNuclearCharge", "ZetaCoefficient"])
                for key, value in self.SHIELDING_CONSTANTS_TABLE.items():
                    writer.writerow([key[0], key[1], key[2], value[0], f"{value[1]:7f}", f"{value[2]:7f}", f"{value[3]:7f}"])
            if self.verbose:
                print(f"Data saved to {csv_file}")
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
            return

    def read_data_from_csv(self, csv_file: str = "shielding_constants.csv"):
        """Reads the data in the shielding constants table from a CSV file"""
        if self.verbose:
            print("=====READING DATA FROM CSV=====")
        try:
            with open(csv_file, mode='r') as file:
                reader = csv.reader(file)
                for row in reader:
                    if self.verbose:
                        print(row)
                    if row[0] == "AtomicNumber":
                        continue
                    else:
                        self.SHIELDING_CONSTANTS_TABLE[(int(row[0]), row[1], row[2])] = (float(row[3]), float(row[4]), float(row[5]))
            if self.verbose:
                print(f"Data read from {csv_file}")

        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
            return  

    def compute_clementi(self, lower_z: int = 1, upper_z: int = 6):
        """Computes the Clementi shielding constants for a given range of atomic numbers"""
        self.LOWER_BOUND_Z = lower_z
        self.UPPER_BOUND_Z = upper_z

        if self.ATOMIC_NUMBER_TO_CONFIGURATION is None:
            self.ATOMIC_NUMBER_TO_CONFIGURATION = {
                i : self.generate_configuration(i).strip() for i in range(self.LOWER_BOUND_Z, self.UPPER_BOUND_Z + 1)
            }

        self.SHIELDING_CONSTANTS_TABLE = {
            (self.LOWER_BOUND_Z + z, configuration, target[0]+target[1]) : (compute_clementi_constant(configuration, target[0]+target[1], verbose=self.verbose))
            for z, configuration in enumerate(self.ATOMIC_NUMBER_TO_CONFIGURATION.values())
                for target in configuration.split(" ")
        }

if __name__ == "__main__":
    calc = ShieldingConstantsCalculator(verbose=True)
    calc.compute(1, 6)
    calc.print_data_to_csv()
    calc.read_data_from_csv()

    calc.compute_clementi(1, 6)