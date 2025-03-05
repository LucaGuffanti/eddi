from atom_descriptor import AtomDescriptor
import os
from config import *
import pandas as pd
from constants import *

atomic_table = pd.read_csv(PERIODIC_TABLE_PATH)
atomic_table = atomic_table[['AtomicNumber', 'Symbol']]
atomic_table['Symbol'] = atomic_table['Symbol'].str.upper()

class PDBReader:
    """Class to read atoms from PDB files."""
    def __init__(self, verbose=False):
        self.atoms = []
        self.verbose = verbose
    
    def read_file(self, filename):
        """Reads a PDB file to extract the atoms that are present"""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File {filename} not found.")
        
        hcount = 0
        acount = 0
        with open(filename, 'r') as f:
            lines = f.readlines()

            for line in lines:
                components = [component.strip() for component in line.split()]
                if components[0] == 'ATOM' or components[0] == 'HETATM':

                    if components[0] == 'ATOM':
                        acount += 1
                    else:
                        hcount += 1

                    number = atomic_table[(atomic_table['Symbol']) == components[-1]]['AtomicNumber'].values[0]   
                    px = float(components[6]) * BOHRS_PER_ANGSTROMS
                    py = float(components[7]) * BOHRS_PER_ANGSTROMS
                    pz = float(components[8]) * BOHRS_PER_ANGSTROMS

                    if self.verbose:
                        print(f"{components[-1]}, Z: {number}, (x, y, z): ({px}, {py}, {pz})")

                    self.atoms.append(
                        AtomDescriptor(
                            number,
                            (px, py, pz)
                        )
                )
                if self.verbose:
                    print(f"Atoms: {acount}, Hetatoms: {hcount}")
            print(f"Read {len(self.atoms)} atoms from {filename} (a={acount}, h={hcount})")
            

if __name__ == "__main__":
    reader = PDBReader(verbose=True)
    reader.read_file('data/ala.pdb')