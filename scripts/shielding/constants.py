import numpy as np

ORBITAL_TO_QUANTUM_NUMBER = {
    "s" : 0,
    "p" : 1,
    "d" : 2,
    "f" : 3
}

QUANTUM_NUMBER_TO_ORBITAL = {
    0 : "s",
    1 : "p",
    2 : "d",
    3 : "f"
}

MAX_ELECTRONS_PER_LEVEL = {
    "s" : 2,
    "p" : 6,
    "d" : 10,
    "f" : 14
}

SLATER_PRINCIPAL_QUANTUM_NUMBER = {
    1 : 1.0,
    2 : 2.0,
    3 : 3.0,
    4 : 3.7,
    5 : 4.0,
    6 : 4.2 
}

SLATER_GROWTH_GROUPS = ['1s', '2s-p', '3s-p', '4s-p', '3d']

# Highest energy level, equal to the highest quantum number we admit
HIGHEST_ENERGY_LEVEL = 5

# Ordering of the orbitals as prescribed by Slater.
SLATER_ORBITALS_ORDERING = np.array([(i, j) for i in range(HIGHEST_ENERGY_LEVEL + 1) for j in range(i)])

# Number of Bohrs per Angstroms
BOHRS_PER_ANGSTROMS = 1.8897259886

AMBER_RADII_ANGSTROMS = {
    1 : 1.4870, # H
    6 : 1.9080, # C
    7 : 1.8240, # N
    8 : 1.6612, # O
    16 : 2.0000 # S
}

AMBER_RADII_BOHRS = {k : v * BOHRS_PER_ANGSTROMS for k, v in AMBER_RADII_ANGSTROMS.items()}

PM_TO_ANGSTROMS = 1e-2