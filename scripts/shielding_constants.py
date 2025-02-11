# shielding_constants.py
# 
# This script computes the shielding constants of electrons in a molecule, following Slater's rules 
# (Slater, J. C. Phys. Rev. 36 (1): 57–64. Bibcode:1930PhRv...36...57S. doi:10.1103/PhysRev.36.57).


import numpy as np

SHIELDING_CONSTANTS_TABLE = {}

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
    "d" : 8,
    "f" : 14
}

# Highest energy level, equal to the highest quantum number we admit
HIGHEST_ENERGY_LEVEL = 5

# Ordering of the orbitals as prescribed by Slater.
ORBITALS_ORDERING = np.array([(i, j) for i in range(HIGHEST_ENERGY_LEVEL + 1) for j in range(i)])

def configuration_to_tuple(configuration: str) -> tuple:
    """Transforms the electron configuration string into a list of tuples that indicate the quantum numbers of the electrons in the
    configuration
    
    Args:
        configuration (str): electronic configuration of the atom

    Returns:
        tuple: a tuple containing a dictionary and the total number of protons. 
        The dictionary is structured in the following way    \\
        {               \\
            (n,l): n_e  \\
        }               \\
        where 
        - n is the principal quantum number
        - l is the secondary quantum number
        - n_e is the number of electrons in the orbital
    """

    print(configuration)
    orbitals_list = configuration.split(" ")

    ionized_flag = False
    configuration_dict = {}
    total_p = 0
    for idx, orbital in enumerate(orbitals_list):
        n   = int(orbital[0])
        l   = ORBITAL_TO_QUANTUM_NUMBER[orbital[1]]
        n_e = int(orbital[2])


        p = 0
        if idx == len(orbitals_list) - 1:
            p = n_e

        else:
            p = MAX_ELECTRONS_PER_LEVEL[orbital[1]]
            if n_e != p:
                ionized_flag = True

        total_p = total_p + p

        print("Orbital:", orbital, "-> n =", n, "l =", l, "n_e =", n_e, "p =", p)
        configuration_dict[(n, l)] = n_e

    print("_______________________________________")
    if (ionized_flag):
        print("-> Atom is ionized.")
    return (configuration_dict, total_p)


def _compute_shielding_constant(configuration: str, target_electron: tuple) -> tuple:
    """Computes the shielding constant given the parameters of the specific electron as well as the effective 
    nuclear charge.

    Args:
        configuration (str): the electron configuration of the atom
        target_electron (tuple): the quantum numbers of the target electron. Expected as the tuple (n,l)
    Returns:
        tuple: a tuple (s, z_eff). s is the shielding constant and z_eff is the effective nuclear charge.
    """

    n,l = target_electron
    
    # Verify that the quantum numbers actually make sense
    assert(n > 0 and l >= 0 and l < n and "ERROR: QUANTUM NUMBERS ARE NOT COHERENT.")

    configuration_dict, total_p = configuration_to_tuple(configuration)

    # Verify that the electron is actually part of the configuration
    assert((n,l) in configuration_dict.keys() and "ERROR: CHOSEN ELECTRON IS NOT PART OF THE ATOM.")
    
    shielding = 0.0

    # The shielding constant is formed, for any group of electrons, from the following contributions.
    # (a) Nothing from any shell outside the one considered.
    # (b) An amount 0.35 from each other electron in the group considered (except the 1s group, where 0.30 is used)
    remaining_electrons = configuration_dict[(n, l)] - 1
    if n == 1:
        shielding = shielding + remaining_electrons * 0.30
    else:
        if l == 0 or l == 1:
            if (n, 1 - l) in configuration_dict.keys():
                remaining_electrons += configuration_dict[(n, 1 - l)]
        shielding = shielding + remaining_electrons * 0.35
        print("Subtracting 0.35 for", remaining_electrons, "electrons in the valence orbital group")

    # (c) If the shell considered is an s, p shell, an amount 0.85 from each electron with a total quantum 
    # less by one, and an amount 1.00 from each electron still further in
    
    # Start by computing the lower orbitals
    lower_orbitals  = [(i, j) for (i, j) in ORBITALS_ORDERING if (i, j) in configuration_dict.keys() and i < n]
    for (i, j) in ORBITALS_ORDERING:
        if (i, j) in configuration_dict.keys() and i == n:
            if j <= 1 and l <= 1:
                pass
            elif j < l:
                lower_orbitals.append((i, j))

    if l == 0 or l == 1:
        lower_n = n - 1
        

        near_orbitals   = [(i, j) for (i, j) in lower_orbitals if i == n-1]
        inner_orbitals  = [(i, j) for (i, j) in lower_orbitals if i < n - 1]

        remaining_electrons = 0

        for orbital in near_orbitals:
            remaining_electrons = remaining_electrons + configuration_dict[orbital]
        
        shielding = shielding + remaining_electrons * 0.85
        print("Target is s or p. Subtracting 0.85 for", remaining_electrons, "in lower orbital")

        remaining_electrons = 0
        for orbital in inner_orbitals:
            remaining_electrons = remaining_electrons + configuration_dict[orbital]
        
        shielding = shielding + remaining_electrons * 1.00
        print("Substracting 1.00 for", remaining_electrons, "in all other orbitals")
    
    else:

        remaining_electrons = 0
        for orbital in lower_orbitals:
            remaining_electrons = remaining_electrons + configuration_dict[orbital]

        shielding = shielding + remaining_electrons * 1.00
        print("Target is not s nor p. Subtracting 1.00 for", remaining_electrons, "in other orbitals")

    return (shielding, total_p - shielding)

def compute_shielding_constant(configuration: str, target_electron: str) -> tuple:
    """Computes the shielding constant given the parameters of the specific electron as well as the effective 
    nuclear charge.


    Args:
        configuration (str): the electron configuration of the atom
        target_electron (tuple): the quantum numbers of the target electron. Expected as the tuple (n,l)
    Returns:
        tuple: a tuple (s, z_eff). s is the shielding constant and z_eff is the effective nuclear charge.
    """
    n = int(target_electron[0])
    l = ORBITAL_TO_QUANTUM_NUMBER[target_electron[1]]
    print(target_electron, "parsed as n =", n, "l =", l)
    return _compute_shielding_constant(configuration, (n, l))



if __name__ == "__main__":
    
    # Testing with the 2 s/p electrons of carbon
    print(compute_shielding_constant("1s2 2s2 2p2", "2s"))
