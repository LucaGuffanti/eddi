
import os

from slater_wavefunction import SlaterWaveFunction
from clementi_wavefunction import ClementiWaveFunction
from constants import *
import plotting.plot_densities as plot_densities
import plotting.plot_electrons as plot_electrons
import plotting.plot_radii as plot_radii
import utils.download_data as download_data
import config 
from multiprocessing import Pool

def create_output_directory():
    try:
        os.mkdir("output")
    except FileExistsError:
        pass

def get_wf(type, z):
    if type=='slater':
        return SlaterWaveFunction(z)
    elif type=='clementi':
        return ClementiWaveFunction(z)

BOHR_RADIUS_ANGSTROMS = 0.52917721067 # Angstroms
BOHR_RADIUS_AU = 1.0 # Atomic Units
HYDROGEN_DENSITY_AT_BOHR_RADIUS = 0.0
BOHR_BASED_RADII = []
BOHR_BASED_EPSILON = 1e-13
BOHR_BASED_ERRORS = []
BOHR_BASED_DENSITIES = []
BOHR_BASED_ITERS = []
BOHR_BASED_ELECTRONS = []
def process_bohr_based_atomic_radius(type, z):
    wf = get_wf(type, z)
    print(f"=== ATOMIC NUMBER {z} ===")
    
    START = 0
    END = 10
    ITERATIONS = 100

    i = 0
    found = False
    density = 0.0
    while i < ITERATIONS and not found:
        possible_radius = (START + END)/2
        density = wf.density(possible_radius)
        print(f"Testing radius {possible_radius}, density {density}")
        if density < HYDROGEN_DENSITY_AT_BOHR_RADIUS:
            END = possible_radius
        else:
            START = possible_radius

        if np.abs(density - HYDROGEN_DENSITY_AT_BOHR_RADIUS)/HYDROGEN_DENSITY_AT_BOHR_RADIUS < BOHR_BASED_EPSILON:
            print(f"Found atomic radius for z = {z}: {possible_radius}. error = {np.abs(density - HYDROGEN_DENSITY_AT_BOHR_RADIUS)/HYDROGEN_DENSITY_AT_BOHR_RADIUS}. it = {i}")
            BOHR_BASED_RADII.append(possible_radius)
            BOHR_BASED_DENSITIES.append(density)
            BOHR_BASED_ERRORS.append(np.abs(density - HYDROGEN_DENSITY_AT_BOHR_RADIUS)/HYDROGEN_DENSITY_AT_BOHR_RADIUS)
            BOHR_BASED_ITERS.append(i+1)
            
            BOHR_BASED_ELECTRONS.append(wf.electrons(possible_radius))

            found = True

        i += 1 
    if not found:
        print(f"No convergence for z = {z}. r = {possible_radius}, density = {density}")
        BOHR_BASED_RADII.append(possible_radius)
        BOHR_BASED_DENSITIES.append(density)
        BOHR_BASED_ERRORS.append(np.abs(density - HYDROGEN_DENSITY_AT_BOHR_RADIUS))
        BOHR_BASED_ITERS.append(i+1)
        BOHR_BASED_ELECTRONS.append(wf.electrons(possible_radius))
    
    return 

DENSITY_CUTOFF_BASED_RADII = []
DENSITY_CUTOFF_BASED_EPSILON = 1e-13
DENSITY_CUTOFF_BASED_ERRORS = []
DENSITY_CUTOFF_BASED_DENSITIES = []
DENSITY_CUTOFF_BASED_ITERS = []
DENSITY_CUTOFF_BASED_ELECTRONS = []
def process_density_cutoff_based_atomic_radius(type, z, cutoff_density):
    wf = get_wf(type, z)

    print(f"=== ATOMIC NUMBER {z} ===")

    START = 0
    END = 10
    ITERATIONS = 100

    i = 0
    found = False
    density = cutoff_density + 1
    possible_radius = 0
    while np.abs(density - cutoff_density)/cutoff_density > DENSITY_CUTOFF_BASED_EPSILON and i < ITERATIONS:
        possible_radius = (START + END)/2
        density = wf.density(possible_radius)
        print(f"Testing radius {possible_radius}, density {density}")

        if density < cutoff_density:
            END = possible_radius
        else:
            START = possible_radius

        if np.abs(density - cutoff_density)/cutoff_density < DENSITY_CUTOFF_BASED_EPSILON:
            print(f"Found atomic radius for z = {z}: {possible_radius}. error = {np.abs(density - cutoff_density)/cutoff_density}. it = {i}")
            DENSITY_CUTOFF_BASED_RADII.append(possible_radius)
            DENSITY_CUTOFF_BASED_DENSITIES.append(density)
            DENSITY_CUTOFF_BASED_ERRORS.append(np.abs(density - cutoff_density)/cutoff_density)
            DENSITY_CUTOFF_BASED_ITERS.append(i+1)
            
            DENSITY_CUTOFF_BASED_ELECTRONS.append(wf.electrons(possible_radius))

            found = True
        
        i += 1
    
    if not found:
        print(f"No convergence for z = {z}. r = {possible_radius}, density = {density}")
        DENSITY_CUTOFF_BASED_RADII.append(possible_radius)
        DENSITY_CUTOFF_BASED_DENSITIES.append(density)
        DENSITY_CUTOFF_BASED_ERRORS.append(np.abs(density - cutoff_density))
        DENSITY_CUTOFF_BASED_ITERS.append(i+1)
        DENSITY_CUTOFF_BASED_ELECTRONS.append(wf.electrons(possible_radius))
            

VDW_BASED_RADII_PM = config.data
VDW_BASED_RADII_BOHRS = [(radius / BOHR_RADIUS_ANGSTROMS) / 100 for radius in VDW_BASED_RADII_PM]
VDW_BASED_DENSITIES = []
VDW_BASED_ELECTRONS = []
def process_vdw_based_atomic_radius(type):
    for i, radius in enumerate(VDW_BASED_RADII_BOHRS):
        wf = get_wf(type, i+1)
        VDW_BASED_DENSITIES.append(wf.density(radius))
        VDW_BASED_ELECTRONS.append(wf.electrons(radius))
    return

    

def process_atomic_number(type, z):
    try:
        os.makedirs(f"output/{type}/z_{z}")
    except FileExistsError:
        pass
    print(f"=== ATOMIC NUMBER {z} ===")
    wf = get_wf(type, z)

    START = 0
    END = 10
    STEP_SIZE = 0.01

    (distances, density_matrix) = wf.compute_single_density_functions(START, END, STEP_SIZE)
    (_, electron_matrix) = wf.compute_single_electron_functions(START, END, STEP_SIZE)
    plot_densities.plot_single_density_functions(type, z, distances, density_matrix)
    plot_electrons.plot_single_electron_functions(type, z, distances, electron_matrix)


def compare_densities(types, lowest_z, highest_z):
    for z in range(lowest_z, highest_z + 1):
        try:
            os.makedirs(f"output/comparison/z_{z}")
        except FileExistsError:
            pass

        wavefunctions = [get_wf(wf_type, z) for wf_type in types]

        distances, densities = None, []
        for wf in wavefunctions:
            distances, density = wf.compute_density_in_interval(0, 10, 0.01)
            densities.append(density)

        plot_densities.plot_comparison(types, z, distances, densities)

def compare_electrons(types, lowest_z, highest_z):
    for z in range(lowest_z, highest_z + 1):
        try:
            os.makedirs(f"output/comparison/z_{z}")
        except FileExistsError:
            pass

        wavefunctions = [get_wf(wf_type, z) for wf_type in types]

        distances, electrons = None, []
        for wf in wavefunctions:
            distances, electron = wf.compute_electrons_in_interval(0, 10, 0.01)
            electrons.append(electron)

        plot_electrons.plot_comparison(types, z, distances, electrons)


if __name__ == "__main__":
    create_output_directory()

    LOWEST_Z = 1
    HIGHEST_Z = 30
   
    # hydrogen_wf = SlaterWaveFunction(1)
    # HYDROGEN_DENSITY_AT_BOHR_RADIUS = hydrogen_wf.density(BOHR_RADIUS_AU)
    # print("Hydrogen density at bohr radius", HYDROGEN_DENSITY_AT_BOHR_RADIUS)

    # for i in range(LOWEST_Z, HIGHEST_Z + 1):
    #     process_atomic_number('slater', i)

    # BOHR_BASED_RADII.append(BOHR_RADIUS_AU)
    # BOHR_BASED_DENSITIES.append(HYDROGEN_DENSITY_AT_BOHR_RADIUS)
    # BOHR_BASED_ERRORS.append(0.0)
    # BOHR_BASED_ITERS.append(0)
    # BOHR_BASED_ELECTRONS.append(hydrogen_wf.electrons(BOHR_RADIUS_AU))

    # DENSITY_CUTOFF_BASED_TARGET_DENSITY = 1e-3

    # VDW_BASED_RADII_BOHRS = VDW_BASED_RADII_BOHRS[:HIGHEST_Z]

    # for z in range(2, HIGHEST_Z + 1):
    #     process_bohr_based_atomic_radius('slater', z)

    # for z in range(1, HIGHEST_Z + 1):    
    #     process_density_cutoff_based_atomic_radius('slater', z, DENSITY_CUTOFF_BASED_TARGET_DENSITY)
    
    # plot_radii.plot_type_based_radii(
    #     'slater',
    #     BOHR_BASED_RADII,
    #     BOHR_BASED_DENSITIES, 
    #     BOHR_BASED_ELECTRONS,
    #     BOHR_BASED_ERRORS, 
    #     BOHR_BASED_ITERS,
    #     BOHR_BASED_EPSILON,
    #     'bohr'
    #     )

    # plot_radii.plot_type_based_radii(
    #     'slater',
    #     DENSITY_CUTOFF_BASED_RADII,
    #     DENSITY_CUTOFF_BASED_DENSITIES,
    #     DENSITY_CUTOFF_BASED_ELECTRONS,
    #     DENSITY_CUTOFF_BASED_ERRORS,
    #     DENSITY_CUTOFF_BASED_ITERS,
    #     DENSITY_CUTOFF_BASED_EPSILON,
    #     'density',
    #     DENSITY_CUTOFF_BASED_TARGET_DENSITY
    # )

    # process_vdw_based_atomic_radius('slater')

    # plot_radii.plot_all(
    #     'slater',
    #     {
    #         'bohr': [BOHR_BASED_RADII, BOHR_BASED_DENSITIES, BOHR_BASED_ELECTRONS],
    #         'vdw' : [VDW_BASED_RADII_BOHRS, VDW_BASED_DENSITIES, VDW_BASED_ELECTRONS],
    #         'density': [DENSITY_CUTOFF_BASED_RADII, DENSITY_CUTOFF_BASED_DENSITIES, DENSITY_CUTOFF_BASED_ELECTRONS, DENSITY_CUTOFF_BASED_TARGET_DENSITY]
    #     }
    # )

    # plot_radii.plot_comparison_separated(
    #     'slater',
    #     {
    #         'bohr': [BOHR_BASED_RADII, BOHR_BASED_DENSITIES, BOHR_BASED_ELECTRONS],
    #         'vdw' : [VDW_BASED_RADII_BOHRS, VDW_BASED_DENSITIES, VDW_BASED_ELECTRONS],
    #         'density': [DENSITY_CUTOFF_BASED_RADII, DENSITY_CUTOFF_BASED_DENSITIES, DENSITY_CUTOFF_BASED_ELECTRONS, DENSITY_CUTOFF_BASED_TARGET_DENSITY]
    #     }
    # )



    # # CLEMENTI
    # print("\n\n\n\n\nCLEMENTI\n\n\n")

    # hydrogen_wf = ClementiWaveFunction(1)
    # HYDROGEN_DENSITY_AT_BOHR_RADIUS = hydrogen_wf.density(BOHR_RADIUS_AU)
    # print("Hydrogen density at bohr radius", HYDROGEN_DENSITY_AT_BOHR_RADIUS)

    # for i in range(LOWEST_Z, HIGHEST_Z + 1):
    #     process_atomic_number('clementi', i)

    # BOHR_BASED_RADII.append(BOHR_RADIUS_AU)
    # BOHR_BASED_DENSITIES.append(HYDROGEN_DENSITY_AT_BOHR_RADIUS)
    # BOHR_BASED_ERRORS.append(0.0)
    # BOHR_BASED_ITERS.append(0)
    # BOHR_BASED_ELECTRONS.append(hydrogen_wf.electrons(BOHR_RADIUS_AU))

    # DENSITY_CUTOFF_BASED_TARGET_DENSITY = 1e-3

    # VDW_BASED_RADII_BOHRS = VDW_BASED_RADII_BOHRS[:HIGHEST_Z]

    # for z in range(2, HIGHEST_Z + 1):
    #     process_bohr_based_atomic_radius('clementi', z)

    # for z in range(1, HIGHEST_Z + 1):    
    #     process_density_cutoff_based_atomic_radius('clementi', z, DENSITY_CUTOFF_BASED_TARGET_DENSITY)
    
    # plot_radii.plot_type_based_radii(
    #     'clementi',
    #     BOHR_BASED_RADII,
    #     BOHR_BASED_DENSITIES, 
    #     BOHR_BASED_ELECTRONS,
    #     BOHR_BASED_ERRORS, 
    #     BOHR_BASED_ITERS,
    #     BOHR_BASED_EPSILON,
    #     'bohr'
    #     )

    # plot_radii.plot_type_based_radii(
    #     'clementi',
    #     DENSITY_CUTOFF_BASED_RADII,
    #     DENSITY_CUTOFF_BASED_DENSITIES,
    #     DENSITY_CUTOFF_BASED_ELECTRONS,
    #     DENSITY_CUTOFF_BASED_ERRORS,
    #     DENSITY_CUTOFF_BASED_ITERS,
    #     DENSITY_CUTOFF_BASED_EPSILON,
    #     'density',
    #     DENSITY_CUTOFF_BASED_TARGET_DENSITY
    # )

    # process_vdw_based_atomic_radius('clementi')

    # plot_radii.plot_all(
    #     'clementi',
    #     {
    #         'bohr': [BOHR_BASED_RADII, BOHR_BASED_DENSITIES, BOHR_BASED_ELECTRONS],
    #         'vdw' : [VDW_BASED_RADII_BOHRS, VDW_BASED_DENSITIES, VDW_BASED_ELECTRONS],
    #         'density': [DENSITY_CUTOFF_BASED_RADII, DENSITY_CUTOFF_BASED_DENSITIES, DENSITY_CUTOFF_BASED_ELECTRONS, DENSITY_CUTOFF_BASED_TARGET_DENSITY]
    #     }
    # )

    # plot_radii.plot_comparison_separated(
    #     'clementi',
    #     {
    #         'bohr': [BOHR_BASED_RADII, BOHR_BASED_DENSITIES, BOHR_BASED_ELECTRONS],
    #         'vdw' : [VDW_BASED_RADII_BOHRS, VDW_BASED_DENSITIES, VDW_BASED_ELECTRONS],
    #         'density': [DENSITY_CUTOFF_BASED_RADII, DENSITY_CUTOFF_BASED_DENSITIES, DENSITY_CUTOFF_BASED_ELECTRONS, DENSITY_CUTOFF_BASED_TARGET_DENSITY]
    #     }
    # )

    compare_densities(['slater', 'clementi'], LOWEST_Z, HIGHEST_Z)

    compare_electrons(['slater', 'clementi'], LOWEST_Z, HIGHEST_Z)