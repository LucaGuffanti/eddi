
import os

from slater_wavefunction import SlaterWaveFunction
from constants import *
import plotting.plot_densities as plot_densities
import plotting.plot_electrons as plot_electrons
import plotting.plot_radii as plot_radii
import utils.download_data as download_data

def create_output_directory():
    try:
        os.mkdir("output")
    except FileExistsError:
        pass


BOHR_RADIUS_ANGSTROMS = 0.52917721067 # Angstroms
BOHR_RADIUS_AU = 1.0 # Atomic Units
HYDROGEN_DENSITY_AT_BOHR_RADIUS = 0.0
BOHR_BASED_RADII = []
BOHR_BASED_EPSILON = 1e-13
BOHR_BASED_ERRORS = []
BOHR_BASED_DENSITIES = []
BOHR_BASED_ITERS = []
BOHR_BASED_ELECTRONS = []
def process_bohr_based_atomic_radius(z):
    slater_wf = SlaterWaveFunction(z)
    print(f"=== ATOMIC NUMBER {z} ===")
    
    START = 0
    END = 10
    ITERATIONS = 100

    i = 0
    found = False
    density = 0.0
    while i < ITERATIONS and not found:
        possible_radius = (START + END)/2
        density = slater_wf.density(possible_radius)
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
            
            BOHR_BASED_ELECTRONS.append(slater_wf.electrons(possible_radius))

            found = True

        i += 1 
    if not found:
        print(f"No convergence for z = {z}. r = {possible_radius}, density = {density}")
        BOHR_BASED_RADII.append(possible_radius)
        BOHR_BASED_DENSITIES.append(density)
        BOHR_BASED_ERRORS.append(np.abs(density - HYDROGEN_DENSITY_AT_BOHR_RADIUS))
        BOHR_BASED_ITERS.append(i+1)
        BOHR_BASED_ELECTRONS.append(slater_wf.electrons(possible_radius))
    
    return 

DENSITY_CUTOFF_BASED_RADII = []
DENSITY_CUTOFF_BASED_EPSILON = 1e-13
DENSITY_CUTOFF_BASED_ERRORS = []
DENSITY_CUTOFF_BASED_DENSITIES = []
DENSITY_CUTOFF_BASED_ITERS = []
DENSITY_CUTOFF_BASED_ELECTRONS = []
def process_density_cutoff_based_atomic_radius(z, cutoff_density):
    wf = SlaterWaveFunction(z)

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
            

VDW_BASED_RADII_PM = download_data.load_vdw_radii(True)
VDW_BASED_RADII_BOHRS = [(radius / BOHR_RADIUS_ANGSTROMS) / 100 for radius in VDW_BASED_RADII_PM]
VDW_BASED_DENSITIES = []
VDW_BASED_ELECTRONS = []
def process_vdw_based_atomic_radius():
    for i, radius in enumerate(VDW_BASED_RADII_BOHRS):
        slater_wf = SlaterWaveFunction(i+1)
        VDW_BASED_DENSITIES.append(slater_wf.density(radius))
        VDW_BASED_ELECTRONS.append(slater_wf.electrons(radius))
    return

    

def process_atomic_number(z):
    try:
        os.mkdir(f"output/z_{z}")
    except FileExistsError:
        pass
    print(f"=== ATOMIC NUMBER {z} ===")
    slater_wf = SlaterWaveFunction(z)

    START = 0
    END = 10
    STEP_SIZE = 0.01

    (distances, density_matrix) = slater_wf.compute_single_density_functions(START, END, STEP_SIZE)
    (_, electron_matrix) = slater_wf.compute_single_electron_functions(START, END, STEP_SIZE)
    plot_densities.plot_single_density_functions(z, distances, density_matrix)
    plot_electrons.plot_single_electron_functions(z, distances, electron_matrix)


if __name__ == "__main__":
    create_output_directory()

    LOWEST_Z = 1
    HIGHEST_Z = 30
   
    hydrogen_wf = SlaterWaveFunction(1)
    HYDROGEN_DENSITY_AT_BOHR_RADIUS = hydrogen_wf.density(BOHR_RADIUS_AU)
    print("Hydrogen density at bohr radius", HYDROGEN_DENSITY_AT_BOHR_RADIUS)

    for z in range(LOWEST_Z, HIGHEST_Z + 1):
        process_atomic_number(z)

    BOHR_BASED_RADII.append(BOHR_RADIUS_AU)
    BOHR_BASED_DENSITIES.append(HYDROGEN_DENSITY_AT_BOHR_RADIUS)
    BOHR_BASED_ERRORS.append(0.0)
    BOHR_BASED_ITERS.append(0)
    BOHR_BASED_ELECTRONS.append(hydrogen_wf.electrons(BOHR_RADIUS_AU))

    DENSITY_CUTOFF_BASED_TARGET_DENSITY = 1e-3

    VDW_BASED_RADII_BOHRS = VDW_BASED_RADII_BOHRS[:HIGHEST_Z]

    for z in range(2, HIGHEST_Z + 1):
        process_bohr_based_atomic_radius(z)

    for z in range(1, HIGHEST_Z + 1):    
        process_density_cutoff_based_atomic_radius(z, DENSITY_CUTOFF_BASED_TARGET_DENSITY)
    
    plot_radii.plot_type_based_radii(
        BOHR_BASED_RADII,
        BOHR_BASED_DENSITIES, 
        BOHR_BASED_ELECTRONS,
        BOHR_BASED_ERRORS, 
        BOHR_BASED_ITERS,
        BOHR_BASED_EPSILON,
        'bohr'
        )

    plot_radii.plot_type_based_radii(
        DENSITY_CUTOFF_BASED_RADII,
        DENSITY_CUTOFF_BASED_DENSITIES,
        DENSITY_CUTOFF_BASED_ELECTRONS,
        DENSITY_CUTOFF_BASED_ERRORS,
        DENSITY_CUTOFF_BASED_ITERS,
        DENSITY_CUTOFF_BASED_EPSILON,
        'density',
        DENSITY_CUTOFF_BASED_TARGET_DENSITY
    )

    process_vdw_based_atomic_radius()

    plot_radii.plot_all(
        {
            'bohr': [BOHR_BASED_RADII, BOHR_BASED_DENSITIES, BOHR_BASED_ELECTRONS],
            'vdw' : [VDW_BASED_RADII_BOHRS, VDW_BASED_DENSITIES, VDW_BASED_ELECTRONS],
            'density': [DENSITY_CUTOFF_BASED_RADII, DENSITY_CUTOFF_BASED_DENSITIES, DENSITY_CUTOFF_BASED_ELECTRONS, DENSITY_CUTOFF_BASED_TARGET_DENSITY]
        }
    )