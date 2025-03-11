from constants import *
import pandas as pd
import config
from slater_wavefunction import SlaterWaveFunction
from plotting.plot_radii import plot_amber
import numpy as np



atoms = pd.read_csv(config.VDW_RADII_PATH)
amber_radii = AMBER_RADII_BOHRS
vdw_radii = [(radius * BOHRS_PER_ANGSTROMS) / 100 for radius in config.data]

atomic_numbers = amber_radii.keys()
radii = amber_radii.values()


amber_radii_densities = []
vdw_radii_densities = []
actual_symbols = atoms[atoms['AtomicNumber'].isin(atomic_numbers)]['Symbol'].values

print(actual_symbols)

for (atomic_number, radius) in amber_radii.items():
    wf = SlaterWaveFunction(atomic_number)
    density = wf.density(radius)
    amber_radii_densities.append(density)
    
    vdw_radius = vdw_radii[atomic_number - 1]
    density = wf.density(vdw_radius)
    print("Density", density)
    vdw_radii_densities.append(density) 

plot_amber(amber_radii_densities, vdw_radii_densities, actual_symbols)
    
print("Average Isodensity for C-N-O", np.mean(amber_radii_densities[1:3]))
    

