import os
import matplotlib.pyplot as plt
import numpy as np

from constants import *


def plot_single_density_functions(atomic_num: int, distances, density_matrix):

    # atomic_radius = RADII[atomic_num]["atomic"] / 100  # Convert from picometers to angstroms
    # vdw_radius = RADII[atomic_num]["vdw"] / 100  # Convert from picometers to angstroms

    # # now from angstroms to bohr
    # atomic_radius_bohr = atomic_radius * 1.8897259886 
    # vdw_radius_bohr = vdw_radius * 1.8897259886

    # Plot each density function with all the others. For each one, highlight it in a single plot
    for i, densities in enumerate(density_matrix):
        plt.figure()
        for j, other_densities in enumerate(density_matrix):
            if j == i:
                plt.plot(distances, other_densities, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}")
            else:
                plt.plot(distances, other_densities, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}", alpha=0.3)

        plt.title(f"Density function - Orbital {SLATER_GROWTH_GROUPS[i]}")
        plt.xlabel("Distance from nucleus [$a_0$]")

        # plt.axvspan(0, atomic_radius_bohr, color='red', linestyle='--', label='Atomic Radius', alpha=0.1)
        # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', linestyle='--', label='Van der Waals Radius', alpha=0.1)

        plt.ylabel("Density [$e^-/a_0^3$]")
        plt.legend()
        plt.savefig(f"output/z_{atomic_num}/density_orbital_{SLATER_GROWTH_GROUPS[i]}.pdf")
        plt.close()
    
    # Plot all density functions together
    plt.figure()
    for i, densities in enumerate(density_matrix):
        plt.plot(distances, densities, label=f"Orbital {SLATER_GROWTH_GROUPS[i]}")
    plt.plot(distances, np.sum(density_matrix, axis=0), label="Sum", linestyle='--', linewidth=2,color='gray' , alpha=0.5)

    # plt.axvspan(0, atomic_radius_bohr, color='red', linestyle='--', label='Atomic Radius', alpha=0.1)
    # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', linestyle='--', label='Van der Waals Radius', alpha=0.1)

    plt.title("Density functions")
    plt.xlabel("Distance from nucleus [$a_0$]")
    plt.ylabel("Density [$e^-/a_0^3$]")
    plt.legend()
    plt.savefig(f"output/z_{atomic_num}/density_functions_all.pdf")
    plt.close()

