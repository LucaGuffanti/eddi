import os
import matplotlib.pyplot as plt
import numpy as np

from constants import *


def plot_single_density_functions(type, atomic_num: int, distances, density_matrix):

    # atomic_radius = RADII[atomic_num]["atomic"] / 100  # Convert from picometers to angstroms
    # vdw_radius = RADII[atomic_num]["vdw"] / 100  # Convert from picometers to angstroms

    # # now from angstroms to bohr
    # atomic_radius_bohr = atomic_radius * 1.8897259886 
    # vdw_radius_bohr = vdw_radius * 1.8897259886

    if type == 'slater':
        GROUP  = SLATER_GROWTH_GROUPS
    else: 
        GROUP  = CLEMENTI_GROWTH_GROUPS

    # Plot each density function with all the others. For each one, highlight it in a single plot
    for i, densities in enumerate(density_matrix):
        plt.figure()
        for j, other_densities in enumerate(density_matrix):
            if j == i:
                plt.plot(distances, other_densities, label=f"Orbital {GROUP[j]}")
            else:
                plt.plot(distances, other_densities, label=f"Orbital {GROUP[j]}", alpha=0.3)

        plt.title(f"Density function - Orbital {GROUP[i]}")
        plt.xlabel("Distance from nucleus [$a_0$]")

        # plt.axvspan(0, atomic_radius_bohr, color='red', linestyle='--', label='Atomic Radius', alpha=0.1)
        # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', linestyle='--', label='Van der Waals Radius', alpha=0.1)

        plt.ylabel("Density [$e^-/a_0^3$]")
        plt.legend()
        plt.savefig(f"output/{type}/z_{atomic_num}/density_orbital_{GROUP[i]}.pdf")
        plt.close()
    
    # Plot all density functions together
    plt.figure()
    for i, densities in enumerate(density_matrix):
        plt.plot(distances, densities, label=f"Orbital {GROUP[i]}")
    plt.plot(distances, np.sum(density_matrix, axis=0), label="Sum", linestyle='--', linewidth=2,color='gray' , alpha=0.5)

    # plt.axvspan(0, atomic_radius_bohr, color='red', linestyle='--', label='Atomic Radius', alpha=0.1)
    # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', linestyle='--', label='Van der Waals Radius', alpha=0.1)

    plt.title("Density functions")
    plt.xlabel("Distance from nucleus [$a_0$]")
    plt.ylabel("Density [$e^-/a_0^3$]")
    plt.legend()
    plt.savefig(f"output/{type}/z_{atomic_num}/density_functions_all.pdf")
    plt.close()


def plot_comparison(wf_types, z, distances, densities):
    os.makedirs(f"output/comparison/z_{z}", exist_ok=True)
    print("Plotting",z)
    fig, axs = plt.subplots(1, 1)
    for i, wf_type in enumerate(wf_types):
        axs.plot(distances, densities[i], label=wf_type, linewidth=2, alpha=0.7)
    axs.set_title(f"Density functions - z = {z}")
    axs.set_xlabel("Distance from nucleus [$a_0$]")
    axs.set_ylabel("Density [$e^-/a_0^3$]")
    axs.legend()
    plt.savefig(f"output/comparison/z_{z}/density_comparison.pdf")
    plt.close()

    # Plot the absolute difference between all pairs of functions
    fig, axs = plt.subplots(1, 1)
    for i in range(len(wf_types)):
        for j in range(i + 1, len(wf_types)):
            abs_diff = np.abs(np.array(densities[i]) - np.array(densities[j]))
            axs.plot(distances, abs_diff, label=f"|{wf_types[i]} - {wf_types[j]}|", alpha=0.7)
    axs.set_title(f"Absolute Differences - z = {z}")
    axs.set_xlabel("Distance from nucleus [$a_0$]")
    axs.set_ylabel("Absolute Difference [$e^-/a_0^3$]")
    axs.legend()
    plt.savefig(f"output/comparison/z_{z}/density_absolute_differences.pdf")
    plt.close()

