import os
import matplotlib.pyplot as plt
from constants import *
import csv as csv



def plot_single_electron_functions(atomic_num: int, distances, electron_matrix):

    # atomic_radius = RADII[atomic_num]["atomic"] / 100  # Convert from picometers to angstroms
    # vdw_radius = RADII[atomic_num]["vdw"] / 100  # Convert from picometers to angstroms

    # # now from angstroms to bohr
    # atomic_radius_bohr = atomic_radius * 1.8897259886 
    # vdw_radius_bohr = vdw_radius * 1.8897259886

    # Plot each electron function with all the others. For each one, highlight it in a single plot
    for i, electrons in enumerate(electron_matrix):
        plt.figure()
        for j, other_electrons in enumerate(electron_matrix):
            if j == i:
                plt.plot(distances, other_electrons, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}")
            else:
                plt.plot(distances, other_electrons, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}", alpha=0.3)
        plt.title(f"Electron function - Orbital {SLATER_GROWTH_GROUPS[i]}")

        # plt.axvspan(0, atomic_radius_bohr, color='red', linestyle='--', label='Atomic Radius', alpha=0.1)
        # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', linestyle='--', label='Van der Waals Radius', alpha=0.1)

        plt.xlabel("Distance from nucleus [$a_0$]")
        plt.ylabel("Electrons")
        plt.legend()
        plt.savefig(f"output/z_{atomic_num}/electron_orbital_{SLATER_GROWTH_GROUPS[i]}.pdf")
        plt.close()

    # Plot all electron functions together
    plt.figure()
    for i, electrons in enumerate(electron_matrix):
        plt.plot(distances, electrons, label=f"Orbital {SLATER_GROWTH_GROUPS[i]}", alpha=0.5)
    plt.plot(distances, np.sum(electron_matrix, axis=0), label="Sum", linestyle='--', color='grey', alpha=0.7, linewidth=2)

    # plt.axvspan(0, atomic_radius_bohr, color='red', alpha=0.1, label='Atomic Radius')
    # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', alpha=0.1, label='Van der Waals Radius')

    plt.title("Electron functions - z = " + str(atomic_num))
    plt.xlabel("Distance from nucleus [$a_0$]")
    plt.ylabel("Electrons")
    plt.legend()
    plt.savefig(f"output/z_{atomic_num}/electron_all.pdf")
    plt.close()
