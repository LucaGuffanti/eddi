import os
import matplotlib.pyplot as plt
from constants import *
import csv as csv



def plot_single_electron_functions(type, atomic_num: int, distances, electron_matrix):

    # atomic_radius = RADII[atomic_num]["atomic"] / 100  # Convert from picometers to angstroms
    # vdw_radius = RADII[atomic_num]["vdw"] / 100  # Convert from picometers to angstroms

    # # now from angstroms to bohr
    # atomic_radius_bohr = atomic_radius * 1.8897259886 
    # vdw_radius_bohr = vdw_radius * 1.8897259886

    # Plot each electron function with all the others. For each one, highlight it in a single plot
    GROUPS = SLATER_GROWTH_GROUPS if type == 'slater' else CLEMENTI_GROWTH_GROUPS
    for i, electrons in enumerate(electron_matrix):
        plt.figure()
        for j, other_electrons in enumerate(electron_matrix):
            if j == i:
                plt.plot(distances, other_electrons, label=f"Orbital {GROUPS[j]}")
            else:
                plt.plot(distances, other_electrons, label=f"Orbital {GROUPS[j]}", alpha=0.3)
        plt.title(f"Electron function - Orbital {GROUPS[i]}")

        # plt.axvspan(0, atomic_radius_bohr, color='red', linestyle='--', label='Atomic Radius', alpha=0.1)
        # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', linestyle='--', label='Van der Waals Radius', alpha=0.1)

        plt.xlabel("Distance from nucleus [$a_0$]")
        plt.ylabel("Electrons")
        plt.legend()
        plt.savefig(f"output/{type}/z_{atomic_num}/electron_orbital_{GROUPS[i]}.pdf")
        plt.close()

    # Plot all electron functions together
    plt.figure()
    for i, electrons in enumerate(electron_matrix):
        plt.plot(distances, electrons, label=f"Orbital {GROUPS[i]}", alpha=0.5)
    plt.plot(distances, np.sum(electron_matrix, axis=0), label="Sum", linestyle='--', color='grey', alpha=0.7, linewidth=2)

    # plt.axvspan(0, atomic_radius_bohr, color='red', alpha=0.1, label='Atomic Radius')
    # plt.axvspan(atomic_radius_bohr, vdw_radius_bohr, color='yellow', alpha=0.1, label='Van der Waals Radius')

    plt.title("Electron functions - z = " + str(atomic_num))
    plt.xlabel("Distance from nucleus [$a_0$]")
    plt.ylabel("Electrons")
    plt.legend()
    plt.savefig(f"output/{type}/z_{atomic_num}/electron_all.pdf")
    plt.close()


def plot_comparison(wf_types, z, distances, electrons):
    os.makedirs(f"output/comparison/z_{z}", exist_ok=True)
    print("Plotting", z)
    fig, axs = plt.subplots(1, 1)
    for i, wf_type in enumerate(wf_types):
        axs.plot(distances, electrons[i], label=wf_type, linewidth=2, alpha=0.7)
    axs.set_title(f"Electrons - z = {z}")
    axs.set_xlabel("Distance from nucleus [$a_0$]")
    axs.set_ylabel("Electrons")
    axs.legend()
    axs.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.savefig(f"output/comparison/z_{z}/electron_comparison.pdf")
    plt.close()

    # Plot absolute differences between all pairs of functions
    fig, axs = plt.subplots(1, 1)
    for i in range(len(wf_types)):
        for j in range(i + 1, len(wf_types)):
            abs_diff = np.abs(np.array(electrons[i]) - np.array(electrons[j]))
            axs.plot(distances, abs_diff, label=f"|{wf_types[i]} - {wf_types[j]}|", alpha=0.7)
    axs.set_title(f"Absolute Differences - z = {z}")
    axs.set_xlabel("Distance from nucleus [$a_0$]")
    axs.set_ylabel("Absolute Difference")
    axs.legend()
    axs.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.savefig(f"output/comparison/z_{z}/electron_absolute_differences.pdf")
    plt.close()