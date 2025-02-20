import os
import matplotlib.pyplot as plt

from slater_wavefunction import SlaterWaveFunction
from shielding_constant_calculator import ORBITALS_GROWTH_ORDERING
from constants import *

def create_output_directory():
    try:
        os.mkdir("output")
    except FileExistsError:
        pass

def plot_single_density_functions(atomic_num: int):
    wf = SlaterWaveFunction(atomic_num)
    (distances, density_matrix) = wf.compute_single_density_functions(START, END, STEP_SIZE)
    
    # Plot each density function with all the others. For each one, highlight it in a single plot
    for i, densities in enumerate(density_matrix):
        plt.figure(figsize=(10, 10))
        for j, other_densities in enumerate(density_matrix):
            if j == i:
                plt.plot(distances, other_densities, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}")
            else:
                plt.plot(distances, other_densities, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}", alpha=0.3)

        # add the sum
        plt.plot(distances, np.sum(density_matrix, axis=0), label="Sum", linestyle='--', color='black', alpha=0.5)
        plt.title(f"Density function - Orbital {SLATER_GROWTH_GROUPS[i]}")
        plt.xlabel("Distance from nucleus")
        plt.ylabel("Density")
        plt.legend()
        plt.savefig(f"output/z_{atomic_num}/density_orbital_{SLATER_GROWTH_GROUPS[i]}.pdf")
        plt.close()
    
    # Plot all density functions together
    plt.figure(figsize=(10, 10))
    for i, densities in enumerate(density_matrix):
        plt.plot(distances, densities, label=f"Orbital {SLATER_GROWTH_GROUPS[i]}")
    plt.plot(distances, np.sum(density_matrix, axis=0), label="Sum", linestyle='--', color='black', alpha=0.5)
    plt.title("Density functions")
    plt.xlabel("Distance from nucleus")
    plt.ylabel("Density")
    plt.legend()
    plt.savefig(f"output/z_{atomic_num}/density_functions_all.pdf")
    plt.close()

    # Plot the sum of all density functions
    plt.figure(figsize=(10, 10))
    plt.plot(distances, np.sum(density_matrix, axis=0), label="Sum", color='black', alpha=1)
    plt.title("Sum of density functions")
    plt.xlabel("Distance from nucleus")
    plt.ylabel("Density")
    plt.legend()
    plt.savefig(f"output/z_{atomic_num}/density_sum.pdf")
    plt.close()

def plot_single_electron_functions(atomic_num: int):
    wf = SlaterWaveFunction(atomic_num)
    (distances, electron_matrix) = wf.compute_single_electron_functions(START, END, STEP_SIZE)

    # Plot each electron function with all the others. For each one, highlight it in a single plot
    for i, electrons in enumerate(electron_matrix):
        plt.figure(figsize=(10, 10))
        for j, other_electrons in enumerate(electron_matrix):
            if j == i:
                plt.plot(distances, other_electrons, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}")
            else:
                plt.plot(distances, other_electrons, label=f"Orbital {SLATER_GROWTH_GROUPS[j]}", alpha=0.3)
        plt.plot(distances, np.sum(electron_matrix, axis=0), label="Sum", linestyle='--', color='black', alpha=0.5)
        plt.title(f"Electron function - Orbital {SLATER_GROWTH_GROUPS[i]}")
        plt.xlabel("Distance from nucleus")
        plt.ylabel("Electrons")
        plt.legend()
        plt.savefig(f"output/z_{atomic_num}/electron_orbital_{SLATER_GROWTH_GROUPS[i]}.pdf")

    # Plot all electron functions together
    plt.figure(figsize=(10, 10))
    for i, electrons in enumerate(electron_matrix):
        plt.plot(distances, electrons, label=f"Orbital {SLATER_GROWTH_GROUPS[i]}")
    plt.plot(distances, np.sum(electron_matrix, axis=0), label="Sum", linestyle='--', color='black', alpha=0.5)
    plt.title("Electron functions")
    plt.xlabel("Distance from nucleus")
    plt.ylabel("Electrons")
    plt.legend()
    plt.savefig(f"output/z_{atomic_num}/electron_all.pdf")
    plt.close()

    # Plot the sum of all electron functions
    plt.figure(figsize=(10, 10))
    plt.plot(distances, np.sum(electron_matrix, axis=0), label="Sum", color='black', alpha=1)
    plt.title("Sum of electron functions")
    plt.xlabel("Distance from nucleus")
    plt.ylabel("Electrons")
    plt.legend()
    plt.savefig(f"output/z_{atomic_num}/electron_sum.pdf")
    plt.close()

def process_atomic_number(z):
    try:
        os.mkdir(f"output/z_{z}")
    except FileExistsError:
        pass
    print(f"=== ATOMIC NUMBER {z} ===")
    slater_wf = SlaterWaveFunction(z)
    (distances, densities) = slater_wf.compute_density_in_interval(START, END, STEP_SIZE)
    (_, electrons) = slater_wf.compute_electrons_in_interval(START, END, STEP_SIZE)

    plot_density_and_electrons(densities, electrons, z)

    return densities, electrons


def plot_density_and_electrons(densities, electrons, z):
    fig, ax1 = plt.subplots(figsize=(10, 10))
    color = 'tab:red'
    ax1.set_xlabel('Distance from nucleus')
    ax1.set_ylabel('Density', color=color)
    ax1.plot(densities, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Electrons', color=color)
    ax2.plot(electrons, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()
    plt.title(f"Density and electrons - Atomic number {z}")
    plt.savefig(f"output/z_{z}/density_electrons.pdf")
    plt.close()

def plot_all_densities(all_densities):
    plt.figure(figsize=(10, 10))
    plt.vlines(0, 0, max(all_densities[0]) * 1.1, colors='gray', linestyles='dashed', label="Nucleus", alpha=0.4)
    for i, densities in enumerate(all_densities):
        plt.plot(densities, label=f"Atomic number {LOWEST_Z + i}")
    plt.title("Densities")
    plt.xlabel("Distance from nucleus")
    plt.xlim([-10, 200])
    plt.legend()
    plt.ylabel("Density")
    plt.savefig("output/densities.pdf")
    plt.close()

def plot_all_electrons(all_electrons):
    plt.figure(figsize=(10, 10))
    plt.vlines(0, 0, max(all_electrons[0]) * 1.1, colors='gray', linestyles='dashed', label="Nucleus", alpha=0.4)
    for i, electrons in enumerate(all_electrons):
        plt.plot(electrons, label=f"Atomic number {LOWEST_Z + i}")
    plt.title("Electrons")
    plt.xlabel("Distance from nucleus")
    plt.legend()
    plt.ylabel("Electrons")
    plt.savefig("output/electrons.pdf")
    plt.close()

def main(atomic_num=8):
    create_output_directory()



    all_densities = []
    all_electrons = []
    for z in range(LOWEST_Z, HIGHEST_Z + 1):
        densities, electrons = process_atomic_number(z)
        all_densities.append(densities)
        all_electrons.append(electrons)
        plot_single_density_functions(z)
        plot_single_electron_functions(z)

    plot_all_densities(all_densities)
    plot_all_electrons(all_electrons)

if __name__ == "__main__":
    LOWEST_Z = 1
    HIGHEST_Z = 26
    START = 0
    END = 10
    STEP_SIZE = 0.01

    main()
