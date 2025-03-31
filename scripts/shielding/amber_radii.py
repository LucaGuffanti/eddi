from constants import *
import pandas as pd
import config
from slater_wavefunction import SlaterWaveFunction
from plotting.plot_radii import plot_amber
import numpy as np


actual_symbols = []

def compare_amber_vdw():

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
    
def check_gradient_at_amber_radius():
    atoms = pd.read_csv(config.VDW_RADII_PATH)
    amber_radii = AMBER_RADII_BOHRS
    atomic_numbers = amber_radii.keys()
    actual_symbols = atoms[atoms['AtomicNumber'].isin(atomic_numbers)]['Symbol'].values

    
    gradients = []
    for am in atomic_numbers:
        sf = SlaterWaveFunction(am)
        radius = amber_radii[am]
        epsilon = 1e-7
        f_r_plus_epsilon = sf.density(radius + epsilon)
        f_r_minus_epsilon = sf.density(radius - epsilon)

        gradient = (f_r_plus_epsilon - f_r_minus_epsilon) / (2 * epsilon)
        gradients.append(gradient)
        print(f"Gradient at {radius} for {am} is {gradient}")

    max_gradient_diff = max(gradients) - min(gradients)
    print(f"Maximum difference in gradient: {max_gradient_diff}")

    import matplotlib.pyplot as plt

    print(actual_symbols)
    plt.figure(figsize=(10, 6))
    plt.plot(gradients, marker='o', linestyle='--', color='b', alpha=0.7, linewidth=1)
    for i, gradient in enumerate(gradients):
        plt.text(i, gradient + 0.0005, f"{gradient:.2e}", fontsize=8, ha='center', va='bottom', color='blue')
    plt.title("Gradient of Density function at AMBER Radius")
    plt.xticks(range(len(gradients)), labels=actual_symbols)
    plt.xlabel("Atomic Number")
    plt.ylabel("Density function gradient (center difference)")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.ylim([min(gradients)-0.01, max(gradients) + 0.01])
    plt.axhline(y=max(gradients), color='r', linestyle='--', label=f"Max Gradient: {max(gradients):.2e}")
    plt.axhline(y=min(gradients), color='g', linestyle='--', label=f"Min Gradient: {min(gradients):.2e}")


    # Annotazione del delta
    plt.text(len(gradients)-0.7, (max(gradients) + min(gradients)) / 2, 
            f"Δ = {max_gradient_diff:.2e}", fontsize=8, color='purple', ha='right', va='bottom', rotation=0)
    
    # Arrow annotation for the gradient difference
    plt.annotate(
        '', 
        xy=(len(gradients) - 0.65, max(gradients)), 
        xytext=(len(gradients) - 0.65, min(gradients)), 
        arrowprops=dict(arrowstyle='<->', color='purple', lw=1.5, alpha=0.7)
    )

    plt.xlim([-0.5, len(gradients)-0.5])
    plt.legend()
    plt.savefig("output/radii/gradient_amber_radius.pdf")

    plt.show()


    

if __name__ == "__main__":
    compare_amber_vdw()
    check_gradient_at_amber_radius()