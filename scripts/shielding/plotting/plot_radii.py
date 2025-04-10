import matplotlib.pyplot as plt
import os
from constants import *

def plot_type_based_radii(wf_type, radii, densities, electrons, errors, iters, epsilon, type, *kwargs):
    
    assert type in ['bohr', 'vdw', 'density'], "Type must be either 'bohr', 'vdw' or 'density'"

    try:
        os.mkdir(f"output/{wf_type}/radii")
    except FileExistsError:
        pass

    # 4 plots one over the other
    fig, axs = plt.subplots(5, 1, figsize=(10, 25))
    
    Z = np.arange(1, len(radii) + 1)
    if type == 'bohr':
        fig.suptitle("Bohr-based atomic radii (w.r.t. hydrogen density at Bohr radius)")
        type_of_error = "Relative error on computed hydrogen density"
        Z_2 = Z[1:]
        errors = errors[1:]
        iters = iters[1:]
    elif type == 'vdw':
        fig.suptitle("Van der Waals radii")
    else:
        Z_2 = Z
        fig.suptitle(f"Density-based atomic radii (w.r.t. density {kwargs[0]})")
        type_of_error = "Relative error on target density"

    axs = axs.flatten()


    periods = [1, 3, 11, 19, 37, 55, 87]  # Starting atomic numbers of new periods
    colors = ['lightgrey', 'lightblue', 'lightgreen', 'lightyellow', 'lightpink', 'lightcoral', 'lightcyan']

    type_of_error = ""


    for i, ax in enumerate(axs):
        ax.plot([Z,Z,Z,Z_2,Z_2][i], [radii, densities, electrons, errors, iters][i],'o-')
        ax.set_xlabel("Atomic number")
        ax.set_ylabel(["Atomic radius [$a_0$]", "Density [$e^-/a_0^3$]", "Electrons", "Relative error", "Iterations"][i])

        ax.set_title(['Atomic radius','Electron density at Atomic Radius','Enclosed electrons at Atomic Radius', type_of_error,'Iterations until convergence'][i])

        for j, period in enumerate(periods):
            if period <= len(radii):
                if j < len(periods) - 1:
                    ax.axvspan(period-0.5, min(periods[j + 1] - 1+0.5, len(radii)+0.5), facecolor=colors[j], alpha=0.2)
                else:
                    ax.axvspan(period-0.5, len(radii)+0.5, facecolor=colors[j], alpha=0.2)
        if (i == 3):
            ax.set_ylim(-2*epsilon, 2*epsilon)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.axhline(y=min(errors), color='gray', linestyle='--', alpha=0.3)
            ax.axhline(y=max(errors), color='gray', linestyle='--', alpha=0.3)
            ax.annotate(
                '', xy=(len(errors)+0.5, min(errors)), xytext=(len(errors)+0.5, max(errors)),
                arrowprops=dict(arrowstyle='<->', color='black')
            )
            ax.text(len(errors) + 0.7, (min(errors) + max(errors)) / 2,
                    f'{max(errors) - min(errors):.2e}', va='center', rotation=90)
        if (i == 4):
            ax.set_ylim(min(iters)-1, max(iters) + 1)
        
        if (i == 2):
            ax.plot(Z, Z, label="Expected electrons", linestyle='--', color='gray', alpha=0.3)
            ax.legend()

    plt.savefig(f"output/{wf_type}/radii/{type}_radii.pdf")
    plt.close()

def plot_all(wf_type, data_dict):
    bohr_data = data_dict['bohr']
    vdw_data = data_dict['vdw']
    density_data = data_dict['density']

    fig, axs = plt.subplots(3, 1, figsize=(10,15))
    axs = axs.flatten()

    fig.suptitle("Comparison")

    bohr_radii = bohr_data[0]


    Z = np.arange(1, len(bohr_radii) + 1)
    periods = [1, 3, 11, 19, 37, 55, 87]  # Starting atomic numbers of new periods
    colors = ['lightgrey', 'lightblue', 'lightgreen', 'lightyellow', 'lightpink', 'lightcoral', 'lightcyan']

    for i, title in enumerate(['Radius [$a_0$]', f'Densities [$e^-/a_0^3$]', 'Electrons']):
        axs[i].plot(Z,bohr_data[i], label="Bohr", linewidth=2, alpha=0.7, marker='o')
        axs[i].plot(Z,density_data[i], label=f"Density {density_data[3]}", linewidth=2, alpha=0.7,marker='s')
        axs[i].plot(Z,vdw_data[i], label="VDW", linewidth=2, alpha=0.7,marker='*')
        axs[i].set_title(['Atomic radius', 'Electron density at Atomic Radius','Enclosed electrons at Atomic Radius'][i])
        axs[i].set_xlabel("Atomic number")
        axs[i].set_ylabel(title)

        for j, period in enumerate(periods):
            if period <= len(bohr_data[i]):
                if j < len(periods) - 1:
                    axs[i].axvspan(period-0.5, min(periods[j + 1] - 1+0.5, len(bohr_data[i])+0.5), facecolor=colors[j], alpha=0.2)
                else:
                    axs[i].axvspan(period-0.5, len(bohr_data[i])+0.5, facecolor=colors[j], alpha=0.2)
        axs[i].legend()    

    plt.savefig(f"output/{wf_type}/radii/comparison.pdf")
    plt.close()

def plot_comparison_separated(wf_type, data_dict):
    bohr_data = data_dict['bohr']
    vdw_data = data_dict['vdw']
    density_data = data_dict['density']

    Z = np.arange(1, len(bohr_data[0]) + 1)
    periods = [1, 3, 11, 19, 37, 55, 87]  # Starting atomic numbers of new periods
    colors = ['lightgrey', 'lightblue', 'lightgreen', 'lightyellow', 'lightpink', 'lightcoral', 'lightcyan']

    titles = ['Atomic radius [$a_0$]', 'Electron density at Atomic Radius [$e^-/a_0^3$]', 'Enclosed electrons at Atomic Radius']
    y_labels = ['Radius [$a_0$]', 'Density [$e^-/a_0^3$]', 'Electrons']

    for i, title in enumerate(titles):
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(Z, bohr_data[i], label="Bohr", linewidth=2, alpha=0.7, marker='o')
        ax.plot(Z, density_data[i], label=f"Density {density_data[3]}", linewidth=2, alpha=0.7, marker='s')
        ax.plot(Z, vdw_data[i], label="VDW", linewidth=2, alpha=0.7, marker='*')
        ax.set_title(title)
        ax.set_xlabel("Atomic number")
        ax.set_ylabel(y_labels[i])

        for j, period in enumerate(periods):
            if period <= len(bohr_data[i]):
                if j < len(periods) - 1:
                    ax.axvspan(period, min(periods[j + 1]-1, len(bohr_data[i])-1), facecolor=colors[j], alpha=0.2)
                else:
                    ax.axvspan(period, len(bohr_data[i]), facecolor=colors[j], alpha=0.2)
        ax.legend()
        plt.savefig(f"output/{wf_type}/radii/comparison_{i}.pdf")
        plt.close()

def plot_amber(wf_type, amber_radii_densities, vdw_radii_densities, symbols):
    fig, ax = plt.subplots()

    atomic_numbers = [i for i in range(len(amber_radii_densities))]
    ax.plot(atomic_numbers, amber_radii_densities, linestyle='--', color='blue', marker='o', label='Densities at AMBER radii', alpha=0.7, linewidth=1)
    ax.plot(atomic_numbers, vdw_radii_densities, linestyle='--', color='red', marker='s', label='Densities at VDW radii', alpha=0.7, linewidth=1)
    plt.xlabel("Atomic number")
    plt.ylabel("Density [$e^-/a_0^3$]")
    plt.title("Comparison of AMBER radii and Slater Wave Function densities")
    ax.set_xticks(atomic_numbers)  # Ensure this is a list of specific atomic numbers you want to display
    ax.set_xticklabels(symbols)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5) 
    for i, txt in enumerate(amber_radii_densities):
        ax.annotate(f'{txt:.4e}', (atomic_numbers[i], amber_radii_densities[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color='blue',
                    bbox=dict(boxstyle="round,pad=0.3", edgecolor='none', facecolor='white', alpha=0.2))

    for i, txt in enumerate(vdw_radii_densities):
        ax.annotate(f'{txt:.4e}', (atomic_numbers[i], vdw_radii_densities[i]), textcoords="offset points", xytext=(0,-15), ha='center', fontsize=8, color='red',
                    bbox=dict(boxstyle="round,pad=0.3", edgecolor='none', facecolor='white', alpha=0.2))
    ax.legend()
    ax.set_xlim([-0.5, len(atomic_numbers)-0.5])
    ax.set_ylim([0, max(max(amber_radii_densities), max(vdw_radii_densities))*1.1])
    
    plt.savefig(f"output/{wf_type}/radii/amber_comparison.pdf")
    plt.close()


def plot_amber_comparison(types, densities, symbols):
    fig, ax = plt.subplots()

    for i, density in enumerate(densities):
        ax.plot(range(len(density)), density, label=types[i], marker='o', linestyle='-', alpha=0.7)

    ax.set_xticks(range(len(symbols)))
    ax.set_xticklabels(symbols)
    ax.set_xlabel("Element")
    ax.set_ylabel("Density [$e^-/a_0^3$]")
    ax.set_title("Density Comparison Across Different WF Types at AMBER Radii")
    ax.legend()
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)

    plt.savefig("output/radii/density_comparison_slater_clementi.pdf")
    plt.show()
    plt.close()
