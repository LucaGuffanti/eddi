import matplotlib.pyplot as plt

def plot_density_computation(sizes, with_cutoff, without_cutoff):
    plt.plot(sizes, with_cutoff, label='With Cutoff', marker='o', linewidth=2)
    plt.plot(sizes, without_cutoff, label='Without Cutoff', marker='s', linewidth=2)
    plt.xlabel('Total grid points')
    plt.ylabel('Time [$s$]')
    plt.title('Comparison between no cutoff and VdW-based cutoff')
    plt.xticks(sizes)
    plt.xscale('log')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.savefig('output/density_computation_times.pdf')
