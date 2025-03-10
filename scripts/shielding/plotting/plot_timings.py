import matplotlib.pyplot as plt

def plot_density_computation(sizes, timings, labels):
    for label, timing in zip(labels, timings):
        plt.plot(sizes, timing, label=label, marker='o', linewidth=2)
    plt.xlabel('Total grid points')
    plt.ylabel('Time [$s$]')
    plt.title('Comparison of different timings')
    plt.xticks(sizes)
    plt.xscale('log')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.savefig('output/density_computation_times.pdf')
