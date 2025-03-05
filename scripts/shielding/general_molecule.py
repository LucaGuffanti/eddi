from atom_descriptor import AtomDescriptor
from constants import *
import plotting.isodensity_plotter as isodensity_plotter
import plotting.plot_timings as plot_timings
import tqdm
import os

def _vdw_calculation(atoms, x, y, z):
    density_field = np.zeros((len(x), len(y), len(z)))
    for i, x_val in tqdm.tqdm(enumerate(x), total=len(x), desc='Computing density field'):
        for j, y_val in enumerate(y):
            for k, z_val in enumerate(z):
                pos = (x_val, y_val, z_val)
                for atom in atoms:
                    radius = np.linalg.norm(np.array(pos) - atom.position)
                    if radius <= atom.vdw_radius:
                        density_field[i, j, k] += atom.radial_coordinate_density(radius)    

    return density_field

def _no_cutoff_calculation(atoms, x, y, z):
    density_field = np.zeros((len(x), len(y), len(z)))
    for i, x_val in tqdm.tqdm(enumerate(x), total=len(x), desc='Computing density field'):
        for j, y_val in enumerate(y):
            for k, z_val in enumerate(z):
                pos = (x_val, y_val, z_val)
                for atom in atoms:
                    density_field[i, j, k] += atom.cartesian_coordinate_density(pos) + atom.cartesian_coordinate_density(pos)

    return density_field

def polyatomic_molecule_density(atoms, x_range, y_range, z_range, mode='none', delta_x=0.1, delta_y=0.1, delta_z=0.1):
    """Computes the electron density field around a given polyatomic molecule.

    Args:
        atoms (list): List of AtomDescriptor objects representing the atoms in the molecule.
        x_range (tuple): Range of x coordinates (min, max).
        y_range (tuple): Range of y coordinates (min, max).
        z_range (tuple): Range of z coordinates (min, max).
        mode (str, optional): Calculation mode ('none' or 'vdw'). Defaults to 'none'.
        delta_x (float, optional): Step size for x coordinates. Defaults to 0.1.
        delta_y (float, optional): Step size for y coordinates. Defaults to 0.1.
        delta_z (float, optional): Step size for z coordinates. Defaults to 0.1.
    """

    x = np.arange(x_range[0], x_range[1], delta_x)
    y = np.arange(y_range[0], y_range[1], delta_y)
    z = np.arange(z_range[0], z_range[1], delta_z)

    print("Computing electron density field for diatomic molecule")
    if mode == 'vdw':
        print("Using Van-der-Waals cutoff Radii")
        for atom in atoms:
            print("Atom:", atom.atomic_number, "Radius:", atom.vdw_radius)
        return _vdw_calculation(atoms, x, y, z)
    else:
        print("No cutoff radii")
        return _no_cutoff_calculation(atoms, x, y, z)

if __name__ == "__main__":
    atom1 = AtomDescriptor(8, (-0.5, 0.5, 0))
    atom2 = AtomDescriptor(8, (0.5, -0.5, 0))
    atom3 = AtomDescriptor(8, (-0.5, -0.5, 0))
    atom4 = AtomDescriptor(8, (0.5, 0.5, 0))

    atoms = [atom1, atom2, atom3, atom4]

    density_field = polyatomic_molecule_density(atoms, [-2, 2], [-2, 2], [-2, 2],  'none', 0.1, 0.1, 0.1)

    plotter = isodensity_plotter.IsodensityPlotter(
        [-2, 2],
        [-2, 2],
        [-2, 2],
        0.1,
        0.1,
        0.1,
        atoms,
        os.path.join("output", "polyatomic")
    )

    # plotter.plot_sliced_isodensity_evolution_2D([density_field], ['density'], 0, 'xy', target_isodensity=10, factor=0.1) 
    plotter.plot_isodensity_3D(density_field, 5)