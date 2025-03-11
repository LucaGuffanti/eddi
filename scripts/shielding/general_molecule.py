from atom_descriptor import AtomDescriptor
from constants import *
import plotting.isodensity_plotter as isodensity_plotter
import plotting.plot_timings as plot_timings
import tqdm
import os
from joblib import Parallel, delayed
from numba import njit

def compute_density(x, y, z, atoms, mode='none'):
    density = 0.0
    for atom in atoms:
        if mode == 'vdw':
            radius = np.linalg.norm(np.array([x, y, z]) - atom.position)
            if radius <= atom.vdw_radius:
                density += atom.radial_coordinate_density(radius)
        elif mode == '2vdw':
            radius = np.linalg.norm(np.array([x, y, z]) - atom.position)
            if radius <= 2*atom.vdw_radius:
                density += atom.radial_coordinate_density(radius)
        else:
            density += atom.cartesian_coordinate_density((x, y, z))
    return density


def polyatomic_molecule_density(atoms, x_range, y_range, z_range, mode='none', delta_x=0.1, delta_y=0.1, delta_z=0.1, jobs=-1):
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
    elif mode == '2vdw':
        print("Using 2*Van-der-Waals cutoff Radii")
        for atom in atoms:
            print("Atom:", atom.atomic_number, "Radius:", atom.vdw_radius)
    else:
        print("No cutoff radii")

    density = np.zeros((len(x), len(y), len(z)))

    def compute_for_x_axis(idx, x):
        slice = np.zeros((len(y), len(z)))
        for j, y_val in enumerate(y):
            for k, z_val in enumerate(z):
                slice[j, k] = compute_density(x, y_val, z_val, atoms, mode)
        return idx, slice
    
    res = Parallel(n_jobs=jobs)(delayed(compute_for_x_axis)(idx, x_val) for idx, x_val in tqdm.tqdm(enumerate(x), total=len(x), desc="Computing density slices"))


    for i, slice in res:
        density[i, :, :] = slice
    return density

if __name__ == "__main__":
    atom1 = AtomDescriptor(8, (-0.5, 0.5, 0))
    atom2 = AtomDescriptor(8, (0.5, -0.5, 0))
    atom3 = AtomDescriptor(8, (-0.5, -0.5, 0))
    atom4 = AtomDescriptor(8, (0.5, 0.5, 0))

    atoms = [atom1, atom2, atom3, atom4]

    density_field = polyatomic_molecule_density(atoms, [-4, 4], [-4, 4], [-4, 4],  'none', 0.1, 0.1, 0.1)
    density_field_2 = polyatomic_molecule_density(atoms, [-4, 4], [-4, 4], [-4, 4],  'vdw', 0.1, 0.1, 0.1)
    density_field_3 = polyatomic_molecule_density(atoms, [-4, 4], [-4, 4], [-4, 4],  '2vdw', 0.1, 0.1, 0.1)


    plotter = isodensity_plotter.IsodensityPlotter(
        [-4, 4],
        [-4, 4],
        [-4, 4],
        0.1,
        0.1,
        0.1,
        atoms,
        os.path.join("output", "polyatomic")
    )

    # plotter.plot_sliced_isodensity_evolution_2D([density_field], ['density'], 0, 'xy', target_isodensity=10, factor=0.1) 
    plotter.plot_sliced_isodensity_evolution_2D([density_field, density_field_2, density_field_3], ['no cutoff', 'vdw', '2vdw'], 0, 'xy', target_isodensity=2, factor=0.01, store_every=10, background_field_id=0)