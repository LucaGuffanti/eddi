from atom_descriptor import AtomDescriptor
from constants import *
import plotting.isodensity_plotter as isodensity_plotter
import plotting.plot_timings as plot_timings
import tqdm

import numpy as np
import time
import csv

def _vdw_calculation(atom_1, atom_2, x, y, z):
    density_field = np.zeros((len(x), len(y), len(z)))
    for i, x_val in tqdm.tqdm(enumerate(x), total=len(x), desc='Computing density field'):
        for j, y_val in enumerate(y):
            for k, z_val in enumerate(z):
                pos = (x_val, y_val, z_val)
                radius_1 = np.linalg.norm(np.array(pos) - atom_1.position)
                radius_2 = np.linalg.norm(np.array(pos) - atom_2.position)
                if radius_1 <= atom_1.vdw_radius:
                    density_field[i, j, k] += atom_1.radial_coordinate_density(radius_1)
                if radius_2 <= atom_2.vdw_radius:
                    density_field[i, j, k] += atom_2.radial_coordinate_density(radius_2)
                    

    return density_field

def _no_cutoff_calculation(atom_1, atom_2, x, y, z):
    density_field = np.zeros((len(x), len(y), len(z)))
    for i, x_val in tqdm.tqdm(enumerate(x), total=len(x), desc='Computing density field'):
        for j, y_val in enumerate(y):
            for k, z_val in enumerate(z):
                pos = (x_val, y_val, z_val)
                density_field[i, j, k] += atom_1.cartesian_coordinate_density(pos) + atom_2.cartesian_coordinate_density(pos)

    return density_field

def diatomic_molecule_density(atom_1, atom_2, x_range, y_range, z_range, mode='none', delta_x=0.1, delta_y=0.1, delta_z=0.1):
    """Computes the electron density field around a given diatomic molecule.

    Args:
        atom_1 (_type_): First atomic species
        atom_2 (_type_): Second atomic species
        grid_refinement (float, optional): Distance of grid points. Defaults to 0.1.
    """

    x = np.arange(x_range[0], x_range[1], delta_x)
    y = np.arange(y_range[0], y_range[1], delta_y)
    z = np.arange(z_range[0], z_range[1], delta_z)

    print("Computing electron density field for diatomic molecule")
    if mode == 'vdw':
        print("Using Van-der-Waals cutoff Radii")
        print("Atom:", atom_1.atomic_number, "Radius:", atom_1.vdw_radius)
        print("Atom:", atom_2.atomic_number, "Radius:", atom_2.vdw_radius)
        return _vdw_calculation(atom_1, atom_2, x, y, z)
    else:
        print("No cutoff radii")
        return _no_cutoff_calculation(atom_1, atom_2, x, y, z)


def from_space_coord_to_grid_coord(space_coord, direction_range, direction_delta):
    return int((space_coord - direction_range[0]) / direction_delta)

def main():
    atom_1 = AtomDescriptor(8, (0, 0, 1.679022))
    atom_2 = AtomDescriptor(8, (0, 0, -1.679022))

    delta_x = 0.1
    delta_y = 0.1
    delta_z = 0.1

    x_range = [-4, 4]
    y_range = [-4, 4]
    z_range = [-4, 4]

    plotter = isodensity_plotter.IsodensityPlotter(
        x_range=x_range,
        y_range=y_range,
        z_range=z_range,
        dx=delta_x,
        dy=delta_y,
        dz=delta_z,
        atoms=[atom_1, atom_2],
        output_dir='output/diatomic_molecule/no_cutoff'
    )

    density_field = diatomic_molecule_density(atom_1, atom_2, x_range, y_range, z_range, 'none', delta_x, delta_y, delta_z)
    isodensity = 2

    plotter.plot_isodensity_3D(density_field, isodensity)
    plotter.plot_cumulated_density_field_2D(density_field, plane='xy')
    plotter.plot_cumulated_density_field_2D(density_field, plane='xz')
    plotter.plot_cumulated_density_field_2D(density_field, plane='yz')

    plotter.plot_cumulated_density_isosurface_2D(density_field, 5, plane='xy')
    plotter.plot_cumulated_density_isosurface_2D(density_field, 5, plane='xz')
    plotter.plot_cumulated_density_isosurface_2D(density_field, 5, plane='yz')
    
    plotter.plot_sliced_density_field_2D(density_field, 0, plane='xy')
    plotter.plot_sliced_density_field_2D(density_field, 0, plane='xz')
    plotter.plot_sliced_density_field_2D(density_field, 0, plane='yz')
   
    plotter.plot_sliced_density_isosurface_2D(density_field, isodensity, 0, plane='xy')
    plotter.plot_sliced_density_isosurface_2D(density_field, isodensity, 0, plane='xz')
    plotter.plot_sliced_density_isosurface_2D(density_field, isodensity, 0, plane='yz')

    plotter.output_dir = 'output/diatomic_molecule/vdw_cutoff'
    plotter.create_dir()

    density_field = diatomic_molecule_density(atom_1, atom_2, x_range, y_range, z_range, 'vdw', delta_x, delta_y, delta_z)
    isodensity = 2

    plotter.plot_isodensity_3D(density_field, isodensity)
    plotter.plot_cumulated_density_field_2D(density_field, plane='xy')
    plotter.plot_cumulated_density_field_2D(density_field, plane='xz')
    plotter.plot_cumulated_density_field_2D(density_field, plane='yz')

    plotter.plot_cumulated_density_isosurface_2D(density_field, 5, plane='xy')
    plotter.plot_cumulated_density_isosurface_2D(density_field, 5, plane='xz')
    plotter.plot_cumulated_density_isosurface_2D(density_field, 5, plane='yz')
    
    plotter.plot_sliced_density_field_2D(density_field, 0, plane='xy')
    plotter.plot_sliced_density_field_2D(density_field, 0, plane='xz')
    plotter.plot_sliced_density_field_2D(density_field, 0, plane='yz')
   
    plotter.plot_sliced_density_isosurface_2D(density_field, isodensity, 0, plane='xy')
    plotter.plot_sliced_density_isosurface_2D(density_field, isodensity, 0, plane='xz')
    plotter.plot_sliced_density_isosurface_2D(density_field, isodensity, 0, plane='yz')

def time_computation():
    atom_1 = AtomDescriptor(8, (-1.14, 0, 0))
    atom_2 = AtomDescriptor(8, (1.14, 0, 0))
    
    no_cutoff = [] 
    cutoff = []
    sizes = []

    delta_x = 0.1
    delta_y = 0.1
    delta_z = 0.1

    LOW = 2
    HIGH = 10

    for i in range(LOW, HIGH + 1):
        x_range = [-i, i]
        y_range = [-i, i]
        z_range = [-i, i]
        size = (2*i)/delta_x * (2*i)/delta_y * (2*i)/delta_z
        sizes.append(size)

        print("Number of grid total points:", size)
        print(" x: ", 2*i/delta_x, " y: ", 2*i/delta_y, " z: ", 2*i/delta_z) 

        start_no_cutoff = time.time()
        density_field = diatomic_molecule_density(atom_1, atom_2, x_range, y_range, z_range, 'none', delta_x, delta_y, delta_z)
        end_no_cutoff = time.time()
        no_cutoff.append(end_no_cutoff - start_no_cutoff)

        start_cutoff = time.time()
        density_field = diatomic_molecule_density(atom_1, atom_2, x_range, y_range, z_range, 'vdw', delta_x, delta_y, delta_z)
        end_cutoff = time.time()
        cutoff.append(end_cutoff - start_cutoff)

    with open('output/density_computation_times.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['GridSize', 'NoCutoffTime', 'CutoffTime'])
        for size, no_cutoff_time, cutoff_time in zip(sizes, no_cutoff, cutoff):
            writer.writerow([size, no_cutoff_time, cutoff_time])
    plot_timings.plot_density_computation(sizes, cutoff, no_cutoff)


if __name__ == "__main__":
    # time_computation()
    main()