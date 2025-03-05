from gaussian_cube_io import GaussianCubeIO   
from plotting.isodensity_plotter import IsodensityPlotter
from diatomic_molecule import diatomic_molecule_density
import matplotlib.pyplot as plt
import os
import numpy as np



def main():

    file_names = [
        'data/o2_1.2158Angstr/o2.cube',
        'data/o2_1.777Angstr/o2.cube',
        'data/o2_3.554Angstr/o2.cube',
        'data/o2_4.554Angstr/o2.cube'
    ]


    for file_name in file_names:
        print(f"Processing file: {file_name}")
        output_dir = file_name.split('/')[1]
        
        reader = GaussianCubeIO()
        print(f"Output directory: {output_dir}")
        reader.read(file_name)

        x_low = reader.origin[0]
        x_high = x_low + reader.npoints[0] * reader.spacing[0]

        y_low = reader.origin[1]
        y_high = y_low + reader.npoints[1] * reader.spacing[1]

        z_low = reader.origin[2]
        z_high = z_low + reader.npoints[2] * reader.spacing[2]

        x_range = [x_low, x_high]
        y_range = [y_low, y_high]
        z_range = [z_low, z_high]

        delta_x = reader.spacing[0]
        delta_y = reader.spacing[1]
        delta_z = reader.spacing[2]

        density_qm = reader.data



        density_field = diatomic_molecule_density(
            reader.atoms[0], reader.atoms[1], x_range, y_range, z_range, 'vdw', delta_x, delta_y, delta_z
        )
        
        plotter = IsodensityPlotter(
            x_range, y_range, z_range, delta_x, delta_y, delta_z, reader.atoms, os.path.join('output', output_dir, "qm")
        )

        plotter.plot_sliced_density_field_2D(density_qm, 0, 'xz')
        plotter.plot_sliced_density_isosurface_2D(density_qm, 0.003, 0, 'xz')
        # plotter.plot_isodensity_3D(density_qm, 0.05)
        plotter.set_output_dir(os.path.join('output', output_dir, "model"))
        plotter.plot_sliced_density_isosurface_2D(density_field, 0.003, 0, 'xz')

        plotter.plot_sliced_density_field_2D(density_field, 0, 'xz')

        plotter.set_output_dir(os.path.join('output', output_dir, 'isodensity_evolution'))
        plotter.plot_sliced_isodensity_evolution_2D([density_field, density_qm], ['model', 'qm'], 0, 'xz', background_field_id=0, target_isodensity=0.2, factor=0.001)
        try:
            plotter.plot_sliced_density_difference_2D(density_field, 'model', density_qm, 'qm', 0, 'xz')
        except Exception as e:
            pass            
        # plotter.plot_isodensity_3D(density_field, 0.05)

        # plotter.output_dir = os.path.join('output', output_dir, "model")

        # plotter.create_dir()
        # plotter.plot_sliced_density_field_2D(density_field, 0, 'xy')

if __name__ == "__main__":
    main()