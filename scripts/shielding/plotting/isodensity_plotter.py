from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import os
from skimage import measure
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from matplotlib.colors import LogNorm

class IsodensityPlotter:
    def __init__(self, x_range, y_range, z_range, dx, dy, dz, atoms, output_dir):
        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        

        self.atoms = atoms
        self.output_dir = output_dir
        self.create_dir()


    def create_dir(self):
        try:
            os.makedirs(self.output_dir)
        except FileExistsError:
            pass

    def from_space_coord_to_grid_coord(self, space_coord, direction):
        if direction not in ['x', 'y', 'z']:
            raise ValueError("Direction must be 'x', 'y', or 'z'")
    
        if direction == 'x':
            return int((space_coord - self.x_range[0]) / self.dx)
        elif direction == 'y':
            return int((space_coord - self.y_range[0]) / self.dy)
        elif direction == 'z':
            return int((space_coord - self.z_range[0]) / self.dz)
        

    def plot_isodensity_3D(self, density_field, isodensity, show_atoms=False):
        print("Rendering 3D Isosurface")

        fig = plt.figure()

        ax = fig.add_subplot(111, projection='3d')

        verts, faces, _, _ = measure.marching_cubes(density_field, level=isodensity)
        ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap='Spectral', lw=1, alpha=0.7, antialiased=True)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.grid(True)
        ax.set_title(f"Isosurface at {isodensity} $e^-/a_0^3$")

        # Add red points for each atom if show_atoms is True
        if show_atoms:
            for atom in self.atoms:
                atom_x = self.from_space_coord_to_grid_coord(atom.position[0], 'x')
                atom_y = self.from_space_coord_to_grid_coord(atom.position[1], 'y')
                atom_z = self.from_space_coord_to_grid_coord(atom.position[2], 'z')
                ax.scatter(atom_x, atom_y, atom_z, color='red', s=50)

        def rotate(angle):
            print("Angle ", angle, end='\r')
            ax.view_init(azim=angle)

        angles = np.arange(0, 360, 10)

        print("Creating animation")
        ani = FuncAnimation(fig, rotate, frames=angles, interval=50)
        ani.save(os.path.join(self.output_dir, f"density_field_i{isodensity}.gif"), writer='imagemagick', fps=5)
        print("Animation completed")
        plt.close()

    def plot_cumulated_density_field_2D(self, density_field, plane='xy', show_atoms=False):
        print("Rendering 2D Cumulated Density Field")
        if plane == 'xy':
            projection = np.sum(density_field, axis=2)
        elif plane == 'yz':
            projection = np.sum(density_field, axis=0)
        elif plane == 'xz':
            projection = np.sum(density_field, axis=1)
        else:
            raise ValueError("Plane must be 'xy', 'yz', or 'xz'")

        print(projection.shape)
        fig, ax = plt.subplots()
        cax = ax.contourf(projection.T, 50, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)

        if show_atoms:
            for atom in self.atoms:
                atom_x = self.from_space_coord_to_grid_coord(atom.position[0], 'x')
                atom_y = self.from_space_coord_to_grid_coord(atom.position[1], 'y')
                atom_z = self.from_space_coord_to_grid_coord(atom.position[2], 'z')

                if plane == 'xy':
                    ax.scatter(atom_x, atom_y, color='red', s=50)
                elif plane == 'yz':
                    ax.scatter(atom_y, atom_z, color='red', s=50)
                elif plane == 'xz':
                    ax.scatter(atom_x, atom_z, color='red', s=50)

        plt.savefig(os.path.join(self.output_dir, f"density_field_cumulated_p{plane}.pdf"))
        ax.set_title(f"Cumulated density field (p={plane}) [$e^-/a_0^3$]")
        plt.close()

    def plot_cumulated_density_isosurface_2D(self, density_field, isodensity, plane='xy', show_atoms=False):
        print("Rendering 2D Cumulated Isosurface")
        if plane == 'xy':
            projection = np.sum(density_field, axis=2)
        elif plane == 'yz':
            projection = np.sum(density_field, axis=0)
        elif plane == 'xz':
            projection = np.sum(density_field, axis=1)
        else:
            raise ValueError("Plane must be 'xy', 'yz', or 'xz'")

        fig, ax = plt.subplots()
        cax = ax.contour(projection.T, levels=[isodensity], colors='red')
        cax = ax.contourf(projection.T, 50, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)


        ax.set_title(f"Cumulated Isosurface at {isodensity} $e^-/a_0^3$ (p={plane})")
        if show_atoms:
            for atom in self.atoms:
                atom_x = self.from_space_coord_to_grid_coord(atom.position[0], 'x')
                atom_y = self.from_space_coord_to_grid_coord(atom.position[1], 'y')
                atom_z = self.from_space_coord_to_grid_coord(atom.position[2], 'z')

                if plane == 'xy':
                    ax.scatter(atom_x, atom_y, color='red', s=50)
                elif plane == 'yz':
                    ax.scatter(atom_y, atom_z, color='red', s=50)
                elif plane == 'xz':
                    ax.scatter(atom_x, atom_z, color='red', s=50)

        plt.savefig(os.path.join(self.output_dir, f"density_field_cumulated_isodensity_p{plane}_i{isodensity}.pdf"))

    def plot_sliced_density_field_2D(self, density_field, slice_coord, plane='xy', show_atoms=False):
        print("Rendering 2D Density Field Slice")

        if plane == 'xy':
            coord = self.from_space_coord_to_grid_coord(slice_coord, 'x')
            projection = density_field[:, :, coord]
        elif plane == 'yz':
            coord = self.from_space_coord_to_grid_coord(slice_coord, 'y')
            projection = density_field[coord, :, :]
        elif plane == 'xz':
            coord = self.from_space_coord_to_grid_coord(slice_coord, 'z')
            projection = density_field[:, coord, :]
        else:
            raise ValueError("Plane must be 'xy', 'yz', or 'xz'")

        fig, ax = plt.subplots()
        ax.set_title(f"Density field slice (p={plane}, s={slice_coord}) [$e^-/a_0^3$]")
        cax = ax.contourf(projection.T, 50, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)


        if show_atoms:
            for atom in self.atoms:
                atom_x = self.from_space_coord_to_grid_coord(atom.position[0], 'x')
                atom_y = self.from_space_coord_to_grid_coord(atom.position[1], 'y')
                atom_z = self.from_space_coord_to_grid_coord(atom.position[2], 'z')

                if plane == 'xy':
                    ax.scatter(atom_x, atom_y, color='red', s=50)
                elif plane == 'yz':
                    ax.scatter(atom_y, atom_z, color='red', s=50)
                elif plane == 'xz':
                    ax.scatter(atom_x, atom_z, color='red', s=50)

        plt.savefig(os.path.join(self.output_dir, f"density_field_slice_p{plane}_s{slice_coord}.pdf"))

    def plot_sliced_density_isosurface_2D(self, density_field, isodensity, slice_coord, plane='xy', show_atoms=False):
        print("Rendering 2D Density Field Slice Isosurface")

        if plane == 'xy':
            coord = self.from_space_coord_to_grid_coord(slice_coord, 'x')
            projection = density_field[:, :, coord]
        elif plane == 'yz':
            coord = self.from_space_coord_to_grid_coord(slice_coord, 'y')
            projection = density_field[coord, :, :]
        elif plane == 'xz':
            coord = self.from_space_coord_to_grid_coord(slice_coord, 'z')
            projection = density_field[:, coord, :]
        else:
            raise ValueError("Plane must be 'xy', 'yz', or 'xz'")

        fig, ax = plt.subplots()
        ax.set_title(f"Density field slice Isosurface at {isodensity} $e^-/a_0^3$ (p={plane}, s={slice_coord})")
        if show_atoms:
            for atom in self.atoms:
                atom_x = self.from_space_coord_to_grid_coord(atom.position[0], 'x')
                atom_y = self.from_space_coord_to_grid_coord(atom.position[1], 'y')
                atom_z = self.from_space_coord_to_grid_coord(atom.position[2], 'z')

                if plane == 'xy':
                    ax.scatter(atom_x, atom_y, color='red', s=50)
                elif plane == 'yz':
                    ax.scatter(atom_y, atom_z, color='red', s=50)
                elif plane == 'xz':
                    ax.scatter(atom_x, atom_z, color='red', s=50)

        cax = ax.contour(projection.T, levels=[isodensity], colors='red')
        cax = ax.contourf(projection.T, 50, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)
        
        plt.savefig(os.path.join(self.output_dir, f"density_field_slice_isodensity_p{plane}_s{slice_coord}_i{isodensity}.pdf"))
