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
        self.set_output_dir(output_dir)

    def set_output_dir(self, output_dir):
        """Set a new output directory and create it if necessary."""
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def from_space_coord_to_grid_coord(self, space_coord, direction):
        """Convert spatial coordinates to grid indices."""
        ranges = {'x': self.x_range, 'y': self.y_range, 'z': self.z_range}
        steps = {'x': self.dx, 'y': self.dy, 'z': self.dz}

        if direction not in ranges:
            raise ValueError("Direction must be 'x', 'y', or 'z'")

        grid_index = round((space_coord - ranges[direction][0]) / steps[direction])
        return int(grid_index)

    def plot_isodensity_3D(self, density_field, isodensity, show_atoms=False):
        """Plot a 3D isosurface of the density field."""
        print("Rendering 3D Isosurface")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        verts, faces, _, _ = measure.marching_cubes(density_field, level=isodensity)
        ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap='Spectral', lw=1, alpha=0.7, antialiased=True)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f"Isosurface at {isodensity} $e^-/a_0^3$")
        ax.grid(True)

        if show_atoms:
            self._plot_atoms(ax)

        self._save_animation(fig, ax, f"density_field_i{isodensity}.gif")

    def _plot_atoms(self, ax, plane=None):
        """Plot atoms in red on the given axis."""
        for atom in self.atoms:
            x, y, z = (
                self.from_space_coord_to_grid_coord(atom.position[0], 'x'),
                self.from_space_coord_to_grid_coord(atom.position[1], 'y'),
                self.from_space_coord_to_grid_coord(atom.position[2], 'z'),
            )

            if plane == 'xy':
                ax.scatter(x, y, color='red', s=50)
            elif plane == 'yz':
                ax.scatter(y, z, color='red', s=50)
            elif plane == 'xz':
                ax.scatter(x, z, color='red', s=50)
            else:
                ax.scatter(x, y, z, color='red', s=50)

    def _save_animation(self, fig, ax, filename):
        """Generate and save a rotating 3D animation."""
        def rotate(angle):
            ax.view_init(azim=angle)
            print(angle)
            if angle % 30 == 0:
                filename = f"frame_{angle}.png"
                fig.savefig(os.path.join(self.output_dir, filename))

        angles = np.arange(0, 360, 10)
        ani = FuncAnimation(fig, rotate, frames=angles, interval=50)
        ani.save(os.path.join(self.output_dir, filename), writer='imagemagick', fps=5)
        plt.close()

    def plot_cumulated_density_field_2D(self, density_field, plane='xy', show_atoms=False):
        """Plot a 2D projection of the cumulated density field."""
        print("Rendering 2D Cumulated Density Field")
        projection = self._get_projection(density_field, plane)

        fig, ax = plt.subplots()
        cax = ax.contourf(projection.T, 50, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)

        if show_atoms:
            self._plot_atoms(ax, plane)

        filename = f"density_field_cumulated_p{plane}.pdf"
        plt.savefig(os.path.join(self.output_dir, filename))
        plt.close()

    def plot_cumulated_density_isosurface_2D(self, density_field, isodensity, plane='xy', show_atoms=False):
        """Plot 2D projection with isodensity contour."""
        print("Rendering 2D Cumulated Isosurface")
        projection = self._get_projection(density_field, plane)

        fig, ax = plt.subplots()
        ax.contour(projection.T, levels=[isodensity], colors='red')
        cax = ax.contourf(projection.T, 50, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)

        if show_atoms:
            self._plot_atoms(ax, plane)

        filename = f"density_field_cumulated_isodensity_p{plane}_i{isodensity}.pdf"
        plt.savefig(os.path.join(self.output_dir, filename))
        plt.close()

    def plot_sliced_density_field_2D(self, density_field, slice_coord, plane='xy', show_atoms=False):
        """Plot a density field slice along a given plane."""
        print("Rendering 2D Density Field Slice")
        projection = self._get_slice(density_field, slice_coord, plane)

        fig, ax = plt.subplots()
        ax.set_title(f"Density field slice (p={plane}, s={slice_coord}) [$e^-/a_0^3$]")
        cax = ax.contourf(projection.T, 100, cmap='viridis', alpha=0.5)
        fig.colorbar(cax, ax=ax)

        if show_atoms:
            self._plot_atoms(ax, plane)

        filename = f"density_field_slice_p{plane}_s{slice_coord}.pdf"
        plt.savefig(os.path.join(self.output_dir, filename))
        plt.close()


    def plot_sliced_density_isosurface_2D(self, density_field, isodensity, slice_coord, plane='xy', show_atoms=False):
        """Plot a 2D isosurface slice along a given plane."""
        print("Rendering 2D Isosurface Slice")
        projection = self._get_slice(density_field, slice_coord, plane)

        fig, ax = plt.subplots()
        ax.set_title(f"Isosurface slice (p={plane}, s={slice_coord}, i={isodensity}) [$e^-/a_0^3$]")
        ax.contour(projection.T, levels=[isodensity], colors='yellow')
        cax = ax.contourf(projection.T, 50, cmap='gray', alpha=0.5)
        fig.colorbar(cax, ax=ax)

        if show_atoms:
            self._plot_atoms(ax, plane)

        filename = f"density_field_slice_isosurface_p{plane}_s{slice_coord}_i{isodensity}.pdf"
        plt.savefig(os.path.join(self.output_dir, filename))
        plt.close()

    def _get_projection(self, density_field, plane):
        """Return a 2D projection of the density field along the specified plane."""
        projections = {'xy': 2, 'yz': 0, 'xz': 1}
        if plane not in projections:
            raise ValueError("Plane must be 'xy', 'yz', or 'xz'")
        return np.sum(density_field, axis=projections[plane])

    def _get_slice(self, density_field, slice_coord, plane):
        """Return a 2D slice of the density field at the given coordinate along the specified plane."""
        slices = {'xy': ('z', 2), 'yz': ('x', 0), 'xz': ('y', 1)}
        if plane not in slices:
            raise ValueError("Plane must be 'xy', 'yz', or 'xz'")

        direction, axis = slices[plane]
        coord = self.from_space_coord_to_grid_coord(slice_coord, direction)
        return np.take(density_field, coord, axis=axis)


    def plot_sliced_isodensity_evolution_2D(self, fields, names, slice_coord, plane, 
                                            background_field_id=-1, target_isodensity=1.0, factor=0.01, store_every=-1, show_atoms=False):
        """Generate a GIF showing the evolution of isodensities for each field in fields."""
        print("Rendering Isodensity Evolution GIF")

        fig, ax = plt.subplots()
        ax.set_title(f"Isodensity Evolution (p={plane}, s={slice_coord}) [$e^-/a_0^3$]")

        colors = plt.cm.viridis(np.linspace(0, 1, len(fields)))
        isodensity_levels = np.arange(0, target_isodensity + factor, factor)

        projections = [self._get_slice(density_field, slice_coord, plane) for density_field in fields]
        background_projection = projections[background_field_id] if background_field_id != -1 else None

        if background_projection is not None:
            cax = ax.contourf(background_projection.T, 50, cmap='gray', alpha=0.5)
            fig.colorbar(cax, ax=ax)

        proxies = [plt.Line2D([0], [0], color=colors[i], lw=2) for i in range(len(fields))]
        labels = names if names else [f"Field {i+1}" for i in range(len(fields))]
        ax.legend(proxies, labels, loc='upper right')

        contour_lines = []
        
        def update(frame):
            nonlocal contour_lines
            
            for line in contour_lines:
                line.remove()
            contour_lines = []

            ax.set_title(f"Isodensity Evolution (p={plane}, s={slice_coord}, i={frame:.3f}) [$e^-/a_0^3$]")
            print(f"Current Isodensity: {frame:.3f}", end='\r')

            for i, projection in enumerate(projections):
                contours = measure.find_contours(projection.T, level=frame)
                for contour in contours:
                    line, = ax.plot(contour[:, 1], contour[:, 0], color=colors[i])
                    contour_lines.append(line)

            if show_atoms:
                self._plot_atoms(ax, plane)

            if store_every > 0 and int(frame / factor) % store_every == 0:
                filename = f"isodensity_evolution_p{plane}_s{slice_coord}_i{frame:.3f}.png"
                plt.savefig(os.path.join(self.output_dir, filename))

            return contour_lines

        ani = FuncAnimation(fig, update, frames=isodensity_levels, repeat=False)

        filename = f"isodensity_evolution_p{plane}_s{slice_coord}.gif"
        os.makedirs(self.output_dir, exist_ok=True)  # Ensure output directory exists
        print("Saving GIF")
        ani.save(os.path.join(self.output_dir, filename), writer='imagemagick', fps=5)
        plt.close()

            

    def plot_sliced_density_difference_2D(self, field_1, name_1, field_2, name_2, slice_coord, plane):
        """Plot the absolute difference between two density fields as a 3D surface plot."""
        
        print("Rendering Density Difference")

        difference = np.abs(field_1 - field_2)
        projection = self._get_slice(difference, slice_coord, plane)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f"Density Difference ({name_1} - {name_2}) (p={plane}, s={slice_coord}) [$e^-/a_0^3$]")
        
        X, Y = np.meshgrid(np.arange(projection.shape[0]), np.arange(projection.shape[1]))
        ax.plot_surface(X, Y, projection.T, cmap='viridis', alpha=0.5)
        fig.colorbar(ax.plot_surface(X, Y, projection.T, cmap='viridis', alpha=0.5), ax=ax)

        filename = f"density_difference_{name_1}_{name_2}_p{plane}_s{slice_coord}.pdf"
        plt.savefig(os.path.join(self.output_dir, filename))
        plt.close()
