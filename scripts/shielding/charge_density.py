from general_molecule import polyatomic_molecule_density
from gaussian_cube_io import GaussianCubeIO
from atom_descriptor import AtomDescriptor
from pdb_reader import PDBReader

import os
import numpy as np
from tqdm import tqdm
from plotting import isodensity_plotter

def compute_charge_density(density_field, atoms, x_range, y_range, z_range, dx, dy, dz):
    """
    Computes the charge density field around a polyatomic molecule.
    The nucleus is modeled as a Gaussian charge distribution, and the electron density is subtracted.
    """
    # Create grid coordinates
    x_vals = np.linspace(x_range[0], x_range[1], density_field.shape[0])
    y_vals = np.linspace(y_range[0], y_range[1], density_field.shape[1])
    z_vals = np.linspace(z_range[0], z_range[1], density_field.shape[2])

    X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

    nucleus_density = np.zeros_like(density_field)

    for atom in tqdm(atoms, desc="Computing charge density"):
        atomic_data: AtomDescriptor = atom
        x, y, z = atomic_data.position
        charge = atomic_data.atomic_number  # Nucleus charge (in elementary charge units)

        # Compute squared distance
        r2 = (X - x)**2 + (Y - y)**2 + (Z - z)**2

        # Use an atom-specific alpha value
        alpha = 1.0 / (0.5 * atomic_data.nuclear_radius**2)  # Example scaling based on atomic radius
        
        # Nuclear charge distribution (normalized Gaussian)
        nucleus_density += charge * np.exp(-alpha * r2) * (alpha / np.pi)**(3/2)

    return nucleus_density - density_field  # Compute total charge density

def preprocess(filename):
    reader : GaussianCubeIO = GaussianCubeIO()
    reader.read(filename)
    atoms = reader.atoms
    density = reader.data
    n_atoms = reader.n_atoms
    origin = reader.origin
    npoints = reader.npoints * 2
    spacing = reader.spacing / 2

    x_range = (origin[0], origin[0] + npoints[0] * spacing[0])
    y_range = (origin[1], origin[1] + npoints[1] * spacing[1])
    z_range = (origin[2], origin[2] + npoints[2] * spacing[2])

    dx = spacing[0]
    dy = spacing[1]
    dz = spacing[2]


    print("Origin:", origin)
    print("Npoints:", npoints)
    print("Spacing:", spacing)

    density = polyatomic_molecule_density(atoms, x_range, y_range, z_range, mode='none', delta_x=dx, delta_y=dy, delta_z=dz, jobs=-1)

    charge = compute_charge_density(density, atoms, (origin[0], origin[0] + npoints[0] * spacing[0]), (origin[1], origin[1] + npoints[1] * spacing[1]), (origin[2], origin[2] + npoints[2] * spacing[2]), spacing[0], spacing[1], spacing[2])
    reader.data = charge
    reader.write_file(output_file)

def preprocess_pdb(filename, output_file):
    reader = PDBReader()
    reader.read_file(filename)
    atoms = reader.atoms
    n_atoms = len(atoms)
    origin = (min(atom.position[0] for atom in atoms), min(atom.position[1] for atom in atoms), min(atom.position[2] for atom in atoms))
    padding = 5.0
    origin = (origin[0] - padding, origin[1] - padding, origin[2] - padding)
    max_x = max(atom.position[0] for atom in atoms) + padding
    max_y = max(atom.position[1] for atom in atoms) + padding
    max_z = max(atom.position[2] for atom in atoms) + padding
    dx = 0.5
    dy = 0.5
    dz = 0.5

    npoints = (int((max_x - origin[0]) / dx), int((max_y - origin[1]) / dy), int((max_z - origin[2]) / dz))
    spacing = (dx, dy, dz)

    x_range = (origin[0], origin[0] + npoints[0] * spacing[0])
    y_range = (origin[1], origin[1] + npoints[1] * spacing[1])
    z_range = (origin[2], origin[2] + npoints[2] * spacing[2])

    density = polyatomic_molecule_density(atoms, x_range, y_range, z_range, mode='none', delta_x=dx, delta_y=dy, delta_z=dz, jobs=-1)
    charge = compute_charge_density(density, atoms, (origin[0], origin[0] + npoints[0] * spacing[0]), (origin[1], origin[1] + npoints[1] * spacing[1]), (origin[2], origin[2] + npoints[2] * spacing[2]), spacing[0], spacing[1], spacing[2])

    writer = GaussianCubeIO()
    writer.atoms = atoms
    writer.data = charge
    writer.n_atoms = n_atoms
    writer.npoints = npoints
    writer.origin = origin
    writer.spacing = spacing
    writer.write_file(output_file)



if __name__ == "__main__":

    filename = "data/h2o.pdb"
    output_file = "output/water/water_charge.cube"

    if not os.path.exists(output_file):
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        preprocess_pdb(filename, output_file)
    reader = GaussianCubeIO()
    reader.read(output_file)
    atoms = reader.atoms
    density = reader.data
    n_atoms = reader.n_atoms
    origin = reader.origin
    npoints = reader.npoints
    spacing = reader.spacing



    plotter = isodensity_plotter.IsodensityPlotter(
        (origin[0], origin[0] + npoints[0] * spacing[0]),
        (origin[1], origin[1] + npoints[1] * spacing[1]),
        (origin[2], origin[2] + npoints[2] * spacing[2]),
        spacing[0],
        spacing[1],
        spacing[2],
        atoms,
        os.path.join("output", "water")
    )
    plotter.plot_sliced_density_field_2D(density, 0, 'xz')