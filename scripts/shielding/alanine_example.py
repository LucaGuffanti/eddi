from atom_descriptor import AtomDescriptor
from pdb_reader import PDBReader
from general_molecule import polyatomic_molecule_density
from plotting.isodensity_plotter import IsodensityPlotter
from gaussian_cube_io import GaussianCubeIO

reader = PDBReader()
# change with actual path
reader.read_file('PATH_TO_ALANINE')

atoms = reader.atoms

x_min = min(atom.position[0] for atom in atoms)
x_max = max(atom.position[0] for atom in atoms)

y_min = min(atom.position[1] for atom in atoms)
y_max = max(atom.position[1] for atom in atoms)

z_min = min(atom.position[2] for atom in atoms)
z_max = max(atom.position[2] for atom in atoms)

delta_x = 0.1
delta_y = 0.1
delta_z = 0.1

print(f"x_min: {x_min}, x_max: {x_max}")
print(f"y_min: {y_min}, y_max: {y_max}")
print(f"z_min: {z_min}, z_max: {z_max}")

x_range = [x_min - 2, x_max + 2]
y_range = [y_min - 2, y_max + 2]
z_range = [z_min - 4, z_max + 4]

density_field = polyatomic_molecule_density(atoms, x_range, y_range, z_range, 'none', delta_x, delta_y, delta_z)

plotter = IsodensityPlotter(
    x_range, y_range, z_range, delta_x, delta_y, delta_z, atoms, 'output/alanine'
)

# plotter.plot_isodensity_3D(density_field, 0.1)

cube_io = GaussianCubeIO()
cube_io.spacing = (delta_x, delta_y, delta_z)
cube_io.origin = (x_range[0], y_range[0], z_range[0])
cube_io.npoints = (len(density_field), len(density_field[0]), len(density_field[0][0]))
cube_io.data = density_field
cube_io.n_atoms = len(atoms)
cube_io.atoms = atoms
cube_io.write_file('output/alanine/alanine.cube')
