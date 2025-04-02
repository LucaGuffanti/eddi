from pdb_reader import PDBReader
from gaussian_cube_io import GaussianCubeIO
from general_molecule import polyatomic_molecule_density


if __name__ == "__main__":

    reader = PDBReader()
    reader.read_file("data/512_atoms.pdb")
    atoms = reader.atoms

    xmin = min(atom.position[0] for atom in atoms)
    xmax = max(atom.position[0] for atom in atoms)
    ymin = min(atom.position[1] for atom in atoms)
    ymax = max(atom.position[1] for atom in atoms)
    zmin = min(atom.position[2] for atom in atoms)
    zmax = max(atom.position[2] for atom in atoms)

    x_range = (xmin - 5, xmax + 5)
    y_range = (ymin - 5, ymax + 5)
    z_range = (zmin - 5, zmax + 5)

    print("X range:", x_range)
    print("Y range:", y_range)
    print("Z range:", z_range)
    
    dx = 1
    dy = 1
    dz = 1

    density = polyatomic_molecule_density('slater', atoms, x_range, y_range, z_range, mode='vdw', delta_x=dx, delta_y=dy, delta_z=dz, jobs=-1)
    
    writer = GaussianCubeIO()
    writer.atoms = atoms
    writer.data = density
    writer.n_atoms = len(atoms)
    writer.npoints = density.shape
    writer.origin = (x_range[0], y_range[0], z_range[0])
    writer.spacing = (dx, dy, dz)

    writer.write_file("output/512_atoms/512_atoms.cube")
    
    



