import numpy as np
from atom_descriptor import AtomDescriptor
import tqdm

class GaussianCubeIO:

    def __init__ (self):
        self.atoms = []
        self.n_atoms = 0
        self.origin = None
        self.npoints = None
        self.spacing = None
        self.filename = None
        self.current_line = 0

    def read (self, filename):
        self.filename = filename
        
        print('Reading file: ', filename)
        with open(filename, 'r') as f:
            lines = f.readlines()
            self.read_header(lines)
            self.read_data(lines)

        print('Done reading file: ', filename)
    
    def read_header(self, lines):
        print("Reading header")
        self.current_line = 2
        line_values = lines[self.current_line].strip().split()
        self.n_atoms = int(line_values[0])
        self.origin = np.array([float(line_values[1]), float(line_values[2]), float(line_values[3])])

        print('Number of atoms: ', self.n_atoms)
        print('Origin: ', self.origin)

        self.current_line += 1
        line_values = [float(i) for i in lines[self.current_line].strip().split()]
        assert len(line_values) == 4 and line_values[2] == float(0) and line_values[3] == float(0)
        x_points = int(line_values[0])
        dx = float(line_values[1])

        self.current_line += 1
        line_values = [float(i) for i in lines[self.current_line].strip().split()]
        assert len(line_values) == 4 and line_values[1] == float(0) and line_values[3] == float(0)
        y_points = int(line_values[0])
        dy = float(line_values[2])

        self.current_line += 1
        line_values = [float(i) for i in lines[self.current_line].strip().split()]
        assert len(line_values) == 4 and line_values[1] == float(0) and line_values[2] == float(0)
        z_points = int(line_values[0])
        dz = float(line_values[3])

        self.npoints = np.array([x_points, y_points, z_points])
        self.spacing = np.array([dx, dy, dz])
        print('Number of points: ', self.npoints)
        self.total_points = x_points * y_points * z_points
        print('Spacing: ', self.spacing)

        self.current_line += 1
        for i in range(self.n_atoms):
            line_values = lines[self.current_line].strip().split()
            atomic_number = int(line_values[0])
            charge = float(line_values[1])
            x = float(line_values[2])
            y = float(line_values[3])
            z = float(line_values[4])
            self.atoms.append(AtomDescriptor(atomic_number, (x, y, z), charge))
            self.current_line += 1

        for atom in self.atoms:
            print(' Atomic number: ', atom.atomic_number)
            print(' Position: ', atom.position)
            print(' Charge: ', atom.charge)


    def read_data(self, lines):
        print("Reading data")
        
        self.data = np.zeros(self.total_points)
        counter = 0
        for i in tqdm.tqdm(range(self.current_line, len(lines))):
            line_values = [float(k) for k in lines[self.current_line].strip().split()]
            for val in line_values:
                self.data[counter] = val
                counter = counter + 1
            self.current_line += 1

        print('Data shape: ', self.data.shape)
        self.data.resize(self.npoints[0], self.npoints[1], self.npoints[2])
        print('Data shape: ', self.data.shape)

    def write_file(self, filename):
        print("Writing file:", filename)
        with open(filename, 'w') as f:
            self._write_header(f)
            self._write_data(f)
        print("Done writing file:", filename)

    def _write_header(self, f):
        f.write(" cube creation\n SCF Total Density\n")
        f.write(f" {self.n_atoms} {self.origin[0]:.6f} {self.origin[1]:.6f} {self.origin[2]:.6f}\n")
        f.write(f" {self.npoints[0]} {self.spacing[0]:.6f} 0.000000 0.000000\n")
        f.write(f" {self.npoints[1]} 0.000000 {self.spacing[1]:.6f} 0.000000\n")
        f.write(f" {self.npoints[2]} 0.000000 0.000000 {self.spacing[2]:.6f}\n")
        for atom in self.atoms:
            f.write(f" {atom.atomic_number} {atom.charge} {atom.position[0]:.6f} {atom.position[1]:.6f} {atom.position[2]:.6f}\n")

    def _write_data(self, f):
        x_points = self.npoints[0]
        y_points = self.npoints[1]
        z_points = self.npoints[2]

        for i in range(x_points):
            for j in range(y_points):
                for k in range(z_points):
                    f.write(f"{self.data[i, j, k]:.5E} ")
                    if k % 6 == 5:
                        f.write("\n")
                f.write("\n")

def main():
    reader = GaussianCubeIO()
    reader.read('data/o2_3.554Angstr/o2.cube')
    reader.write_file('output/my_file.cube')
    reader.read('output/my_file.cube')

if __name__ == '__main__':
    main()