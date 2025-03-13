import numpy as np
from gaussian_cube_io import GaussianCubeIO
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from joblib import Parallel, delayed
from tqdm import tqdm
from matplotlib import cm
from constants import BOHRS_PER_ANGSTROMS

def volume_computation(density_field: np.ndarray, isodensity, dx, dy, dz, accept_exact_isodensity=False):
    """Computes the volume of a given isodensity region in a density field. Space
    is subdivided in a cubic grid with each cube of size dx * dy * dz. The volume
    of the isodensity region is computed by counting the number 

    Args:
        density_field (_type_): 3D density field
        isodensity (_type): isodensity value
        dx (_type_): voxel dx
        dy (_type_): voxel dy
        dz (_type_): voxel dz
    """

    n_points_x = density_field.shape[0]
    n_points_y = density_field.shape[1]
    n_points_z = density_field.shape[2]

    n_cubes_x = n_points_x - 1
    n_cubes_y = n_points_y - 1
    n_cubes_z = n_points_z - 1

    volume = 0

    def compute_volume_slice(i, accept_exact_isodensity=False):
        volume = 0
        for j in range(n_cubes_y):
            for z in range(n_cubes_z):
                if accept_exact_isodensity:
                    if (density_field[i, j, z] >= isodensity and
                        density_field[i + 1, j, z] >= isodensity and
                        density_field[i, j + 1, z] >= isodensity and
                        density_field[i + 1, j + 1, z] >= isodensity and
                        density_field[i, j, z + 1] >= isodensity and
                        density_field[i + 1, j, z + 1] >= isodensity and
                        density_field[i, j + 1, z + 1] >= isodensity and
                        density_field[i + 1, j + 1, z + 1] >= isodensity):
                            volume += 1
                else:
                    if (density_field[i, j, z] > isodensity and
                        density_field[i + 1, j, z] > isodensity and
                        density_field[i, j + 1, z] > isodensity and
                        density_field[i + 1, j + 1, z] > isodensity and
                        density_field[i, j, z + 1] > isodensity and
                        density_field[i + 1, j, z + 1] > isodensity and
                        density_field[i, j + 1, z + 1] > isodensity and
                        density_field[i + 1, j + 1, z + 1] > isodensity):
                            volume += 1
        return volume

    volume = sum(Parallel(n_jobs=-1)(delayed(compute_volume_slice)(i, accept_exact_isodensity) for i in tqdm(range(n_cubes_x), desc="Processing slices")))
    return volume * dx * dy * dz

def represent_volume_units(density_field, isodensity, dx, dy, dz, accept_exact_isodensity=False):
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    n_points_x = density_field.shape[0]
    n_points_y = density_field.shape[1]
    n_points_z = density_field.shape[2]

    n_cubes_x = n_points_x - 1
    n_cubes_y = n_points_y - 1
    n_cubes_z = n_points_z - 1

    def process_cube(i, accept_exact_isodensity=False):
        verts=[]
        for j in range(n_cubes_y):
            for k in range(n_cubes_z):
                if accept_exact_isodensity:
                    if (density_field[i, j, k] >= isodensity and
                        density_field[i + 1, j, k] >= isodensity and
                        density_field[i, j + 1, k] >= isodensity and
                        density_field[i + 1, j + 1, k] >= isodensity and
                        density_field[i, j, k + 1] >= isodensity and
                        density_field[i + 1, j, k + 1] >= isodensity and
                        density_field[i, j + 1, k + 1] >= isodensity and
                        density_field[i + 1, j + 1, k + 1] >= isodensity):
                        
                        # Check if the cube is on the surface
                        if (i == 0 or i == n_cubes_x - 1 or
                            j == 0 or j == n_cubes_y - 1 or
                            k == 0 or k == n_cubes_z - 1 or
                            density_field[i - 1, j, k] < isodensity or
                            density_field[i + 2, j, k] < isodensity or
                            density_field[i, j - 1, k] < isodensity or
                            density_field[i, j + 2, k] < isodensity or
                            density_field[i, j, k - 1] < isodensity or
                            density_field[i, j, k + 2] < isodensity):
                            
                            # Define the vertices of the cube
                            vertices = [
                                [i * dx, j * dy, k * dz],
                                [(i + 1) * dx, j * dy, k * dz],
                                [(i + 1) * dx, (j + 1) * dy, k * dz],
                                [i * dx, (j + 1) * dy, k * dz],
                                [i * dx, j * dy, (k + 1) * dz],
                                [(i + 1) * dx, j * dy, (k + 1) * dz],
                                [(i + 1) * dx, (j + 1) * dy, (k + 1) * dz],
                                [i * dx, (j + 1) * dy, (k + 1) * dz]
                            ]
                            verts.append(vertices)
                else:
                    if (density_field[i, j, k] > isodensity and
                        density_field[i + 1, j, k] > isodensity and
                        density_field[i, j + 1, k] > isodensity and
                        density_field[i + 1, j + 1, k] > isodensity and
                        density_field[i, j, k + 1] > isodensity and
                        density_field[i + 1, j, k + 1] > isodensity and
                        density_field[i, j + 1, k + 1] > isodensity and
                        density_field[i + 1, j + 1, k + 1] > isodensity):
                        
                        # Check if the cube is on the surface
                        if (i == 0 or i == n_cubes_x - 1 or
                            j == 0 or j == n_cubes_y - 1 or
                            k == 0 or k == n_cubes_z - 1 or
                            density_field[i - 1, j, k] < isodensity or
                            density_field[i + 2, j, k] < isodensity or
                            density_field[i, j - 1, k] < isodensity or
                            density_field[i, j + 2, k] < isodensity or
                            density_field[i, j, k - 1] < isodensity or
                            density_field[i, j, k + 2] < isodensity):
                            
                            # Define the vertices of the cube
                            vertices = [
                                [i * dx, j * dy, k * dz],
                                [(i + 1) * dx, j * dy, k * dz],
                                [(i + 1) * dx, (j + 1) * dy, k * dz],
                                [i * dx, (j + 1) * dy, k * dz],
                                [i * dx, j * dy, (k + 1) * dz],
                                [(i + 1) * dx, j * dy, (k + 1) * dz],
                                [(i + 1) * dx, (j + 1) * dy, (k + 1) * dz],
                                [i * dx, (j + 1) * dy, (k + 1) * dz]
                            ]
                            verts.append(vertices)
        if verts:
            return verts
        return None 

    results = Parallel(n_jobs=-1)(
        delayed(process_cube)(i, accept_exact_isodensity) 
        for i in tqdm(range(n_cubes_x), desc="Processing X-axis") 
        
    )

    res = []
    for r in results:
        if r:
            res.extend(r)
    results = res

    unique_vertices = set()
    for vertices in tqdm(results, desc="Processing cubes"):
        if vertices:
            for vertex in vertices:
                unique_vertices.add(tuple(vertex))
    
    unique_vertices = np.array(list(unique_vertices))
    # Normalize the z values for the color mapping
    norm = plt.Normalize(unique_vertices[:, 2].min(), unique_vertices[:, 2].max())
    colors = cm.viridis(norm(unique_vertices[:, 2]))

    print("Total points", len(unique_vertices))

    ax.scatter(unique_vertices[:, 0], unique_vertices[:, 1], unique_vertices[:, 2], color=colors, alpha=0.6, s=0.05)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

def volume_computation_angstroms(density_field, isodensity, dx, dy, dz, accept_exact_isodensity=False):
    return volume_computation(density_field, isodensity, dx, dy, dz, accept_exact_isodensity) / (BOHRS_PER_ANGSTROMS ** 3)

def compute_isodensity_for_volume(density_field, target_volume, dx, dy, dz, mode='angstroms', accept_exact_isodensity=False):
    starting_volume = volume_computation(density_field, 0.1, dx, dy, dz)
    max_it = 100
    it = 0
    eps = 1e-3
    low_volume_isodensity = max(density_field.flatten())
    high_volume_isodensity = 0.0

    while it < max_it and np.abs(starting_volume - target_volume)/target_volume > eps:
        
        isodensity = (low_volume_isodensity + high_volume_isodensity) / 2
        volume = 0
        if mode == 'angstrom':
            volume = volume_computation_angstroms(density_field, isodensity, dx, dy, dz, accept_exact_isodensity)
        else:
            volume = volume_computation(density_field, isodensity, dx, dy, dz, accept_exact_isodensity)
        

        # if the computed volume is too big than the target, we need to get to a smaller volume, which is done by 
        # increasing the isodensity of the computation
        if volume > target_volume:
            high_volume_isodensity = isodensity
        elif volume < target_volume:
            low_volume_isodensity = isodensity
        else:
            break
        print("it:", it+1, "isodensity:", isodensity, "volume:", volume, "target_volume:", target_volume)
        it += 1

    return isodensity

if __name__ == "__main__":

    reader = GaussianCubeIO()
    reader.read('output/alanine/alanine.cube')

    represent_volume_units(reader.data, 0.00274, reader.spacing[0], reader.spacing[1], reader.spacing[2], accept_exact_isodensity=False)