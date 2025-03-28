/**
 * @file eddy_density_field.c
 * @author Luca Guffanti
 * @brief Implementation of the density field 
 */

#include "eddi_density_field.h"

bool eddi_new_density_field(
    eddi_density_field_t* density_field, 
    const eddi_real_t dx,
    const eddi_real_t dy, 
    const eddi_real_t dz,
    const eddi_point_t* origin,
    const eddi_size_t n_x,
    const eddi_size_t n_y,
    const eddi_size_t n_z
)
{

    density_field->field = (eddi_array_t) malloc(sizeof(eddi_real_t) * n_x * n_y * n_z); 
    if (!density_field->field)
    {
        EDDI_DEBUG_PRINT("[ERROR] Density field object instatiation failed (internal field)\n");
        return EDDI_RETURN_FAILURE;
    }
    memset(density_field->field, 0, sizeof(eddi_real_t) * n_x * n_y * n_z);

    density_field->dx = dx;
    density_field->dy = dy;
    density_field->dz = dz;
    density_field->x_size = n_x;
    density_field->y_size = n_y;
    density_field->z_size = n_z;
    memcpy(&(density_field->origin), origin, 3 * sizeof(eddi_real_t));

    return EDDI_RETURN_SUCCESS; // Allocation successful
}

bool eddi_init_field_from_molecule(
    eddi_density_field_t* density_field, 
    const eddi_molecule_t* molecule,
    const eddi_real_t padding,
    const eddi_real_t dx,
    const eddi_real_t dy,
    const eddi_real_t dz
)
{

    // Extract the mininum and maximum positions of the atoms along each direction
    eddi_real_t min_x = molecule->atoms_x[0];
    eddi_real_t max_x = molecule->atoms_x[0];

    eddi_real_t min_y = molecule->atoms_y[0];
    eddi_real_t max_y = molecule->atoms_y[0];

    eddi_real_t min_z = molecule->atoms_z[0];
    eddi_real_t max_z = molecule->atoms_z[0];

    for (eddi_size_t i = 1; i < molecule->n_atoms; ++i)
    {
        if (molecule->atoms_x[i] < min_x) min_x = molecule->atoms_x[i];
        if (molecule->atoms_x[i] > max_x) max_x = molecule->atoms_x[i];

        if (molecule->atoms_y[i] < min_y) min_y = molecule->atoms_y[i];
        if (molecule->atoms_y[i] > max_y) max_y = molecule->atoms_y[i];

        if (molecule->atoms_z[i] < min_z) min_z = molecule->atoms_z[i];
        if (molecule->atoms_z[i] > max_z) max_z = molecule->atoms_z[i];
    }

    // Add padding to the min and max values
    min_x -= padding;
    max_x += padding;

    min_y -= padding;
    max_y += padding;

    min_z -= padding;
    max_z += padding;

    // Calculate the number of grid points along each direction
    eddi_size_t n_x = (eddi_size_t)ceil((max_x - min_x) / dx);
    eddi_size_t n_y = (eddi_size_t)ceil((max_y - min_y) / dy);
    eddi_size_t n_z = (eddi_size_t)ceil((max_z - min_z) / dz);

    // Set the origin point
    eddi_point_t origin = {min_x, min_y, min_z};

    // Initialize the density field
    return eddi_new_density_field(density_field, dx, dy, dz, &origin, n_x, n_y, n_z);
}


void eddi_compute_density_field(eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    const eddi_real_t dx = density_field->dx;
    const eddi_real_t dy = density_field->dy;
    const eddi_real_t dz = density_field->dz;

    const eddi_real_t origin_x = density_field->origin.x; 
    const eddi_real_t origin_y = density_field->origin.y;
    const eddi_real_t origin_z = density_field->origin.z;

    const eddi_size_t nx = density_field->x_size;
    const eddi_size_t ny = density_field->y_size;
    const eddi_size_t nz = density_field->z_size;

    EDDI_DEBUG_PRINT("dx: %f, dy: %f, dz: %f\n", dx, dy, dz);
    EDDI_DEBUG_PRINT("origin_x: %f, origin_y: %f, origin_z: %f\n", origin_x, origin_y, origin_z);
    EDDI_DEBUG_PRINT("nx: %zu, ny: %zu, nz: %zu\n", nx, ny, nz);

    const eddi_size_t atoms = molecule->n_atoms;

// #ifdef _OPENMP
//     #pragma omp parallel for collapse(3) shared(density_field, molecule)
// #else
// #endif
    // for (eddi_size_t it_x = 0; it_x < nx; ++it_x)
    // {
    //     for (eddi_size_t it_y = 0; it_y < ny; ++it_y)
    //     {
    //         for (eddi_size_t it_z = 0; it_z < nz; ++it_z)
    //         {
    //             eddi_real_t local_density = 0.0;
    //             eddi_real_t radius;
    //             for (eddi_size_t atom_idx = 0; atom_idx < atoms; ++atom_idx)
    //             {
    //                 const eddi_real_t cx = origin_x + it_x * dx;
    //                 const eddi_real_t cy = origin_y + it_y * dy;
    //                 const eddi_real_t cz = origin_z + it_z * dz;
                    
    //                 const eddi_real_t DeltaX = cx - molecule->atoms_x[atom_idx];
    //                 const eddi_real_t DeltaY = cy - molecule->atoms_y[atom_idx];
    //                 const eddi_real_t DeltaZ = cz - molecule->atoms_z[atom_idx];
    //                 radius = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);
                    
    //                 if (radius < dx) 
    //                 {
    //                     local_density += 0.0;
    //                 } 
    //                 else
    //                 {
    //                     local_density += molecule->density[atom_idx](radius, 0.0, 0.0);
    //                 }

    //             }
    //             density_field->field[it_x * ny * nz + it_y * nz + it_z] = local_density;
    //         }
    //     }
    // }

    const eddi_real_t dx_2 = dx * dx;

    eddi_array_t x = molecule->atoms_x;
    eddi_array_t y = molecule->atoms_y;
    eddi_array_t z = molecule->atoms_z;

    eddi_array_t field = density_field->field;

    const eddi_size_t n_pts = nx * ny * nz;
    #pragma omp parallel for 
    for (eddi_size_t idx = 0; idx < n_pts; ++idx)
    {
        const eddi_size_t it_z = idx % nz;
        const eddi_size_t it_y = (idx / nz) % ny;
        const eddi_size_t it_x = (idx / nz) / ny;

        const eddi_real_t cx = origin_x + it_x * dx;
        const eddi_real_t cy = origin_y + it_y * dy;
        const eddi_real_t cz = origin_z + it_z * dz;

        eddi_real_t local_density = 0.0;
        for (eddi_size_t atom_idx = 0; atom_idx < atoms; ++atom_idx)
        {
            
            const eddi_real_t DeltaX = cx - x[atom_idx];
            const eddi_real_t DeltaY = cy - y[atom_idx];
            const eddi_real_t DeltaZ = cz - z[atom_idx];
            const eddi_real_t radius_2 = DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ;
            
            if (radius_2 < dx_2) 
            {
                continue;
            } 
            else
            {
                local_density += molecule->density[atom_idx](sqrt(radius_2), 0.0, 0.0);
            }

        }
        field[idx] = local_density;
    }
}


void eddi_free_density_field(eddi_density_field_t* density_field)
{
    // Only need to deallocate the field
    free(density_field->field);
}