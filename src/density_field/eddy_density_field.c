/**
 * @file eddy_density_field.c
 * @author Luca Guffanti
 * @brief Implementation of the density field 
 */

#include "eddi_density_field.h"



bool eddi_new_density_field(eddi_density_field_t* density_field, 
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

    printf("dx: %f, dy: %f, dz: %f\n", dx, dy, dz);
    printf("origin_x: %f, origin_y: %f, origin_z: %f\n", origin_x, origin_y, origin_z);
    printf("nx: %zu, ny: %zu, nz: %zu\n", nx, ny, nz);

    const eddi_size_t atoms = molecule->n_atoms;

    eddi_real_t cx ;
    eddi_real_t cy;
    eddi_real_t cz;
    
    const double target_iso = 0.106329;

    cx = origin_x;
    for (eddi_size_t it_x = 0; it_x < nx; ++it_x)
    {
        cy = origin_y;
        for (eddi_size_t it_y = 0; it_y < ny; ++it_y)
        {
            cz = origin_z;
            for (eddi_size_t it_z = 0; it_z < nz; ++it_z)
            {
                eddi_real_t local_density = 0.0;
                eddi_real_t radius;
                for (eddi_size_t atom_idx = 0; atom_idx < atoms; ++atom_idx)
                {
                    const eddi_real_t DeltaX = cx - molecule->atoms_x[atom_idx];
                    const eddi_real_t DeltaY = cy - molecule->atoms_y[atom_idx];
                    const eddi_real_t DeltaZ = cz - molecule->atoms_z[atom_idx];
                    radius = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);
                    
                    if (radius < dx) 
                    {
                        local_density += 0.0;
                    } 
                    else
                    {
                        local_density += molecule->density[atom_idx](radius, 0.0, 0.0);
                    }

                }
                density_field->field[it_x * ny * nz + it_y * nz + it_z] = local_density;
                
                cz += dz;
            }
            cy += dy;
        }
        cx += dx;
    }
}

void eddi_free_density_field(eddi_density_field_t* density_field)
{
    // Only need to deallocate the field
    free(density_field->field);
    density_field->field = NULL;
}