/**
 * @file eddi_density_writers.c
 * @author Luca Guffanti
 * @brief Implementation of the density writers
 */

#include "eddi_density_writers.h"

bool eddi_write_binary(const char* filename, const eddi_density_field_t* density_field)
{
    FILE* fp;

    if (!(fp = fopen(filename, "wb")))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not open file");
        return EDDI_RETURN_FAILURE;
    }
    
    // Start by writing the grid deltas
    fwrite((void*) &(density_field->dx), sizeof(eddi_real_t), 1, fp);
    fwrite((void*) &(density_field->dy), sizeof(eddi_real_t), 1, fp);
    fwrite((void*) &(density_field->dz), sizeof(eddi_real_t), 1, fp);

    // Then the origin
    fwrite((void*) &(density_field->origin.x), sizeof(eddi_real_t), 1, fp);
    fwrite((void*) &(density_field->origin.y), sizeof(eddi_real_t), 1, fp);
    fwrite((void*) &(density_field->origin.z), sizeof(eddi_real_t), 1, fp);


    // Then the number of elements along each dimension
    fwrite((void*) &(density_field->x_size), sizeof(eddi_size_t), 1, fp);
    fwrite((void*) &(density_field->y_size), sizeof(eddi_size_t), 1, fp);
    fwrite((void*) &(density_field->z_size), sizeof(eddi_size_t), 1, fp);

    // And finally the entire field
    fwrite((void*) (density_field->field), sizeof(eddi_real_t), density_field->x_size * density_field->y_size* density_field->z_size, fp);
    fclose(fp);

    return EDDI_RETURN_SUCCESS;
}

bool eddi_write_gaussian_cube(const char* filename, const eddi_density_field_t* density_field, const eddi_molecule_t* molecule)
{
    FILE* fp;
    if (!(fp = fopen(filename, "wb")))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not open file");
        return EDDI_RETURN_FAILURE;
    }

    // The first two lines are not used
    fprintf(fp, "Eddi Electron Density Distribution\n");
    fprintf(fp, "Outer: X, Middle: Y, Inner: Z\n");

    // The third line has the number of atoms and the origin of the domain
    fprintf(fp, "%ld %.6lf %.6lf %.6lf\n", molecule->n_atoms, density_field->origin.x, density_field->origin.y, density_field->origin.z);
    
    // In the following three lines there are the extensions of each dimension and the resolutions of the voxels
    fprintf(fp, "%zu %.6lf 0.000000 0.000000\n", density_field->x_size, density_field->dx);
    fprintf(fp, "%zu 0.000000 %.6lf 0.000000\n", density_field->y_size, density_field->dy);
    fprintf(fp, "%zu 0.000000 0.000000 %.6lf\n", density_field->z_size, density_field->dz);

    // Then, a list of atoms with atomic number, charge and position vector
    for (eddi_size_t idx = 0; idx < molecule->n_atoms; ++idx)
        fprintf(fp, "%d 0.000000 %.6lf %.6lf %.6lf\n", molecule->atomic_numbers[idx], molecule->atoms_x[idx], molecule->atoms_y[idx], molecule->atoms_z[idx]);

    const eddi_size_t nx = density_field->x_size;
    const eddi_size_t ny = density_field->y_size;
    const eddi_size_t nz = density_field->z_size;

    // Finally, the list of electron densities
    for (eddi_size_t ix = 0; ix < nx; ++ix)
    {
        for (eddi_size_t iy = 0; iy < ny; ++iy)
        {
            for (eddi_size_t iz = 0; iz < nz; ++iz)
            {
                fprintf(fp, "%.5lE ", density_field->field[ix * nz * ny + iy * nz + iz]);
                if (iz % 6 == 5)
                    fprintf(fp, "\n");
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
    return EDDI_RETURN_SUCCESS;
}
