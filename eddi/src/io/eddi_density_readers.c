/**
 * @file eddi_density_readers.c
 * @author Luca Guffanti
 * @brief Implementation of the density readers
 */

#include "eddi_density_readers.h"

bool eddi_read_binary(const char* filename, eddi_density_field_t* density_field)
{
    FILE* fp;
    
    fp = fopen(filename, "rb");
    if (!fp)
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not open file\n");
        return EDDI_RETURN_FAILURE;
    }

    eddi_real_t dx;
    eddi_real_t dy;
    eddi_real_t dz;

    eddi_point_t origin;

    eddi_size_t x_size;
    eddi_size_t y_size;
    eddi_size_t z_size;

    // ... First we read the three voxel dimensions: dx dy dz
    fread((void*) &(dx), sizeof(eddi_real_t), 1, fp);
    fread((void*) &(dy), sizeof(eddi_real_t), 1, fp);
    fread((void*) &(dz), sizeof(eddi_real_t), 1, fp);
    
    // ... Then the origin coordinates: origin_x, origin_y, origin_z
    fread((void*) &(origin), sizeof(eddi_point_t), 1, fp);
    
    // ... Then the number of number of points along each dimension
    fread((void*) &(x_size), sizeof(eddi_size_t), 1, fp);
    fread((void*) &(y_size), sizeof(eddi_size_t), 1, fp);
    fread((void*) &(z_size), sizeof(eddi_size_t), 1, fp);
    EDDI_DEBUG_PRINT("[INFO] Reading density field header\n");

    
    if (!eddi_new_density_field(density_field, dx, dy, dz, &origin, x_size, y_size, z_size))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not allocate density field\n");
        fclose(fp);
        return EDDI_RETURN_FAILURE;
    }
    
    // Now the density field object is constructed and all its internals are ready to contain the electron density values
    fread((void*) (density_field->field), sizeof(eddi_real_t), x_size * y_size * z_size, fp);
    EDDI_DEBUG_PRINT("[INFO] Density field was created\n");

    // Now the density field is ready for use: close the file and return success.
    fclose(fp);

    return EDDI_RETURN_SUCCESS;
}

bool eddi_read_gaussian_cube(const char* filename, eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    FILE* fp;

    fp = fopen(filename, "r");
    if (!fp)
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not open file\n");
        return EDDI_RETURN_FAILURE;
    }

    eddi_real_t dx;
    eddi_real_t dy;
    eddi_real_t dz;

    eddi_point_t origin;

    eddi_size_t n_atoms;

    int x_size;
    int y_size;
    int z_size;

    // The first two lines are comments
    fscanf(fp,"%*[^\n]%*c");
    fscanf(fp,"%*[^\n]%*c");

    // Then, we have the number of atoms and position of the origin
#ifdef EDDI_HIGH_PRECISION
    fscanf(fp, " %zu %lf %lf %lf", &n_atoms, &origin.x, &origin.y, &origin.z);
    EDDI_DEBUG_PRINT("Atoms: %ld  Origin: %f %f %f\n", n_atoms, origin.x, origin.y, origin.z);

    // Then, the number of points along each direction and the resolution
    fscanf(fp, " %d %lf %*lf %*lf", &x_size, &dx);
    fscanf(fp, " %d %*lf %lf %*lf", &y_size, &dy);
    fscanf(fp, " %d %*lf %*lf %lf", &z_size, &dz);
#else
    fscanf(fp, " %zu %f %f %f", &n_atoms, &origin.x, &origin.y, &origin.z);
    EDDI_DEBUG_PRINT("Atoms: %ld  Origin: %f %f %f\n", n_atoms, origin.x, origin.y, origin.z);

    // Then, the number of points along each direction and the resolution
    fscanf(fp, " %d %f %*f %*f", &x_size, &dx);
    fscanf(fp, " %d %*f %f %*f", &y_size, &dy);
    fscanf(fp, " %d %*f %*f %f", &z_size, &dz);
#endif

    // TODO: add support for angstrom to bohrs conversion.
    assert(x_size > 0 && y_size > 0 && z_size > 0 && "ANGSTROMS TO BOHRS CONVERSION NOT YET SUPPORTED");

    
    EDDI_DEBUG_PRINT("nx: %d dx: %f\n", x_size, dx);
    EDDI_DEBUG_PRINT("ny: %d dy: %f\n", y_size, dy);
    EDDI_DEBUG_PRINT("nz: %d dz: %f\n", z_size, dz);

    // Then a list of atoms
    
    eddi_atomic_number_t atomic_numbers[n_atoms];
    eddi_real_t x_positions[n_atoms];
    eddi_real_t y_positions[n_atoms];
    eddi_real_t z_positions[n_atoms];
    
    // TODO: add support for the charges if necessary
    for (eddi_size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx)
    {
#ifdef EDDI_HIGH_PRECISION
        fscanf(fp, " %hhd %*lf %lf %lf %lf\n", &atomic_numbers[atom_idx], &x_positions[atom_idx], &y_positions[atom_idx], &z_positions[atom_idx]);
#else
    fscanf(fp, " %hhd %*lf %f %f %f\n", &atomic_numbers[atom_idx], &x_positions[atom_idx], &y_positions[atom_idx], &z_positions[atom_idx]);
#endif
        EDDI_DEBUG_PRINT("Atom %zu: Atomic Number: %d, Position: (%f, %f, %f)\n", 
                 atom_idx, atomic_numbers[atom_idx], 
                 x_positions[atom_idx], y_positions[atom_idx], z_positions[atom_idx]);
    }


    if (!eddi_new_molecule(molecule, n_atoms, x_positions, y_positions, z_positions, atomic_numbers))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not create molecule\n");
        fclose(fp);
        return EDDI_RETURN_FAILURE;
    }

    EDDI_DEBUG_PRINT("[INFO] Creating Molecule\n");
    // Then, construct the density field and prepare it to hold the values

    if (!eddi_new_density_field(density_field, dx, dy, dz, &origin, x_size, y_size, z_size))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not allocate density field\n");
        fclose(fp);
        return EDDI_RETURN_FAILURE;
    }

    EDDI_DEBUG_PRINT("[INFO] Reading density field\n");
    eddi_size_t all_points = x_size * y_size * z_size;
    for (eddi_size_t idx = 0; idx < all_points; ++idx)
    {
        fscanf(fp, " %lE", &density_field->field[idx]);
    }


    fclose(fp);


    return EDDI_RETURN_SUCCESS;
}
