#include "eddi_molecule.h"
#include "eddi_density_field.h"
#include "eddi_density_writers.h"

// Doing a simple simulation of an He molecule centered in 0.0
int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_real_t x[] = {0.0, 0.5};
    eddi_real_t y[] = {0.0, 0.0};
    eddi_real_t z[] = {0.0, 0,0};
    eddi_atomic_number_t atomic_numbers[] = {2, 2};

    eddi_new_molecule(&molecule, 2, x, y, z, atomic_numbers);

    eddi_density_field_t density_field;

    eddi_init_field_from_molecule(&density_field, &molecule, 2.0, 0.1, 0.1, 0.1);
    
    eddi_compute_density_field(&density_field, &molecule);
    
    eddi_write_gaussian_cube(argv[1], &density_field, &molecule);
    eddi_write_binary(argv[2], &density_field);

    eddi_free_molecule(&molecule);
    eddi_free_density_field(&density_field);

    return EDDI_TEST_SUCCESS;
    
}