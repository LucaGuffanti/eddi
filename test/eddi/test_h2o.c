#include "eddi_molecule.h"
#include "eddi_density_field.h"
#include "eddi_density_writers.h"
#include "eddi_molecule_readers.h"

// Doing a simple simulation of an He molecule centered in 0.0
int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t density_field;

    eddi_read_pdb(argv[1], &molecule);
    eddi_init_field_from_molecule(&density_field, &molecule, 6.0, 0.05, 0.05, 0.05);

    eddi_cl_info_t info = {.cx = 10.0, .cy = 10.0, .cz = 10.0};
    eddi_compute_density_field_cl(&density_field, &molecule, &info);
    
    eddi_write_gaussian_cube(argv[2], &density_field, &molecule);
    eddi_write_binary(argv[3], &density_field);

    eddi_free_molecule(&molecule);
    eddi_free_density_field(&density_field);
    eddi_free_atomic_data();

    return EDDI_TEST_SUCCESS;
    
}