#include "eddi_molecule.h"
#include "eddi_density_field.h"
#include "eddi_density_writers.h"
#include "eddi_molecule_readers.h"

// Doing a simple simulation of an Hydrogren molecule centered in 0.0
int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t density_field;

    eddi_read_pdb(argv[1], &molecule);
    
    eddi_init_field_from_molecule(&density_field, &molecule, 1.0, 0.5, 0.5, 0.5);
    eddi_compute_density_field(&density_field, &molecule);
    
    eddi_write_gaussian_cube(argv[2], &density_field, &molecule);
    eddi_write_binary(argv[3], &density_field);

    printf("Free data\n");
    eddi_free_molecule(&molecule);
    eddi_free_density_field(&density_field);
    eddi_free_atomic_data();

    return EDDI_TEST_SUCCESS;
    
}