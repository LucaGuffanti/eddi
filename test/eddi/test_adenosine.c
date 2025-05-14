#include "eddi_molecule_readers.h"
#include "eddi_density_writers.h"

int main(int argc, char** argv)
{

    eddi_molecule_t molecule;
    eddi_density_field_t density_field;
    
    puts("Reading molecule from file");
    eddi_read_mol(argv[1], &molecule);
    // eddi_print_molecule(&molecule);
    puts("Initializing field");
    bool res1 = !eddi_init_field_from_molecule(&density_field, &molecule, 5.0, 1.0, 1.0, 1.0);

    puts("Computing density");
    eddi_compute_density_field(&density_field, &molecule);

    puts("Printing to gaussian");
    bool res2 = !eddi_write_gaussian_cube(argv[2], &density_field, &molecule);

    eddi_write_binary(argv[3], &density_field);

    eddi_free_atomic_data();
    eddi_free_molecule(&molecule);
    eddi_free_density_field(&density_field);

    return res1 || res2;
}