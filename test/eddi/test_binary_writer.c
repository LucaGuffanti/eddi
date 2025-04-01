#include "stdio.h"
#include "eddi_density_readers.h"
#include "eddi_density_field.h"
#include "eddi_density_writers.h"
#include "eddi_molecule.h"

int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t field;

    assert(argc == 2);

    bool result1 = eddi_read_gaussian_cube(argv[1], &field, &molecule);

    eddi_free_molecule(&molecule);
    bool result2;
    if (result1)
        result2 = !eddi_write_binary("data.bin", &field);

    eddi_free_density_field(&field);
    return result2;
}