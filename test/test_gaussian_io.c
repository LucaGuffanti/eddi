#include "stdio.h"
#include "eddi_density_readers.h"
#include "eddi_density_field.h"
#include "eddi_density_writers.h"
#include "eddi_molecule.h"

int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t field;

    assert(argc == 3);
    if (!eddi_read_gaussian_cube(argv[1], &field, &molecule))
    {
        eddi_free_density_field(&field);
        eddi_free_molecule(&molecule);
        return EDDI_TEST_FAILURE;
    }  

    bool res = eddi_write_gaussian_cube(argv[2], &field, &molecule);
    eddi_free_density_field(&field);
    eddi_free_molecule(&molecule);

    return !res;

}