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
    if (!eddi_read_gaussian_cube(argv[1], &field, &molecule)) return EDDI_TEST_FAILURE;
    return !eddi_write_binary("data.bin", &field);

}