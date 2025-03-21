#include "stdio.h"
#include "eddi_density_readers.h"
#include "eddi_density_field.h"
#include "eddi_density_writers.h"
#include "eddi_molecule.h"

int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t field;


    printf("Argc: %d\n", argc);
    assert(argc == 5);

    // Start by reading a gaussian cube
    if (!eddi_read_gaussian_cube(argv[1], &field, &molecule))
        return EDDI_TEST_FAILURE;
    
    // Rewrite it to make sure that the implementation is correct
    if (!eddi_write_gaussian_cube(argv[2], &field, &molecule))
        return EDDI_TEST_FAILURE;
    
    // Construct the binary file from the gaussian cube
    if (!eddi_write_binary(argv[3], &field))
        return EDDI_TEST_FAILURE;

    // Read the just written binary file
    if (!eddi_read_binary(argv[3], &field))
        return EDDI_TEST_FAILURE;

    // And output it as another gaussian cube to check for correctness
    if (!eddi_write_gaussian_cube(argv[4], &field, &molecule))
        return EDDI_TEST_FAILURE;

    return EDDI_TEST_SUCCESS;

}