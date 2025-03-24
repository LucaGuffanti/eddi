#include "eddi_molecule_readers.h"




int main(int argc, char** argv)
{
    eddi_molecule_t molecule;

    assert(argc == 2);
    bool final_result = !eddi_read_pdb(argv[1], &molecule);

    EDDI_DEBUG_CALL(eddi_print_molecule(&molecule));
    assert(final_result == EDDI_TEST_SUCCESS);

    return EDDI_TEST_SUCCESS;
}