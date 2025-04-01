#include "eddi_molecule_readers.h"




int main(int argc, char** argv)
{
    assert(argc == 2);
 
    eddi_molecule_t molecule;
    bool final_result = !eddi_read_pdb(argv[1], &molecule);

    eddi_free_molecule(&molecule);
    eddi_free_atomic_data();

    return EDDI_TEST_SUCCESS;
}