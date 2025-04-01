#include "eddi_molecule_readers.h"

int main(int argc, char** argv)
{
    eddi_molecule_t molecule;

    bool res = !eddi_read_mol(argv[1], &molecule);
    
    eddi_free_molecule(&molecule);
    
    return res;
}
