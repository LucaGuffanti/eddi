/**
 * @file eddi_molecule_readers.c
 * @author Luca Guffanti
 * @brief Implementation of the molecule readers
 */

#include "eddi_molecule_readers.h"



bool eddi_read_pdb(const char* filename, eddi_molecule_t* molecule)
{
    FILE* fp;

    if (!(fp = fopen(filename, "r")))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not open file\n");
        return EDDI_TEST_FAILURE;
    }

    EDDI_DEBUG_PRINT("[INFO] Reading pdb file\n");

    eddi_size_t count = 0;
    
    eddi_input_atom_list_t* atom_list = NULL;

    char type_of_record[10];
    fscanf(fp, "%[^' ']", type_of_record);
    while(!feof(fp))
    {
        // If either an ATOM or HETATM is read than manage it
        if((type_of_record[0] == 'A' && type_of_record[1] == 'T') || (type_of_record[0] == 'H' && type_of_record[3] == 'A'))
        {
            double x;
            double y;
            double z;
            char name[3];
            
            fscanf(fp, " %*d %*s %*s %*s %*d %lf %lf %lf %*lf %*lf %s", &x, &y, &z, name);
            // EDDI_DEBUG_PRINT("%lf %lf %lf %s\n", x, y, z, name);

            // PDB files provide atom information in Angstroms!
            x = x * EDDI_ANGSTROM_TO_BOHR;
            y = y * EDDI_ANGSTROM_TO_BOHR;
            z = z * EDDI_ANGSTROM_TO_BOHR;

            if (!eddi_add_atom_to_list_symbol(&atom_list, x, y, z, name))
            {
                EDDI_DEBUG_PRINT("[ERROR] Atom was not added to list\n");
                return EDDI_RETURN_FAILURE;
            }
            count++;
        }

        fscanf(fp, "%*[^\n]%*c");
        fscanf(fp, "%[^' ']", type_of_record);
    }
    
    fclose(fp);

    // Now transform the atomic list to a molecule object.
    eddi_atom_list_to_molecule(atom_list, molecule, count);
    eddi_free_input_atom_list(atom_list);
    
    return EDDI_RETURN_SUCCESS;
}

