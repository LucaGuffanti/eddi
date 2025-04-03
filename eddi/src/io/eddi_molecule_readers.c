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
    fscanf(fp, "%6c", type_of_record);
    type_of_record[6] = '\0';

    while(!feof(fp))
    {
        // If either an ATOM or HETATM is read than manage it
        if((type_of_record[0] == 'A' && type_of_record[1] == 'T') /*|| (type_of_record[0] == 'H' && type_of_record[3] == 'A')*/)
        {
            eddi_real_t x;
            eddi_real_t y;
            eddi_real_t z;
            char name[3];
            
            // atom serial number - 5 int
            // atom name - 4 char
            // alternate location indicator - 1 char
            // residue name - 4 char
            // chain identifies - 1 char 
            // residue sequence number - 4 int
            // code for insertion of residues - 1 char
            // x - 8.3 real
            // y - 8.3 real
            // z - 8.3 real
            // occupancy  - 6.2 real
            // temperature factor - 6.2 real
            // segment identifier - 4 char
            
            char x_s[9];
            char y_s[9];
            char z_s[9];
            fscanf(fp, "%*5d%*c%*4c%*c%*3c%*c%*c%*4c%*c%*c%*c%*c%8c%8c%8c%*6c%*6c%*10c%2s", x_s, y_s, z_s, name);
        
            sscanf(x_s, "%lf", &x);
            sscanf(y_s, "%lf", &y);
            sscanf(z_s, "%lf", &z);

            EDDI_DEBUG_PRINT("%lf ", x);
            EDDI_DEBUG_PRINT("%lf ", y);
            EDDI_DEBUG_PRINT("%lf ", z);
            
            EDDI_DEBUG_PRINT("%s\n", name);

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
            EDDI_DEBUG_PRINT("%d\n", count);
        }

        fscanf(fp, "%*[^\n]%*c");
        fscanf(fp, "%6c", type_of_record);
    }

    fclose(fp);

    // Now transform the atomic list to a molecule object.
    eddi_atom_list_to_molecule(atom_list, molecule, count);
    eddi_free_input_atom_list(atom_list);
    
    return EDDI_RETURN_SUCCESS;
}

bool eddi_read_mol(const char* filename, eddi_molecule_t* molecule)
{
    FILE *fp;

    if (!(fp = fopen(filename, "r")))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not open file\n");
        return EDDI_RETURN_FAILURE;
    }

    fscanf(fp, "%*[^\n]%*c");
    fscanf(fp, "%*[^\n]%*c");
    fscanf(fp, "%*[^\n]%*c");

    eddi_size_t counts = 0;

    fscanf(fp, " %zu%*[^\n]%*c", &counts);

    eddi_real_t x[counts];
    eddi_real_t y[counts];
    eddi_real_t z[counts];
    eddi_atomic_number_t atomic_numbers[counts];

    char symbol[3];

    
    for (eddi_size_t i = 0; i < counts; ++i)
    {
        fscanf(fp, " %lf %lf %lf %[^' '] %*[^\n]%*c", &x[i], &y[i], &z[i], symbol);
        EDDI_DEBUG_PRINT("[INFO] Atom %zu: %lf %lf %lf %s\n", i + 1, x[i], y[i], z[i], symbol);

        x[i] *= EDDI_ANGSTROM_TO_BOHR; 
        y[i] *= EDDI_ANGSTROM_TO_BOHR;
        z[i] *= EDDI_ANGSTROM_TO_BOHR;

        atomic_numbers[i] = eddi_symbol_to_number(symbol);
    }

    fclose(fp);

    return eddi_new_molecule(molecule, counts, x, y, z, atomic_numbers);
}