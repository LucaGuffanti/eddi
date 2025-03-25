/**
 * @file eddi_input_atom_list.c
 * @author Luca Guffanti
 * @brief Implements the input atom list as defined in the header file
 */

#include "eddi_input_atom_list.h"

bool eddi_add_atom_to_list_symbol(eddi_input_atom_list_t** list, const eddi_real_t x, const eddi_real_t y, const eddi_real_t z, const char symbol[3])
{
    return eddi_add_atom_to_list_atomic_number(list, x, y, z, eddi_symbol_to_number(symbol));
}

bool eddi_add_atom_to_list_atomic_number(eddi_input_atom_list_t** list, const eddi_real_t x, const eddi_real_t y, const eddi_real_t z, const eddi_atomic_number_t atomic_number)
{
    eddi_input_atom_list_t* atom;
    if (!(atom = (eddi_input_atom_list_t*) malloc(sizeof(eddi_input_atom_list_t))))
    {
        EDDI_DEBUG_PRINT("[ERROR] Could not allocate atom object.");
        eddi_free_input_atom_list(*list);
        return EDDI_RETURN_FAILURE;
    }

    atom->x = x;
    atom->y = y;
    atom->z = z;
    atom->atomic_number = atomic_number;
    atom->next = NULL;

    if (*list == NULL)
    {
        *list = atom;
    }
    else 
    {
        atom->next = *list;
        *list = atom;
    }
    EDDI_DEBUG_PRINT("Added atom %lf %lf %lf %d\n", x, y, z, atomic_number);
    
    return EDDI_RETURN_SUCCESS;
}

void eddi_free_input_atom_list(eddi_input_atom_list_t* list)
{
    eddi_input_atom_list_t* it = list;
    while (list != NULL)
    {
        list = list->next;
        free(it);
        it = list;
    }
    free(list);
}
