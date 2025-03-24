/**
 * @file eddi_molecule.c
 * @author Luca Guffanti
 * @brief Implementation of the molecule functionalities
 */

#include "eddi_molecule.h"


bool eddi_new_molecule(eddi_molecule_t* molecule, eddi_size_t n_atoms, eddi_array_t x, eddi_array_t y, eddi_array_t z, eddi_atomic_number_t* atomic_numbers)
{
    molecule->n_atoms = n_atoms;
    
    molecule->atoms_x = (eddi_array_t) malloc(sizeof(eddi_real_t) * n_atoms);
    molecule->atoms_y = (eddi_array_t) malloc(sizeof(eddi_real_t) * n_atoms);
    molecule->atoms_z = (eddi_array_t) malloc(sizeof(eddi_real_t) * n_atoms);
    molecule->atomic_numbers = (eddi_atomic_number_t*) malloc(sizeof(eddi_atomic_number_t) * n_atoms);

    memcpy(molecule->atoms_x, x, sizeof(eddi_real_t) * n_atoms);
    memcpy(molecule->atoms_y, y, sizeof(eddi_real_t) * n_atoms);
    memcpy(molecule->atoms_z, z, sizeof(eddi_real_t) * n_atoms);
    memcpy(molecule->atomic_numbers, atomic_numbers, sizeof(char) * n_atoms);


    // TODO: associate a set of wavefunctions to each atom
    
    return EDDI_RETURN_SUCCESS;

}

bool eddi_atom_list_to_molecule(eddi_input_atom_list_t* list, eddi_molecule_t* molecule, eddi_size_t n_atoms)
{

    eddi_real_t x[n_atoms];
    eddi_real_t y[n_atoms];
    eddi_real_t z[n_atoms];
    eddi_atomic_number_t atomic_numbers[n_atoms];

    eddi_size_t it = 0;
    while (list != NULL)
    {
        x[it] = list->x;
        y[it] = list->y;
        z[it] = list->z;
        atomic_numbers[it] = list->atomic_number;

        list = list->next;
        ++it;
    }

    return eddi_new_molecule(molecule, n_atoms, x, y, z, atomic_numbers);
}

void eddi_print_molecule(eddi_molecule_t* molecule)
{
    for (eddi_size_t it = 0; it < molecule->n_atoms; ++it)
        printf("Atom %zu: (%lf, %lf, %lf), Atomic Number: %d\n", 
               it, 
               molecule->atoms_x[it], 
               molecule->atoms_y[it], 
               molecule->atoms_z[it], 
               molecule->atomic_numbers[it]);
}

