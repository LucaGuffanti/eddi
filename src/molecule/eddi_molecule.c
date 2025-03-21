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

    memmove(molecule->atoms_x, x, sizeof(eddi_real_t) * n_atoms);
    memmove(molecule->atoms_y, y, sizeof(eddi_real_t) * n_atoms);
    memmove(molecule->atoms_z, z, sizeof(eddi_real_t) * n_atoms);
    memmove(molecule->atomic_numbers, atomic_numbers, sizeof(char) * n_atoms);


    // TODO: associate a set of wavefunctions to each atom
    
    return EDDI_RETURN_SUCCESS;

}
