/**
 * @file eddi_moleucle.h
 * @author Luca Guffanti
 * @brief Contains functions to initialize, manipulate and destroy molecules
 */

#ifndef __EDDI_MOLECULE_H__
#define __EDDI_MOLECULE_H__

#include "eddi_base_includes.h"
#include "eddi_input_atom_list.h"
#include "eddi_density_functions.h"

/**
 * @brief of a molecule (Structure of Arrays).
 */
typedef struct {

    /**
     * @brief x coordinates of the atoms.
     */
    eddi_real_t* atoms_x;

    /**
     * @brief y coordinates of the atoms.
     */
    eddi_real_t* atoms_y;

    /**
     * @brief z coordinates of the atoms.
     */
    eddi_real_t* atoms_z;

    /**
     * @brief density functions of the atoms.
     */
    eddi_real_t (**density)(eddi_real_t, eddi_real_t, eddi_real_t);
    
    /**
     * @brief List of atomic numbers;
     */
    eddi_atomic_number_t* atomic_numbers;

    /**
     * @brief Number of atoms in the molecule.
     */
    eddi_size_t n_atoms;

} eddi_molecule_t;

/**
 * @brief Allocates a molecule object by instantiating all the components
 * of the
 * @param n_atoms number of atoms in the molecule
 * @return eddi_molecule_t* pointer to the allocated molecule
 */
bool eddi_new_molecule(eddi_molecule_t* molecule, eddi_size_t n_atoms, eddi_array_t x, eddi_array_t y, eddi_array_t z, eddi_atomic_number_t* atomic_numbers);

/**
 * @brief Translates an atom list to a molecule
 * @param list input atom list
 * @param molecule molecule object
 */
bool eddi_atom_list_to_molecule(eddi_input_atom_list_t* list, eddi_molecule_t* molecule, eddi_size_t n_atoms);


/**
 * @brief Prints the atoms in a molecule
 */
void eddi_print_molecule(eddi_molecule_t* molecule);

/**
 * @brief Frees the memory allocated to a molecule
 */
void eddi_free_molecule(eddi_molecule_t* molecule);

#endif // __EDDI_MOLECULE_H__
