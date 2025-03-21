/**
 * @file eddi_moleucle.h
 * @author Luca Guffanti
 * @brief Contains functions to initialize, manipulate and destroy molecules
 */

#ifndef __EDDI_MOLECULE_H__
#define __EDDI_MOLECULE_H__

#include "eddi_base_includes.h"

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
    eddi_real_t (**density)(double, double, double);
    
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
bool eddi_destroy_molecule(eddi_molecule_t* molecule);


#endif // __EDDI_MOLECULE_H__
