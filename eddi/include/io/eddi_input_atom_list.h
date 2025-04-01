/**
 * @file eddi_input_atom_list.h
 * @author Luca Guffanti
 * @brief Contains function prototypes for reading molecules from a list of atoms.
 * Consists in an intermediate step for reading molecules from different file formats.
 * 
 */

#ifndef __EDDI_INPUT_ATOM_LIST_H__
#define __EDDI_INPUT_ATOM_LIST_H__

#include "eddi_base_includes.h"
#include "eddi_atomic_data.h"

/**
 * @brief List representing atomic species in input files
 * 
 */
typedef struct eddi_input_atom_list
{
    eddi_real_t x;
    eddi_real_t y;
    eddi_real_t z;
    eddi_atomic_number_t atomic_number;

    struct eddi_input_atom_list* next;

} eddi_input_atom_list_t;


/**
 * @brief Adds an atom to the list. The atom is specified by its chemical symbol.
 * 
 * @param list input atom list
 * @param x x position of the atom
 * @param y y position of the atom
 * @param z z position of the atom
 * @param symbol chemical symbol of the atom (2-character null terminated string)
 * @return true if the addition occurs without any errors
 * @return false if any error occurs during the addition of the atom
 */
bool eddi_add_atom_to_list_symbol(eddi_input_atom_list_t** list, const eddi_real_t x, const eddi_real_t y, const eddi_real_t z, const char symbol[3]);

/**
 * @brief Adds an atom to the list. The atom is specified by the atomic number
 * 
 * @param list input atom list
 * @param x x position of the atom
 * @param y y position of the atom
 * @param z z position of the atom
 * @param atomic_number atomic number of the atom
 * @return true if the addition occurs without any errors
 * @return false if any error occurs during the addition of the atom
 */
bool eddi_add_atom_to_list_atomic_number(eddi_input_atom_list_t** list, const eddi_real_t x, const eddi_real_t y, const eddi_real_t z, const eddi_atomic_number_t atomic_number);

/**
 * @brief Frees the atom list memory
 * 
 * @param list input atom list
 */
void eddi_free_input_atom_list(eddi_input_atom_list_t* list);

#endif // __EDDI_INPUT_ATOM_LIST_H__