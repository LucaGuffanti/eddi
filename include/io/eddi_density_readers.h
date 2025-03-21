/**
 * @file eddi_density_readers.h
 * @author Luca Guffanti
 * @brief Contains function prototypes for reading density fields from different file formats
 */

#ifndef __EDDI_DENSITY_READERS_H__
#define __EDDI_DENSITY_READERS_H__

#include "eddi_density_field.h"
#include "eddi_molecule.h"

/**
 * @brief Reads a binary file containing a density field, initializing the density field
 * 
 * @param filename filename
 * @param density_field density field object
 * @return true if the file is read correctly
 * @return false if any error occurs while the file is being read
 */
bool eddi_read_binary(const char* filename, eddi_density_field_t* density_field);

/**
 * @brief Reads a gaussian cube file containing information about a density field and the molecule that generates it
 * 
 * @param filename filename
 * @param density_field density field object
 * @param molecule molecule object
 * @return true if the file is read correctly 
 * @return false if any error occurs while the file is being read
 */
bool eddi_read_gaussian_cube(const char* filename, eddi_density_field_t* density_field, eddi_molecule_t* molecule);


#endif // __EDDI_DENSITY_READERS_H__

