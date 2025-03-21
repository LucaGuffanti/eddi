/**
 * @file density_writers.h
 * @author Luca Guffanti
 * @brief Contains function prototypes for writing density fields to different file formats
 */

#ifndef __EDDI_DENSITY_WRITERS_H__
#define __EDDI_DENSITY_WRITERS_H__

#include "eddi_density_field.h"
#include "eddi_molecule.h"

/**
 * @brief Writes the density field to a binary file with the following format
 * dx dy dz
 * origin_x origin_y origin_z
 * n_points_x n_points_y n_points_z
 * for x in 0, n_points_x-1
 *    for y in 0, n_points_y-1
 *      for z in 0, n_points_z-1
 *        density(x,y,z)
 * 
 * @param filename name of the output file
 * @param density_field density field
 * 
 * @returns true if the file is written successfully
 * @returns false if any problem occurrs during the writing operation.
 * 
 */
bool eddi_write_binary(const char* filename, const eddi_density_field_t* density_field);

/**
 * @brief Writes the density field to a Gaussian cube file

 * @param filename name of the output file
 * @param density_field density field
 * @param molecule molecule
 * 
 * @return true if the file is written successfully
 * @return false if any problem occurs during the writing operation.
 */
bool eddi_write_gaussian_cube(const char* filename, const eddi_density_field_t* density_field, const eddi_molecule_t* molecule);


#endif // __EDDI_DENSITY_WRITERS_H__