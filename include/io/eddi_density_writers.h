/**
 * @file density_writers.h
 * @author Luca Guffanti
 * @brief Contains function prototypes for writing density fields to different file formats
 */

#ifndef __EDDI_DENSITY_WRITERS_H__
#define __EDDI_DENSITY_WRITERS_H__

#include "eddi_base_includes.h"

/**
 * @brief Writes the density field to a binary file with the following format
 * 
 * n_points_x n_points_y n_points_z
 * for z in 0, n_points_z-1
 *    for y in 0, n_points_y-1
 *      for x in 0, n_points_x-1
 *        density(z,y,x)
 * 
 * @param filename name of the output file
 * @param density_field_descriptor descriptor of the density field
 * @param molecule_descriptor descriptor of the molecule 
 * 
 * @returns true if the file is written successfully
 * @returns false if any problem occurrs during the writing operation.
 * 
 */
bool eddi_write_binary(const char* filename, const eddi_density_field_descriptor_t* density_field_descriptor, const eddi_molecule_descriptor_t* molecule_descriptor);

/**
 * @brief Writes the density field to a Gaussian cube file

 * @param filename name of the output file
 * @param density_field_descriptor descriptor of the density field
 * @param molecule_descriptor descriptor of the molecule
 * 
 * @return true if the file is written successfully
 * @return false if any problem occurs during the writing operation.
 */
bool eddi_write_gaussian(const char* filename, const eddi_density_field_descriptor_t* density_field_descriptor, const eddi_molecule_descriptor_t* molecule_descriptor);


#endif // __EDDI_DENSITY_WRITERS_H__