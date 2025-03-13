/**
 * @file density_readers.h
 * @author Luca Guffanti
 * @brief Contains function prototypes for reading density fields from different file formats
 */

#ifndef __EDDI_DENSITY_READERS_H__
#define __EDDI_DENSITY_READERS_H__

#include "eddi_base_includes.h"

bool eddi_read_binary(const char* filename, eddi_density_field_descriptor_t* density_field_descriptor, eddi_molecule_descriptor_t* molecule_descriptor);
bool eddi_read_gaussian_cube(const char* filename, eddi_density_field_descriptor_t* density_field_descriptor, eddi_molecule_descriptor_t* molecule_descriptor);



#endif // __EDDI_DENSITY_READERS_H__

