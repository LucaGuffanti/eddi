/**
 * @file eddi_moleucle.h
 * @author Luca Guffanti
 * @brief Contains functions to initialize, manipulate and destroy molecule descriptors
 */

#ifndef __EDDI_MOLECULE_H__
#define __EDDI_MOLECULE_H__

#include "eddi_base_includes.h"

void eddi_init_molecule(eddi_molecule_descriptor_t* molecule_descriptor);
void eddi_destroy_molecule(eddi_molecule_descriptor_t* molecule_descriptor);


#endif // __EDDI_MOLECULE_H__
