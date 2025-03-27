/**
 * @file eddi_wavefunctions.h
 * @author Luca Guffanti
 * @brief Contains Density functions. Atoms are mapped to the appropriate wavefunction based
 * on the atomic number. 
 */

#ifndef __EDDI_DENSITY_FUNCTIONS__
#define __EDDI_DENSITY_FUNCTIONS__

#include "eddi_base_includes.h"
#include "math.h"

#define EDDI_N_SLATER_DENSITIES 26

extern eddi_real_t (*eddi_densities[EDDI_N_SLATER_DENSITIES])(eddi_real_t, eddi_real_t, eddi_real_t);
    
#endif // __EDDI_DENSITY_FUNCTIONS__