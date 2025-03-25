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

#define EDDI_NUM_DENSITIES 18


eddi_real_t eddi_H_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_He_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Li_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Be_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_B_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_C_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_N_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_O_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_F_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Ne_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Na_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Mg_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Al_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Si_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_P_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_S_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Cl_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
eddi_real_t eddi_Ar_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi);
// TODO: Add the remaining density functions

extern eddi_real_t (*eddi_densities[EDDI_NUM_DENSITIES])(eddi_real_t, eddi_real_t, eddi_real_t);


#endif // __EDDI_DENSITY_FUNCTIONS__