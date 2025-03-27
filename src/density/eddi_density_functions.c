/**
 * @file eddi_wavefunctions.h
 * @author Luca Guffanti
 * @brief Implements the wavefunctions based on Slater's type orbitals
 */

#include "eddi_density_functions.h"

eddi_real_t eddi_H_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    // 1s1
    return 0.3183098861837909 * exp(-2.0 * r);
}

eddi_real_t eddi_He_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    //1s2
    return 1.56385647082292 * exp(-2.0 * 1.700000 * r);
}

eddi_real_t eddi_Li_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Be_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_B_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_C_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_N_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_O_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 2 * 145.3189672691407 * exp(-2 * 7.7 * r) + 6 * 6.46600263938176 * pow(r, 2.0) * exp(-2 * 2.275 * r);
}

eddi_real_t eddi_F_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Ne_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Na_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Mg_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Al_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Si_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_P_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_S_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Cl_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t eddi_Ar_slater(const eddi_real_t r, const eddi_real_t theta, const eddi_real_t phi)
{
    return 0.0;
}

eddi_real_t (*eddi_densities[EDDI_NUM_DENSITIES])(eddi_real_t, eddi_real_t, eddi_real_t) = {
    eddi_H_slater,
    eddi_He_slater,
    eddi_Li_slater,
    eddi_Be_slater,
    eddi_B_slater,
    eddi_C_slater,
    eddi_N_slater,
    eddi_O_slater,
    eddi_F_slater,
    eddi_Ne_slater,
    eddi_Na_slater,
    eddi_Mg_slater,
    eddi_Al_slater,
    eddi_Si_slater,
    eddi_P_slater,
    eddi_S_slater,
    eddi_Cl_slater,
    eddi_Ar_slater
};