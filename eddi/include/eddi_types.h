/**
 * @file eddi_types.h
 * @author Luca Guffanti
 * @brief Contains the definition of the types used in the eddy library
 * 
 */

#ifndef __EDDI_TYPES_H__
#define __EDDI_TYPES_H__

#include "stddef.h" 
#include "stdint.h"
#include "stdbool.h"
#include "complex.h"


/**
 * Definition of the type of floating point number
 */

#define EDDI_HIGH_PRECISION

#ifdef EDDI_HIGH_PRECISION 
#define EDDI_FLOATING_POINT_TYPE double
#else
#define EDDI_FLOATING_POINT_TYPE float
#endif

/**
 * @brief Type of the real numbers.
 */
typedef EDDI_FLOATING_POINT_TYPE eddi_real_t;

/**
 * @brief Type of the complex numbers.
 */
typedef EDDI_FLOATING_POINT_TYPE complex eddi_complex_t;

/**
 * @brief Type of the size of containers.
 */
typedef size_t eddi_size_t;

/**
 * @brief Array of real numbers.
 */
typedef eddi_real_t* eddi_array_t;

/**
 * @brief Type of atomic number. 
 */
typedef char eddi_atomic_number_t;

/**
 * @brief Type of point in 3D space. 
 */
typedef struct {
    eddi_real_t x;
    eddi_real_t y;
    eddi_real_t z;
} eddi_point_t;

typedef struct {
    eddi_real_t cx;
    eddi_real_t cy;
    eddi_real_t cz;
} eddi_cl_info_t;

typedef struct {
    eddi_real_t high_cutoff_radius;
    eddi_real_t low_cutoff_radius;
} eddi_cutoff_data_t;

#endif // __EDDI_TYPES_H__