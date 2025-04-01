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
#define EDDI_FLOATING_POINT_TYPE double

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



#endif // __EDDI_TYPES_H__