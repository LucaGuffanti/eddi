/**
 * @file eddi_base_includes.h
 * @author Luca Guffanti
 * @brief Contains the basic includes for the EDDI library
 */

#ifndef __EDDI_BASE_INCLUDES_H__
#define __EDDI_BASE_INCLUDES_H__

#include "eddi_constants.h"
#include "eddi_macro.h"
#include "eddi_types.h"

#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "string.h"

#ifdef _OPENMP
    #include "omp.h"
#endif 

#endif // __EDDI_BASE_INCLUDES_H__