/**
 * @file eddi_atomic_data.h
 * @author Luca Guffanti
 * @brief Contains the atomic data used in the library, and functions to map numbers to symbols
 * 
 */

#ifndef __EDDI_ATOMIC_DATA_H__
#define __EDDI_ATOMIC_DATA_H__

#define MAX_ATOMIC_NUMBER 93

#include "uthash.h"
#include "eddi_base_includes.h"

/**
 * @brief Mapping between symbol and number to be added as an element of
 * a hashmap to ease data retrieval.
 * 
 */
typedef struct {
    /**
     * @brief Null terminated atomic symbol
     */
    char symbol[3];

    /**
     * @brief Associated atomic number
     */
    eddi_atomic_number_t number;

    /**
     * @brief Handle for the hashmap
     */
    UT_hash_handle hh;

} eddi_atomic_mapping_t;

extern eddi_atomic_mapping_t* symbol_to_number_map;

/**
 * @brief Initializes the hashmap
*/
void eddi_init_atomic_data();

/**
 * @brief Extracts the atomic number from the symbol, using a hashmap
 * @param symbol atomic symbol
 * @return atomic number associated to the symbol
 */
eddi_atomic_number_t eddi_symbol_to_number(const char* symbol);


#endif // __EDDI_ATOMIC_DATA_H__