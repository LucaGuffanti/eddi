/**
 * @file eddi_atomic_data.c
 * @author Luca Guffanti
 * @brief Implements the functions defined in eddi_atomic_data.h 
 * 
 */

#include "eddi_atomic_data.h"

eddi_symbol_number_mapping_t* symbol_to_number_map = NULL;
eddi_symbol_number_mapping_t* mappings_space = NULL;

void eddi_init_symbol_mapping()
{
    eddi_symbol_number_mapping_t* mappings_space = (eddi_symbol_number_mapping_t*) malloc(sizeof(eddi_symbol_number_mapping_t) * MAX_ATOMIC_NUMBER);
    char* names[] = {
        "H\0", "HE\0", "LI\0", "BE\0", "B\0", "C\0", "N\0", "O\0", "F\0", "NE\0",
        "NA\0", "MG\0", "AL\0", "SI\0", "P\0", "S\0", "CL\0", "AR\0", "K\0", "CA\0",
        "SC\0", "TI\0", "V\0", "CR\0", "MN\0", "FE\0", "CO\0", "NI\0", "CU\0", "ZN\0",
        "GA\0", "GE\0", "AS\0", "SE\0", "BR\0", "KR\0", "RB\0", "SR\0", "Y\0", "ZR\0",
        "NB\0", "MO\0", "TC\0", "RU\0", "RH\0", "PD\0", "AG\0", "CD\0", "IN\0", "SN\0",
        "SB\0", "TE\0", "I\0", "XE\0", "CS\0", "BA\0", "LA\0", "CE\0", "PR\0", "ND\0",
        "PM\0", "SM\0", "EU\0", "GD\0", "TB\0", "DY\0", "HO\0", "ER\0", "TM\0", "YB\0",
        "LU\0", "HF\0", "TA\0", "W\0", "RE\0", "OS\0", "IR\0", "PT\0", "AU\0", "HG\0",
        "TL\0", "PB\0", "BI\0", "PO\0", "AT\0", "RN\0", "FR\0", "RA\0", "AC\0", "TH\0",
        "PA\0", "U\0", "NP\0"
    };

    for (size_t i = 0; i < MAX_ATOMIC_NUMBER; ++i)
    {
        eddi_symbol_number_mapping_t* current_mapping = mappings_space + i;
        current_mapping->number = i+1;
        memcpy(current_mapping->symbol, names[i], 3);
    
        HASH_ADD_STR(symbol_to_number_map, symbol, current_mapping);
    }
}

void eddi_init_atomic_data()
{

    // Initialize the symbol mapping
    EDDI_DEBUG_PRINT("Initializing symbol to number mapping\n");
    eddi_init_symbol_mapping();
}

eddi_atomic_number_t eddi_symbol_to_number(const char* symbol)
{
    if (symbol_to_number_map == NULL)
    {
        eddi_init_atomic_data();
        EDDI_DEBUG_PRINT("Initializing data\n");
    }

    eddi_symbol_number_mapping_t* mapping;
    EDDI_DEBUG_PRINT("Searching atomic number for %s\n", symbol);
    HASH_FIND_STR( symbol_to_number_map, symbol, mapping);
    EDDI_DEBUG_PRINT("Retrieved %d for %s\n", mapping->number, symbol);
    return mapping->number;
}

void eddi_free_atomic_data()
{
    free(mappings_space);
}