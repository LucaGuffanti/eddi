/**
 * @file eddi_macro.h
 * @author Luca Guffanti
 * @brief Contains a number of preprocessor macros
 */

#ifndef __EDDI_MACRO_H__
#define __EDDI_MACRO_H__

// Enable and disable prints throughout the library
#ifdef EDDI_DEBUG
    #define EDDI_DEBUG_PRINT(...) printf(__VA_ARGS__)
    #define EDDI_DEBUG_CALL(...) __VA_ARGS__
#else
    #define EDDI_DEBUG_PRINT(...)
    #define EDDI_DEBUG_CALL(...)
#endif

#endif // __EDDI_MACRO_H__