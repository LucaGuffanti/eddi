/**
 * @file eddi_types.h
 * @author Luca Guffanti
 * @brief Contains the definition of the types used in the eddy library
 * 
 */

#ifndef __EDDI_TYPES_H__
#define __EDDI_TYPES_H__

#include <stddef.h> 
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

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
 * @brief Type of the density field (3D array of real numbers, linearized for memory efficiency).
 */
typedef eddi_real_t* eddi_density_field_t;

/**
 * @brief Type of atomic number. 
 */
typedef uint8_t eddi_atomic_number_t;

/**
 * @brief Type of point in 3D space. 
 */
typedef struct {
    eddi_real_t x;
    eddi_real_t y;
    eddi_real_t z;
} eddi_point_t;

/**
 * @brief Type of the descriptor of a 3D density field. 
 */
typedef struct {
    /**
     * @brief 3d density field.
     */
    eddi_density_field_t density_field;
    
    /**
     * @brief Size of the density field in the x direction in terms of number of points.
     */
    eddi_size_t x_size;

    /**
     * @brief Size of the density field in the y direction in terms of number of points.
     */
    eddi_size_t y_size;

    /**
     * @brief Size of the density field in the z direction in terms of number of points.
     */
    eddi_size_t z_size;
    
    /**
     * @brief Resolution along the x axis in units of Bohrs.
     */
    eddi_real_t dx;

    /**
     * @brief Resolution along the y axis in units of Bohrs.
     */
    eddi_real_t dy;

    /**
     * @brief Resolution along the z axis in units of Bohrs.
     */
    eddi_real_t dz;

    /**
     * @brief Coordinates of the origin of the density field in 3D space in units of Bohrs.
     * 
     */
    eddi_point_t origin;

} eddi_density_field_descriptor_t;

/**
 * @brief Descriptor of an atom. 
 */
typedef struct {
    /**
     * @brief Position of the atom in 3D in units of Bohrs. 
     */
    eddi_point_t position;

    /**
     * @brief Atomic number. 
     */
    eddi_atomic_number_t atomic_number;

    /**
     * @brief Position Wavefuntion, as a function pointer, taking a position expressed in 
     * spherical coordinates as an input
     */
    eddi_real_t (* position_wavefunction)(double, double, double);
    
    
    /**
     * @brief Electron Density function, as a function pointer 
     * taking a position expressed in spherical coordinates as an input
     */
    eddi_real_t (* electron_density_function)(double, double, double);

} eddi_atom_descriptor_t;

/**
 * @brief Descriptor of an atom.
 */
typedef struct {
    /**
     * @brief List of atom descriptors.
     */
    eddi_atom_descriptor_t* atoms;

    /**
     * @brief Number of atoms in the molecule.
     */
    eddi_size_t n_atoms;

} eddi_molecule_descriptor_t;

#endif // __EDDI_TYPES_H__