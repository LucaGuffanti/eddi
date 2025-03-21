/**
 * @file eddi_density_field.h
 * @author Luca Guffanti
 * @brief Contains functions for the initilization, manipulation and destruction of the density field
 */

#ifndef __EDDI_DENSITY_FIELD_H__
#define __EDDI_DENSITY_FIELD_H__

#include "eddi_base_includes.h"

/**
 * @brief Type of a 3D density field. 
 */
typedef struct {
    /**
     * @brief 3d density field.
     */
    eddi_array_t field;
    
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

} eddi_density_field_t;


/**
 * @brief Constructs a new density field object given its fundamental parameters
 * 
 * @param density_field density field object
 * @param dx x extension of a voxel
 * @param dy y extension of a voxel
 * @param dz z extension of a voxel
 * @param origin origin point of the cubic domain
 * @param n_x number of points along the x direction
 * @param n_y number of points along the y direction
 * @param n_z number of points along the z direction
 * @return true if the density field is create correctly 
 * @return false if any error occurs during instantiation of the density field
 */
bool eddi_new_density_field(eddi_density_field_t* density_field, 
    const eddi_real_t dx,
    const eddi_real_t dy, 
    const eddi_real_t dz,
    const eddi_point_t* origin,
    const eddi_size_t n_x,
    const eddi_size_t n_y,
    const eddi_size_t n_z
);

#endif // __EDDI_DENSITY_FIELD_H__
