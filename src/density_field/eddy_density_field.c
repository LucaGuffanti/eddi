/**
 * @file eddy_density_field.c
 * @author Luca Guffanti
 * @brief Implementation of the density field 
 */

#include "eddi_density_field.h"



bool eddi_new_density_field(eddi_density_field_t* density_field, 
    const eddi_real_t dx,
    const eddi_real_t dy, 
    const eddi_real_t dz,
    const eddi_point_t* origin,
    const eddi_size_t n_x,
    const eddi_size_t n_y,
    const eddi_size_t n_z
)
{

    density_field->field = (eddi_array_t) malloc(sizeof(eddi_real_t) * n_x * n_y * n_z); 
    if (!density_field->field)
    {
        EDDI_DEBUG_PRINT("[ERROR] Density field object instatiation failed (internal field)\n");
        return EDDI_RETURN_FAILURE;
    }

    density_field->dx = dx;
    density_field->dy = dy;
    density_field->dz = dz;
    density_field->x_size = n_x;
    density_field->y_size = n_y;
    density_field->z_size = n_z;
    memmove(&(density_field->origin), origin, 3 * sizeof(eddi_real_t));

    return EDDI_RETURN_SUCCESS; // Allocation successful
}
