/**
 * @file eddy_density_field.c
 * @author Luca Guffanti
 * @brief Implementation of the density field 
 */

#include "eddi_density_field.h"

bool eddi_new_density_field(
    eddi_density_field_t* density_field, 
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
    memset(density_field->field, 0, sizeof(eddi_real_t) * n_x * n_y * n_z);

    density_field->dx = dx;
    density_field->dy = dy;
    density_field->dz = dz;
    density_field->x_size = n_x;
    density_field->y_size = n_y;
    density_field->z_size = n_z;
    memcpy(&(density_field->origin), origin, 3 * sizeof(eddi_real_t));

    return EDDI_RETURN_SUCCESS; // Allocation successful
}

bool eddi_init_field_from_molecule(
    eddi_density_field_t* density_field, 
    const eddi_molecule_t* molecule,
    const eddi_real_t padding,
    const eddi_real_t dx,
    const eddi_real_t dy,
    const eddi_real_t dz
)
{

    // Extract the mininum and maximum positions of the atoms along each direction
    eddi_real_t min_x = molecule->atoms_x[0];
    eddi_real_t max_x = molecule->atoms_x[0];

    eddi_real_t min_y = molecule->atoms_y[0];
    eddi_real_t max_y = molecule->atoms_y[0];

    eddi_real_t min_z = molecule->atoms_z[0];
    eddi_real_t max_z = molecule->atoms_z[0];

    for (eddi_size_t i = 1; i < molecule->n_atoms; ++i)
    {
        if (molecule->atoms_x[i] < min_x) min_x = molecule->atoms_x[i];
        if (molecule->atoms_x[i] > max_x) max_x = molecule->atoms_x[i];

        if (molecule->atoms_y[i] < min_y) min_y = molecule->atoms_y[i];
        if (molecule->atoms_y[i] > max_y) max_y = molecule->atoms_y[i];

        if (molecule->atoms_z[i] < min_z) min_z = molecule->atoms_z[i];
        if (molecule->atoms_z[i] > max_z) max_z = molecule->atoms_z[i];
    }

    // Add padding to the min and max values
    min_x -= padding;
    max_x += padding;

    min_y -= padding;
    max_y += padding;

    min_z -= padding;
    max_z += padding;

    // Calculate the number of grid points along each direction
    eddi_size_t n_x = (eddi_size_t)ceil((max_x - min_x) / dx);
    eddi_size_t n_y = (eddi_size_t)ceil((max_y - min_y) / dy);
    eddi_size_t n_z = (eddi_size_t)ceil((max_z - min_z) / dz);

    // Set the origin point
    eddi_point_t origin = {min_x, min_y, min_z};

    // Initialize the density field
    return eddi_new_density_field(density_field, dx, dy, dz, &origin, n_x, n_y, n_z);
}



void eddi_compute_density_field(eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    const eddi_real_t dx = density_field->dx;
    const eddi_real_t dy = density_field->dy;
    const eddi_real_t dz = density_field->dz;

    const eddi_real_t origin_x = density_field->origin.x; 
    const eddi_real_t origin_y = density_field->origin.y;
    const eddi_real_t origin_z = density_field->origin.z;

    const eddi_size_t nx = density_field->x_size;
    const eddi_size_t ny = density_field->y_size;
    const eddi_size_t nz = density_field->z_size;

    EDDI_DEBUG_PRINT("dx: %f, dy: %f, dz: %f\n", dx, dy, dz);
    EDDI_DEBUG_PRINT("origin_x: %f, origin_y: %f, origin_z: %f\n", origin_x, origin_y, origin_z);
    EDDI_DEBUG_PRINT("nx: %zu, ny: %zu, nz: %zu\n", nx, ny, nz);

    const eddi_size_t atoms = molecule->n_atoms;
    const eddi_real_t dx_2 = dx * dx;

    eddi_array_t x = molecule->atoms_x;
    eddi_array_t y = molecule->atoms_y;
    eddi_array_t z = molecule->atoms_z;

    eddi_array_t field = density_field->field;

    const eddi_real_t cutoff_radius = 10.0;
    const eddi_real_t max_radius_squared = cutoff_radius * cutoff_radius;
    const size_t n_radii = (size_t) (max_radius_squared);
    eddi_real_t radii_sqrts[n_radii];
    
    // printf("Maximum radius %lf\n", max_radius_squared);
    // printf("Computing radii\n");
    // #pragma omp parallel for
    for (int i = 0; i < n_radii; ++i)
    {
        radii_sqrts[i] = sqrt((eddi_real_t)i); 
    }



    const eddi_size_t n_pts = nx * ny * nz;
    eddi_real_t radius;
    #pragma omp parallel for private(radius)
    for (eddi_size_t idx = 0; idx < n_pts; ++idx)
    {
        const eddi_size_t it_z = idx % nz;
        const eddi_size_t it_y = (idx / nz) % ny;
        const eddi_size_t it_x = (idx / nz) / ny;

        const eddi_real_t cx = origin_x + it_x * dx;
        const eddi_real_t cy = origin_y + it_y * dy;
        const eddi_real_t cz = origin_z + it_z * dz;


        eddi_real_t local_density = 0.0;

        for (eddi_size_t atom_idx = 0; atom_idx < atoms; ++atom_idx)
        {
            
            const eddi_real_t DeltaX = cx - x[atom_idx];
            const eddi_real_t DeltaY = cy - y[atom_idx];
            const eddi_real_t DeltaZ = cz - z[atom_idx];
            const eddi_real_t radius_2 = DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ;
            
            if (radius_2 < dx_2 || radius_2 > max_radius_squared) 
            {
                continue;
            } 
            else
            {
                radius = radii_sqrts[(size_t)radius_2];
                // printf("Radius sq %lf goes to %zu %lf\n", radius_2, (size_t) radius_2, radius);
                local_density += molecule->density[atom_idx](radius, 0.0, 0.0);
            }

        }
        field[idx] = local_density;
    }
}

void eddi_compute_density_field_cl(eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    // Here we also construct the cell list because this is just experimental

    const eddi_real_t cutoff_radius =6.0;
    const eddi_real_t max_radius_squared = cutoff_radius * cutoff_radius;
    const size_t n_radii = (size_t) (max_radius_squared);
    eddi_real_t radii_sqrts[n_radii];
    
    // printf("Maximum radius %lf\n", max_radius_squared);
    // printf("Computing radii\n");
    // #pragma omp parallel for
    for (int i = 0; i < n_radii; ++i)
    {
        radii_sqrts[i] = sqrt((eddi_real_t)i); 
    }


    const eddi_real_t c_x = 10.0; 
    const eddi_real_t c_y = 10.0;
    const eddi_real_t c_z = 10.0;

    const eddi_size_t nc_x = density_field->x_size / c_x;
    const eddi_size_t nc_y = density_field->y_size / c_y;
    const eddi_size_t nc_z = density_field->z_size / c_z;
    const eddi_size_t nc_tot = nc_x * nc_y * nc_z;

    const eddi_real_t multiplier = molecule->n_atoms * 0.1;
    const eddi_size_t c_n_atoms = multiplier + 1;


    eddi_molecule_t cells[nc_tot];
    eddi_size_t occupancy[nc_tot];
    
    memset(occupancy, sizeof(eddi_size_t) * nc_tot, 0);
    for (eddi_size_t i = 0; i < nc_tot; ++i)
    {
        cells[i].atoms_x = (eddi_real_t*) malloc(sizeof(eddi_real_t) * c_n_atoms);
        cells[i].atoms_y = (eddi_real_t*) malloc(sizeof(eddi_real_t) * c_n_atoms);
        cells[i].atoms_z = (eddi_real_t*) malloc(sizeof(eddi_real_t) * c_n_atoms);
        cells[i].density = malloc(sizeof(ptrdiff_t) * c_n_atoms);
    }

    // Now copy the atom data in the cells
    for (eddi_size_t i = 0; i < molecule->n_atoms; ++i)
    {
        const eddi_size_t cell_x = (eddi_size_t)((molecule->atoms_x[i] - density_field->origin.x) / c_x);
        const eddi_size_t cell_y = (eddi_size_t)((molecule->atoms_y[i] - density_field->origin.y) / c_y);
        const eddi_size_t cell_z = (eddi_size_t)((molecule->atoms_z[i] - density_field->origin.z) / c_z);

        if (cell_x >= nc_x || cell_y >= nc_y || cell_z >= nc_z)
        {
            printf("[WARNING] Atom is outside the domain boundaries and will be ignored.\n");
            continue;
        }

        const eddi_size_t cell_idx = cell_x * nc_y * nc_z + cell_y * nc_z + cell_z;
        cells[cell_idx].atoms_x[occupancy[cell_idx]] = molecule->atoms_x[i];
        cells[cell_idx].atoms_y[occupancy[cell_idx]] = molecule->atoms_y[i];
        cells[cell_idx].atoms_z[occupancy[cell_idx]] = molecule->atoms_z[i];
        cells[cell_idx].density[occupancy[cell_idx]] = molecule->density[i];
        occupancy[cell_idx]++;
    } 

    printf("Cell list with %zu cells. %zu each\n", nc_tot, c_n_atoms);

    const eddi_size_t n_points = density_field->x_size * density_field->y_size * density_field->z_size;
    
    #pragma omp parallel for
    for (eddi_size_t p_idx = 0; p_idx < n_points; ++p_idx)
    {

        const eddi_size_t i_z = p_idx % density_field->z_size;
        const eddi_size_t i_y = (p_idx / density_field->z_size) % density_field->y_size;
        const eddi_size_t i_x = (p_idx / density_field->z_size) / density_field->y_size;

        const eddi_real_t pz = density_field->origin.z + i_z * density_field->dz;
        const eddi_real_t py = density_field->origin.y + i_y * density_field->dy;
        const eddi_real_t px = density_field->origin.x + i_x * density_field->dx;


        const eddi_size_t cell_x = (eddi_size_t)((px - density_field->origin.x) / c_x);
        const eddi_size_t cell_y = (eddi_size_t)((py - density_field->origin.y) / c_y);
        const eddi_size_t cell_z = (eddi_size_t)((pz - density_field->origin.z) / c_z);

        // printf("Point (%f, %f, %f) is in cell (%zu, %zu, %zu)\n", px, py, pz, cell_x, cell_y, cell_z);

        if (cell_x >= nc_x || cell_y >= nc_y || cell_z >= nc_z)
        {
            printf("[WARNING] Point is outside the domain boundaries and will be ignored.\n");
            continue;
        }

        const eddi_size_t cell_id = cell_x * nc_y * nc_z + cell_y * nc_z + cell_z;
        
        // Access the 27 cells that contain and neighbor the point
        eddi_real_t local = 0.0;
        for (int x = -1; x <= 1; ++x)
        {
            // printf("%zu\n", x);
            for (int y = -1; y <= 1; ++y)
            {
                // printf("%zu\n", y);
                for (int z = -1; z <= 1; ++z)
                {

                    // printf("%zu\n", z);
                    
                    const eddi_size_t neighbor_x = cell_x + x;
                    const eddi_size_t neighbor_y = cell_y + y;
                    const eddi_size_t neighbor_z = cell_z + z;
                    
                    // printf("Neighbor cell coordinates: (%zu, %zu, %zu)\n", neighbor_x, neighbor_y, neighbor_z);

                    if (neighbor_x >= nc_x || neighbor_y >= nc_y || neighbor_z >= nc_z)
                        {
                            continue; // Skip out-of-bound neighbors
                        }
                        
                    const eddi_size_t neighbor_id = neighbor_x * nc_y * nc_z + neighbor_y * nc_z + neighbor_z;
                    // printf("Neighbor cell ID: %zu\n", neighbor_id);

                    // Access the atoms in the neighboring cell
                    if (occupancy[neighbor_id] == 0)
                    {
                        continue; // Skip empty cells
                    }

                    for (eddi_size_t atom_idx = 0; atom_idx < occupancy[neighbor_id]; ++atom_idx)
                    {
                        // printf("Atom ID: %zu\n", atom_idx);
                        const eddi_real_t atom_x = cells[neighbor_id].atoms_x[atom_idx];
                        const eddi_real_t atom_y = cells[neighbor_id].atoms_y[atom_idx];
                        const eddi_real_t atom_z = cells[neighbor_id].atoms_z[atom_idx];

                        // printf("Atom coordinates: (%f, %f, %f)\n", atom_x, atom_y, atom_z);

                        // Perform computations with the atom data
                        const eddi_real_t dx = px - atom_x;
                        const eddi_real_t dy = py - atom_y;
                        const eddi_real_t dz = pz - atom_z;
                        const eddi_real_t distance_squared = dx * dx + dy * dy + dz * dz;
                        
                        if (distance_squared < max_radius_squared && distance_squared > 4.0)
                            local += cells[neighbor_id].density[atom_idx](radii_sqrts[(eddi_size_t) distance_squared], 0.0, 0.0); 
                    }
                }
            }
        }
        density_field->field[p_idx] = local;
    }
}

void eddi_compute_density_field_cl_opt(eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    // Here we also construct the cell list because this is just experimental

    const eddi_real_t cutoff_radius = 8.0;
    const eddi_real_t max_radius_squared = cutoff_radius * cutoff_radius;
    const size_t n_radii = (size_t) (max_radius_squared);
    // eddi_real_t radii_sqrts[n_radii];
    
    // // printf("Maximum radius %lf\n", max_radius_squared);
    // // printf("Computing radii\n");
    // // #pragma omp parallel for
    // for (int i = 0; i < n_radii; ++i)
    // {
    //     radii_sqrts[i] = sqrt((eddi_real_t)i); 
    // }

    EDDI_DEBUG_PRINT("Constructed radii\n");
    const eddi_real_t c_x = 10.0; 
    const eddi_real_t c_y = 10.0;
    const eddi_real_t c_z = 10.0;

    const eddi_size_t nc_x = density_field->x_size / c_x;
    const eddi_size_t nc_y = density_field->y_size / c_y;
    const eddi_size_t nc_z = density_field->z_size / c_z;
    const eddi_size_t nc_tot = nc_x * nc_y * nc_z;

    const eddi_real_t multiplier = molecule->n_atoms * 0.1;
    const eddi_size_t c_n_atoms = multiplier + 1;


    eddi_molecule_t cells[nc_tot];
    eddi_size_t occupancy[nc_tot];
    EDDI_DEBUG_PRINT("Allocated cells and occupancy\n");

    memset(occupancy, 0, sizeof(eddi_size_t) * nc_tot);
    for (eddi_size_t i = 0; i < nc_tot; ++i)
    {
        EDDI_DEBUG_PRINT("Initializing cell %d\n", i);
        cells[i].atoms_x = (eddi_real_t*) malloc(sizeof(eddi_real_t) * c_n_atoms);
        cells[i].atoms_y = (eddi_real_t*) malloc(sizeof(eddi_real_t) * c_n_atoms);
        cells[i].atoms_z = (eddi_real_t*) malloc(sizeof(eddi_real_t) * c_n_atoms);
        cells[i].density = malloc(sizeof(ptrdiff_t) * c_n_atoms);
    }

    EDDI_DEBUG_PRINT("Initialized cells\n");
    EDDI_DEBUG_PRINT("Adding atoms to cells\n");

    // Now copy the atom data in the cells
    for (eddi_size_t i = 0; i < molecule->n_atoms; ++i)
    {
        EDDI_DEBUG_PRINT("Atom %d\n", i);
        const eddi_size_t cell_x = (eddi_size_t)((molecule->atoms_x[i] - density_field->origin.x) / c_x);
        const eddi_size_t cell_y = (eddi_size_t)((molecule->atoms_y[i] - density_field->origin.y) / c_y);
        const eddi_size_t cell_z = (eddi_size_t)((molecule->atoms_z[i] - density_field->origin.z) / c_z);

        if (cell_x >= nc_x || cell_y >= nc_y || cell_z >= nc_z)
        {
            printf("[WARNING] Atom is outside the domain boundaries and will be ignored.\n");
            continue;
        }

        const eddi_size_t cell_idx = cell_x * nc_y * nc_z + cell_y * nc_z + cell_z;
        cells[cell_idx].atoms_x[occupancy[cell_idx]] = molecule->atoms_x[i];
        cells[cell_idx].atoms_y[occupancy[cell_idx]] = molecule->atoms_y[i];
        cells[cell_idx].atoms_z[occupancy[cell_idx]] = molecule->atoms_z[i];
        cells[cell_idx].density[occupancy[cell_idx]] = molecule->density[i];
        occupancy[cell_idx]++;
    } 

    printf("Cell list with %zu cells. %zu each\n", nc_tot, c_n_atoms);

    const eddi_size_t n_points = density_field->x_size * density_field->y_size * density_field->z_size;
    
    #pragma omp parallel for
    for (eddi_size_t p_idx = 0; p_idx < n_points; ++p_idx)
    {

        const eddi_size_t i_z = p_idx % density_field->z_size;
        const eddi_size_t i_y = (p_idx / density_field->z_size) % density_field->y_size;
        const eddi_size_t i_x = (p_idx / density_field->z_size) / density_field->y_size;

        const eddi_real_t pz = density_field->origin.z + i_z * density_field->dz;
        const eddi_real_t py = density_field->origin.y + i_y * density_field->dy;
        const eddi_real_t px = density_field->origin.x + i_x * density_field->dx;


        const eddi_size_t cell_i = (eddi_size_t)((px - density_field->origin.x) / c_x);
        const eddi_size_t cell_j = (eddi_size_t)((py - density_field->origin.y) / c_y);
        const eddi_size_t cell_k = (eddi_size_t)((pz - density_field->origin.z) / c_z);

        // printf("Point (%f, %f, %f) is in cell (%zu, %zu, %zu)\n", px, py, pz, cell_x, cell_y, cell_z);

        if (cell_i >= nc_x || cell_j >= nc_y || cell_k >= nc_z)
        {
            printf("[WARNING] Point is outside the domain boundaries and will be ignored.\n");
            continue;
        }
        const eddi_size_t cell_id = cell_i * nc_y * nc_z + cell_j * nc_z + cell_k;
        
        const eddi_real_t d_i_m1 = (px - density_field->origin.x - (cell_i) * c_x) * (px - density_field->origin.x -(cell_i) * c_x);
        const eddi_real_t d_i_p1 = (px - density_field->origin.x - (cell_i + 1) * c_x) * (px - density_field->origin.x -(cell_i + 1) * c_x); 
        const eddi_real_t d_j_m1 = (py - density_field->origin.y - (cell_j) * c_y) * (py - density_field->origin.y - (cell_j) * c_y); 
        const eddi_real_t d_j_p1 = (py - density_field->origin.y - (cell_j + 1) * c_y) * (py - density_field->origin.y - (cell_j + 1) * c_y); 
        const eddi_real_t d_k_m1 = (pz - density_field->origin.z - (cell_k) * c_z) * (pz - density_field->origin.z - (cell_k) * c_z); 
        const eddi_real_t d_k_p1 = (pz - density_field->origin.z - (cell_k + 1) * c_z) * (pz - density_field->origin.z - (cell_k + 1) * c_z); 

        const eddi_real_t d_m1_m1_m1 = d_i_m1 + d_j_m1 + d_k_m1;
        const eddi_real_t d_m1_m1_c0 = d_i_m1 + d_j_m1         ;
        const eddi_real_t d_m1_m1_p1 = d_i_m1 + d_j_m1 + d_k_p1;
        const eddi_real_t d_m1_c0_m1 = d_i_m1          + d_k_m1; 
        const eddi_real_t d_m1_c0_c0 = d_i_m1                  ;
        const eddi_real_t d_m1_c0_p1 = d_i_m1          + d_k_p1; 
        const eddi_real_t d_m1_p1_m1 = d_i_m1 + d_j_p1 + d_k_m1;
        const eddi_real_t d_m1_p1_c0 = d_i_m1 + d_j_p1         ;
        const eddi_real_t d_m1_p1_p1 = d_i_m1 + d_j_p1 + d_k_p1;


        const eddi_real_t d_c0_m1_m1 =          d_j_m1 + d_k_m1;
        const eddi_real_t d_c0_m1_c0 =          d_j_m1;
        const eddi_real_t d_c0_m1_p1 =          d_j_m1 + d_k_p1;
        const eddi_real_t d_c0_c0_m1 =                   d_k_m1;
        const eddi_real_t d_c0_c0_c0 =                      0.0;
        const eddi_real_t d_c0_c0_p1 =                   d_k_p1;
        const eddi_real_t d_c0_p1_m1 =          d_j_p1 + d_k_m1;
        const eddi_real_t d_c0_p1_c0 =          d_j_p1;
        const eddi_real_t d_c0_p1_p1 =          d_j_p1 + d_k_p1;


        const eddi_real_t d_p1_m1_m1 = d_i_p1 + d_j_m1 + d_k_m1;
        const eddi_real_t d_p1_m1_c0 = d_i_p1 + d_j_m1         ;
        const eddi_real_t d_p1_m1_p1 = d_i_p1 + d_j_m1 + d_k_p1;
        const eddi_real_t d_p1_c0_m1 = d_i_p1          + d_k_m1; 
        const eddi_real_t d_p1_c0_c0 = d_i_p1                  ;
        const eddi_real_t d_p1_c0_p1 = d_i_p1          + d_k_p1; 
        const eddi_real_t d_p1_p1_m1 = d_i_p1 + d_j_p1 + d_k_m1;
        const eddi_real_t d_p1_p1_c0 = d_i_p1 + d_j_p1         ;
        const eddi_real_t d_p1_p1_p1 = d_i_p1 + d_j_p1 + d_k_p1;

        const eddi_real_t arr[] = {
            d_m1_m1_m1, d_m1_m1_c0, d_m1_m1_p1,
            d_m1_c0_m1, d_m1_c0_c0, d_m1_c0_p1,
            d_m1_p1_m1, d_m1_p1_c0, d_m1_p1_p1,
            d_c0_m1_m1, d_c0_m1_c0, d_c0_m1_p1,
            d_c0_c0_m1, d_c0_c0_c0, d_c0_c0_p1,
            d_c0_p1_m1, d_c0_p1_c0, d_c0_p1_p1,
            d_p1_m1_m1, d_p1_m1_c0, d_p1_m1_p1,
            d_p1_c0_m1, d_p1_c0_c0, d_p1_c0_p1,
            d_p1_p1_m1, d_p1_p1_c0, d_p1_p1_p1
        };

        eddi_real_t partial = 0.0;
        for (int inner_i = 0; inner_i <= 2; ++inner_i)
        {
            for (int inner_j = 0; inner_j <= 2; ++inner_j)
            {
                for (int inner_k = 0; inner_k <= 2; ++inner_k)
                {
                    // Move the computation of the distance here.

                    const size_t idx = inner_i * 3 * 3 + inner_j * 3 + inner_k;
                    if (arr[idx] < max_radius_squared)
                    {
                        const size_t neighbor_i = cell_i + inner_i - 1;
                        const size_t neighbor_j = cell_j + inner_j - 1;
                        const size_t neighbor_k = cell_k + inner_k - 1;


                        if (neighbor_i < nc_x && neighbor_j < nc_y && neighbor_k < nc_z)
                        {
                            const size_t neighbor_index = neighbor_i * nc_y * nc_z + neighbor_j * nc_z + neighbor_k;
                            if (occupancy[neighbor_index] > 0) {
                                for (size_t idx = 0; idx < occupancy[neighbor_index]; ++idx) {
                                    const eddi_real_t u_x = cells[neighbor_index].atoms_x[idx];
                                    const eddi_real_t u_y = cells[neighbor_index].atoms_y[idx];
                                    const eddi_real_t u_z = cells[neighbor_index].atoms_z[idx];

                                    const eddi_real_t delta_x = px - u_x;
                                    const eddi_real_t delta_y = py - u_y;
                                    const eddi_real_t delta_z = pz - u_z;

                                    const eddi_real_t distance_squared = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

                                    if (distance_squared < max_radius_squared && distance_squared > 4.0) {
                                        partial += cells[neighbor_index].density[idx](sqrt(distance_squared), 0.0, 0.0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        density_field->field[p_idx] = partial;
    }
}


void eddi_compute_density_field_atom(eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    const eddi_real_t cutoff_radius = 8.0;
    const eddi_real_t max_radius_squared = cutoff_radius * cutoff_radius;
    const size_t n_radii = (size_t) (max_radius_squared);
    eddi_real_t radii_sqrts[n_radii];
    
    // printf("Maximum radius %lf\n", max_radius_squared);
    // printf("Computing radii\n");
    
    const eddi_size_t atoms = molecule->n_atoms;
    const eddi_size_t n_x = density_field->x_size;
    const eddi_size_t n_y = density_field->y_size;
    const eddi_size_t n_z = density_field->z_size;

    const eddi_real_t dx = density_field->dx;
    const eddi_real_t dy = density_field->dy;
    const eddi_real_t dz = density_field->dz;

    const eddi_real_t ox = density_field->origin.x;
    const eddi_real_t oy = density_field->origin.y;
    const eddi_real_t oz = density_field->origin.z;

    #pragma omp parallel for shared(density_field)
    for (eddi_size_t atom_id = 0; atom_id < atoms; ++atom_id)
    {
        const eddi_real_t starting_x = molecule->atoms_x[atom_id] - cutoff_radius;
        const eddi_real_t starting_y = molecule->atoms_y[atom_id] - cutoff_radius;
        const eddi_real_t starting_z = molecule->atoms_z[atom_id] - cutoff_radius;

        const eddi_real_t end_x = molecule->atoms_x[atom_id] + cutoff_radius;
        const eddi_real_t end_y = molecule->atoms_y[atom_id] + cutoff_radius;
        const eddi_real_t end_z = molecule->atoms_z[atom_id] + cutoff_radius;

        for (eddi_real_t x = starting_x;  x <= end_x; x += dx)
        {
            for (eddi_real_t y = starting_y; y <= end_y; y += dy)
            {
                for (eddi_real_t z = starting_z; z <= end_z; z += dz)
                {
                    const eddi_size_t idx_x = (eddi_size_t)((x - ox) / dx);
                    const eddi_size_t idx_y = (eddi_size_t)((y - oy) / dy);
                    const eddi_size_t idx_z = (eddi_size_t)((z - oz) / dz);

                    if (idx_x < n_x && idx_y < n_y && idx_z < n_z)
                    {
                        const eddi_size_t idx = idx_x * n_y * n_z + idx_y * n_z + idx_z;
                        #pragma omp atomic
                        density_field->field[idx] += molecule->density[atom_id](sqrt((x - molecule->atoms_x[atom_id]) * (x - molecule->atoms_x[atom_id]) +
                                                                                    (y - molecule->atoms_y[atom_id]) * (y - molecule->atoms_y[atom_id]) +
                                                                                    (z - molecule->atoms_z[atom_id]) * (z - molecule->atoms_z[atom_id])), 0.0, 0.0);
                    }
                }
            }
        }
        
    }
}

void eddi_free_density_field(eddi_density_field_t* density_field)
{
    // Only need to deallocate the field
    free(density_field->field);
}