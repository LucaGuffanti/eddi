/**
 * @file eddy_density_field.c
 * @author Luca Guffanti
 * @brief Implementation of the density field 
 */

#include "eddi_density_field.h"


eddi_real_t slater(eddi_atomic_number_t z, eddi_real_t r, eddi_real_t theta, eddi_real_t phi)
{
    switch (z)
    {
    case 1:
    {
        // Slater wavefunction for atomic number 1
        return 0.3183098861837909 * exp(-2 * 1.0 * r);
    }

    case 2:
    {
        // Slater wavefunction for atomic number 2
        return 2 * 1.56385647082292 * exp(-2 * 1.7 * r);
    }

    case 3:
    {
        // Slater wavefunction for atomic number 3
        return 2 * 6.26529348975552 * exp(-2 * 2.7 * r) + 0.01231106589279565 * pow(r, 2.0) * exp(-2 * 0.65 * r);
    }

    case 4:
    {
        // Slater wavefunction for atomic number 4
        return 2 * 16.12335066486755 * exp(-2 * 3.7 * r) + 2 * 0.09348715662341694 * pow(r, 2.0) * exp(-2 * 0.9750000000000001 * r);
    }

    case 5:
    {
        // Slater wavefunction for atomic number 5
        return 2 * 33.04788731336516 * exp(-2 * 4.7 * r) + 3 * 0.39395410856946067 * pow(r, 2.0) * exp(-2 * 1.3 * r);
    }

    case 6:
    {
        // Slater wavefunction for atomic number 6
        return 2 * 58.94876275201326 * exp(-2 * 5.7 * r) + 4 * 1.2022525285933259 * pow(r, 2.0) * exp(-2 * 1.625 * r);
    }

    case 7:
    {
        // Slater wavefunction for atomic number 7
        return 2 * 95.73583629830304 * exp(-2 * 6.7 * r) + 5 * 2.9915890119493427 * pow(r, 2.0) * exp(-2 * 1.9500000000000002 * r);
    }

    case 8:
    {
        // Slater wavefunction for atomic number 8
        return 2 * 145.3189672691407 * exp(-2 * 7.7 * r) + 6 * 6.46600263938176 * pow(r, 2.0) * exp(-2 * 2.275 * r);
    }

    case 9:
    {
        // Slater wavefunction for atomic number 9
        return 2 * 209.6080149816861 * exp(-2 * 8.7 * r) + 7 * 12.606531474226806 * pow(r, 2.0) * exp(-2 * 2.6 * r);
    }

    case 10:
    {
        // Slater wavefunction for atomic number 10
        return 2 * 290.51283875301874 * exp(-2 * 9.7 * r) + 8 * 22.71737905948798 * pow(r, 2.0) * exp(-2 * 2.9250000000000003 * r);
    }

    case 11:
    {
        // Slater wavefunction for atomic number 11
        return 2 * 389.94329790024904 * exp(-2 * 10.7 * r) + 8 * 50.007109436657835 * pow(r, 2.0) * exp(-2 * 3.4250000000000003 * r) + 0.0016135320335943066 * pow(r, 4.0) * exp(-2 * 0.7333333333333331 * r);
    }

    case 12:
    {
        // Slater wavefunction for atomic number 12
        return 2 * 509.80925174048366 * exp(-2 * 11.7 * r) + 8 * 98.8387349088708 * pow(r, 2.0) * exp(-2 * 3.9250000000000003 * r) + 2 * 0.00987945178833327 * pow(r, 4.0) * exp(-2 * 0.9500000000000005 * r);
    }

    case 13:
    {
        // Slater wavefunction for atomic number 13
        return 2 * 652.0205595908079 * exp(-2 * 12.7 * r) + 8 * 180.00942628732454 * pow(r, 2.0) * exp(-2 * 4.425000000000001 * r) + 3 * 0.041619334990451416 * pow(r, 4.0) * exp(-2 * 1.1666666666666667 * r);
    }

    case 14:
    {
        // Slater wavefunction for atomic number 14
        return 2 * 818.487080768349 * exp(-2 * 13.7 * r) + 8 * 307.4397701413907 * pow(r, 2.0) * exp(-2 * 4.925000000000001 * r) + 4 * 0.13713727109030824 * pow(r, 4.0) * exp(-2 * 1.3833333333333335 * r);
    }

    case 15:
    {
        // Slater wavefunction for atomic number 15
        return 2 * 1011.1186745901936 * exp(-2 * 14.7 * r) + 8 * 498.57165616526265 * pow(r, 2.0) * exp(-2 * 5.425000000000001 * r) + 5 * 0.37975848643135146 * pow(r, 4.0) * exp(-2 * 1.6000000000000003 * r);
    }

    case 16:
    {
        // Slater wavefunction for atomic number 16
        return 2 * 1231.825200373448 * exp(-2 * 15.7 * r) + 8 * 774.7661644779737 * pow(r, 2.0) * exp(-2 * 5.925000000000001 * r) + 6 * 0.9238348519132616 * pow(r, 4.0) * exp(-2 * 1.8166666666666664 * r);
    }

    case 17:
    {
        // Slater wavefunction for atomic number 17
        return 2 * 1482.5165174352162 * exp(-2 * 16.7 * r) + 8 * 1161.7014531331229 * pow(r, 2.0) * exp(-2 * 6.425000000000001 * r) + 7 * 2.032954558631853 * pow(r, 4.0) * exp(-2 * 2.0333333333333337 * r);
    }

    case 18:
    {
        // Slater wavefunction for atomic number 18
        return 2 * 1765.1024850926005 * exp(-2 * 17.7 * r) + 8 * 1689.7706452760394 * pow(r, 2.0) * exp(-2 * 6.925000000000001 * r) + 8 * 4.129954204672846 * pow(r, 4.0) * exp(-2 * 2.25 * r);
    }

    case 19:
    {
        // Slater wavefunction for atomic number 19
        return 2 * 2081.4929626627027 * exp(-2 * 18.7 * r) + 8 * 2394.4797166230687 * pow(r, 2.0) * exp(-2 * 7.425000000000001 * r) + 8 * 10.862523867423388 * pow(r, 4.0) * exp(-2 * 2.5833333333333335 * r) + 2.9906737611881975e-05 * pow(r, 5.4) * exp(-2 * 0.5945945945945944 * r);
    }

    case 20:
    {
        // Slater wavefunction for atomic number 20
        return 2 * 2433.5978094626257 * exp(-2 * 19.7 * r) + 8 * 3316.845382815414 * pow(r, 2.0) * exp(-2 * 7.925000000000001 * r) + 8 * 25.402426141625856 * pow(r, 4.0) * exp(-2 * 2.9166666666666665 * r) + 2 * 0.00026309666655790993 * pow(r, 5.4) * exp(-2 * 0.7702702702702706 * r);
    }

    case 21:
    {
        // Slater wavefunction for atomic number 21
        return 2 * 2823.3268848094744 * exp(-2 * 20.7 * r) + 8 * 4503.792986750791 * pow(r, 2.0) * exp(-2 * 8.425 * r) + 8 * 54.18151395524859 * pow(r, 4.0) * exp(-2 * 3.25 * r) + 2 * 0.00040479679115536026 * pow(r, 5.4) * exp(-2 * 0.8108108108108107 * r) + 0.01414710605261292 * pow(r, 4.0) * exp(-2 * 1.0 * r);
    }

    case 22:
    {
        // Slater wavefunction for atomic number 22
        return 2 * 3252.590048020349 * exp(-2 * 21.7 * r) + 8 * 6008.554385940561 * pow(r, 2.0) * exp(-2 * 8.925 * r) + 8 * 107.31936045217444 * pow(r, 4.0) * exp(-2 * 3.5833333333333335 * r) + 2 * 0.0006098558208557976 * pow(r, 5.4) * exp(-2 * 0.8513513513513509 * r) + 2 * 0.05583016062783766 * pow(r, 4.0) * exp(-2 * 1.2166666666666661 * r);
    }

    case 23:
    {
        // Slater wavefunction for atomic number 23
        return 2 * 3723.2971584123543 * exp(-2 * 22.7 * r) + 8 * 7891.06583987383 * pow(r, 2.0) * exp(-2 * 9.425 * r) + 8 * 200.02482190951136 * pow(r, 4.0) * exp(-2 * 3.9166666666666665 * r) + 2 * 0.0009014374909556084 * pow(r, 5.4) * exp(-2 * 0.891891891891892 * r) + 3 * 0.17583204016482312 * pow(r, 4.0) * exp(-2 * 1.4333333333333336 * r);
    }

    case 24:
    {
        // Slater wavefunction for atomic number 24
        return 2 * 4237.358075302591 * exp(-2 * 23.7 * r) + 8 * 10218.365897376594 * pow(r, 2.0) * exp(-2 * 9.925 * r) + 8 * 354.31547390255463 * pow(r, 4.0) * exp(-2 * 4.25 * r) + 2 * 0.0013094801520034804 * pow(r, 5.4) * exp(-2 * 0.9324324324324331 * r) + 4 * 0.47103721296316764 * pow(r, 4.0) * exp(-2 * 1.6499999999999997 * r);
    }

    case 25:
    {
        // Slater wavefunction for atomic number 25
        return 2 * 4796.682658008162 * exp(-2 * 24.7 * r) + 8 * 13064.993283967953 * pow(r, 2.0) * exp(-2 * 10.425 * r) + 8 * 601.0875231006424 * pow(r, 4.0) * exp(-2 * 4.583333333333333 * r) + 2 * 0.0018722314140677398 * pow(r, 5.4) * exp(-2 * 0.9729729729729734 * r) + 5 * 1.117210516657861 * pow(r, 4.0) * exp(-2 * 1.8666666666666671 * r);
    }

    case 26:
    {
        // Slater wavefunction for atomic number 26
        return 2 * 5403.180765846168 * exp(-2 * 25.7 * r) + 8 * 16513.384789216943 * pow(r, 2.0) * exp(-2 * 10.925 * r) + 8 * 982.5687970790284 * pow(r, 4.0) * exp(-2 * 4.916666666666667 * r) + 2 * 0.0026380387491179032 * pow(r, 5.4) * exp(-2 * 1.0135135135135134 * r) + 6 * 2.409788611299073 * pow(r, 4.0) * exp(-2 * 2.0833333333333335 * r);
    }

    default:
        return 0.0;
    }
}


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

    max_x = min_x + ceil((max_x - min_x) / dx) * dx;


    min_y -= padding;
    max_y += padding;

    max_y = min_y + ceil((max_y - min_y) / dy) * dy;


    min_z -= padding;
    max_z += padding;
    max_z = min_z + ceil((max_z - min_z) / dz) * dz;


    // Calculate the number of grid points along each direction
    eddi_size_t n_x = (eddi_size_t)((max_x - min_x) / dx + 0.5);
    eddi_size_t n_y = (eddi_size_t)((max_y - min_y) / dy + 0.5);
    eddi_size_t n_z = (eddi_size_t)((max_z - min_z) / dz + 0.5);

    // Set the origin point
    eddi_point_t origin = {min_x, min_y, min_z};
    // printf("Field\n");
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

    const eddi_real_t cutoff_radius = 12.0;
    const eddi_real_t max_radius_squared = cutoff_radius * cutoff_radius;



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
            
            if (radius_2 > max_radius_squared) 
            {
                continue;
            } 
            else
            {
                local_density += molecule->density[atom_idx](sqrt(radius_2), 0.0, 0.0);
            }

        }
        field[idx] = local_density;
    }
}

void eddi_compute_density_field_cl(eddi_density_field_t* density_field, eddi_molecule_t* molecule, eddi_cl_info_t* info)
{
    // Here we also construct the cell list because this is just experimental

    // const size_t n_radii = (size_t) (max_radius_squared);
    // eddi_real_t radii_sqrts[n_radii];
    
    // // #pragma omp parallel for
    // for (int i = 0; i < n_radii; ++i)
    // {
    //     radii_sqrts[i] = sqrt((eddi_real_t)i); 
    // }


    const eddi_real_t c_x = info->cx; 
    const eddi_real_t c_y = info->cy;
    const eddi_real_t c_z = info->cz;

    const eddi_size_t nc_x = density_field->x_size / c_x;
    const eddi_size_t nc_y = density_field->y_size / c_y;
    const eddi_size_t nc_z = density_field->z_size / c_z;
    const eddi_size_t nc_tot = nc_x * nc_y * nc_z;

    const eddi_real_t multiplier = molecule->n_atoms * 0.2;
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
        const eddi_size_t cell_x = (eddi_size_t)floor(((molecule->atoms_x[i] - density_field->origin.x) / c_x));
        const eddi_size_t cell_y = (eddi_size_t)floor(((molecule->atoms_y[i] - density_field->origin.y) / c_y));
        const eddi_size_t cell_z = (eddi_size_t)floor(((molecule->atoms_z[i] - density_field->origin.z) / c_z));

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

    // printf("Cell list with %zu cells. %zu each\n", nc_tot, c_n_atoms);

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


        const eddi_size_t cell_x = (eddi_size_t)floor(((px - density_field->origin.x) / c_x + 1e-6));
        const eddi_size_t cell_y = (eddi_size_t)floor(((py - density_field->origin.y) / c_y + 1e-6));
        const eddi_size_t cell_z = (eddi_size_t)floor(((pz - density_field->origin.z) / c_z + 1e-6));

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
                        local += cells[neighbor_id].density[atom_idx](sqrt(distance_squared), 0.0, 0.0); 
                    }
                }
            }
        }
        density_field->field[p_idx] = local;
    }

    // Free allocated memory for cells
    for (eddi_size_t i = 0; i < nc_tot; ++i)
    {
        free(cells[i].atoms_x);
        free(cells[i].atoms_y);
        free(cells[i].atoms_z);
        free(cells[i].density);
    }
}



void eddi_compute_density_field_cl_opt(eddi_density_field_t* density_field, eddi_molecule_t* molecule, eddi_cl_info_t* info)
{
    // Here we also construct the cell list because this is just experimental
    const eddi_real_t max_radius_squared = 12.0 * 12.0;

    EDDI_DEBUG_PRINT("Constructed radii\n");
    const eddi_real_t c_x = info->cx; 
    const eddi_real_t c_y = info->cy;
    const eddi_real_t c_z = info->cz;

    const eddi_size_t nc_x = density_field->x_size / c_x;
    const eddi_size_t nc_y = density_field->y_size / c_y;
    const eddi_size_t nc_z = density_field->z_size / c_z;
    const eddi_size_t nc_tot = nc_x * nc_y * nc_z;

    const eddi_real_t multiplier = molecule->n_atoms * 0.6;
    const eddi_size_t c_n_atoms = multiplier + 1;
    
    eddi_molecule_t cells[nc_tot];
    eddi_size_t occupancy[nc_tot];

    EDDI_DEBUG_PRINT("Allocated cells and occupancy\n");

    memset(occupancy, 0, sizeof(eddi_size_t) * nc_tot);
    // Allocate per-cell atom arrays
    eddi_real_t* x_values = malloc(sizeof(eddi_real_t) * c_n_atoms * nc_tot);
    eddi_real_t* y_values = malloc(sizeof(eddi_real_t) * c_n_atoms * nc_tot);
    eddi_real_t* z_values = malloc(sizeof(eddi_real_t) * c_n_atoms * nc_tot);
    eddi_atomic_number_t* atomic_numbers_values = malloc(sizeof(eddi_atomic_number_t) * c_n_atoms * nc_tot);

    for (eddi_size_t i = 0; i < nc_tot; ++i) {
        cells[i].atoms_x = x_values + (i * c_n_atoms);
        cells[i].atoms_y = y_values + (i * c_n_atoms);
        cells[i].atoms_z = z_values + (i * c_n_atoms);
        cells[i].atomic_numbers = atomic_numbers_values + (i * c_n_atoms);
    }

    EDDI_DEBUG_PRINT("Initialized cells\n");
    EDDI_DEBUG_PRINT("Adding atoms to cells\n");

    // Now copy the atom data in the cells
    for (eddi_size_t i = 0; i < molecule->n_atoms; ++i)
    {
        EDDI_DEBUG_PRINT("Atom %d\n", i);
        const eddi_size_t cell_x = (eddi_size_t)floor(((molecule->atoms_x[i] - density_field->origin.x) / c_x));
        const eddi_size_t cell_y = (eddi_size_t)floor(((molecule->atoms_y[i] - density_field->origin.y) / c_y));
        const eddi_size_t cell_z = (eddi_size_t)floor(((molecule->atoms_z[i] - density_field->origin.z) / c_z));

        if (cell_x >= nc_x || cell_y >= nc_y || cell_z >= nc_z)
        {
            printf("[WARNING] Atom is outside the domain boundaries and will be ignored.\n");
            exit(1);
            continue;
        }

        const eddi_size_t cell_idx = cell_x * nc_y * nc_z + cell_y * nc_z + cell_z;
        cells[cell_idx].atoms_x[occupancy[cell_idx]] = molecule->atoms_x[i];
        cells[cell_idx].atoms_y[occupancy[cell_idx]] = molecule->atoms_y[i];
        cells[cell_idx].atoms_z[occupancy[cell_idx]] = molecule->atoms_z[i];
        // cells[cell_idx].density[occupancy[cell_idx]] = molecule->density[i];
        cells[cell_idx].atomic_numbers[occupancy[cell_idx]] = molecule->atomic_numbers[i];
        occupancy[cell_idx]++;
    } 

    const eddi_size_t n_points = density_field->x_size * density_field->y_size * density_field->z_size;
    
    #pragma omp parallel for shared(density_field)
    for (eddi_size_t p_idx = 0; p_idx < n_points; ++p_idx)
    {

        const eddi_size_t i_z = p_idx % density_field->z_size;
        const eddi_size_t i_y = (p_idx / density_field->z_size) % density_field->y_size;
        const eddi_size_t i_x = (p_idx / density_field->z_size) / density_field->y_size;

        const eddi_real_t pz = density_field->origin.z + i_z * density_field->dz;
        const eddi_real_t py = density_field->origin.y + i_y * density_field->dy;
        const eddi_real_t px = density_field->origin.x + i_x * density_field->dx;


        const eddi_size_t cell_i = (eddi_size_t)floor(((px - density_field->origin.x) / c_x + 1e-6));
        const eddi_size_t cell_j = (eddi_size_t)floor(((py - density_field->origin.y) / c_y + 1e-6));
        const eddi_size_t cell_k = (eddi_size_t)floor(((pz - density_field->origin.z) / c_z + 1e-6));

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

                    const size_t idx = inner_i * 3 * 3 + inner_j * 3 + inner_k;
                    if (arr[idx] <= max_radius_squared)
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
                                    if (distance_squared <= 12*12.0)
                                    {
                                        partial += slater(cells[neighbor_index].atomic_numbers[idx], sqrt(distance_squared), 0.0, 0.0);

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

    free(x_values);
    free(y_values);
    free(z_values);
    free(atomic_numbers_values);
}


#define SQR(x) ((x) * (x))

void eddi_compute_density_field_cl_opt_v2(eddi_density_field_t* density_field, eddi_molecule_t* molecule, eddi_cl_info_t* info)
{
    const eddi_real_t max_radius_squared = 12.0 * 12.0;
    const eddi_real_t c_x = info->cx;
    const eddi_real_t c_y = info->cy;
    const eddi_real_t c_z = info->cz;

    const eddi_size_t nc_x = density_field->x_size / c_x;
    const eddi_size_t nc_y = density_field->y_size / c_y;
    const eddi_size_t nc_z = density_field->z_size / c_z;
    const eddi_size_t nc_tot = nc_x * nc_y * nc_z;

    const eddi_size_t points_per_block_x = (eddi_size_t)ceil(c_x / density_field->dx);
    const eddi_size_t points_per_block_y = (eddi_size_t)ceil(c_y / density_field->dy);
    const eddi_size_t points_per_block_z = (eddi_size_t)ceil(c_z / density_field->dz);

    const eddi_size_t c_n_atoms = (eddi_size_t)(molecule->n_atoms * 0.6) + 1;

    eddi_molecule_t cells[nc_tot];
    eddi_size_t occupancy[nc_tot];
    memset(occupancy, 0, sizeof(eddi_size_t) * nc_tot);

    // Allocate per-cell atom arrays
    eddi_real_t* x_values = malloc(sizeof(eddi_real_t) * c_n_atoms * nc_tot);
    eddi_real_t* y_values = malloc(sizeof(eddi_real_t) * c_n_atoms * nc_tot);
    eddi_real_t* z_values = malloc(sizeof(eddi_real_t) * c_n_atoms * nc_tot);
    eddi_atomic_number_t* atomic_numbers_values = malloc(sizeof(eddi_atomic_number_t) * c_n_atoms * nc_tot);

    for (eddi_size_t i = 0; i < nc_tot; ++i) {
        cells[i].atoms_x = x_values + (i * c_n_atoms);
        cells[i].atoms_y = y_values + (i * c_n_atoms);
        cells[i].atoms_z = z_values + (i * c_n_atoms);
        cells[i].atomic_numbers = atomic_numbers_values + (i * c_n_atoms);
    }

    // Assign atoms to cells
    for (eddi_size_t i = 0; i < molecule->n_atoms; ++i) {
        eddi_size_t cell_x = (eddi_size_t)((molecule->atoms_x[i] - density_field->origin.x) / c_x);
        eddi_size_t cell_y = (eddi_size_t)((molecule->atoms_y[i] - density_field->origin.y) / c_y);
        eddi_size_t cell_z = (eddi_size_t)((molecule->atoms_z[i] - density_field->origin.z) / c_z);

        if (cell_x >= nc_x || cell_y >= nc_y || cell_z >= nc_z) {
            fprintf(stderr, "[WARNING] Atom %zu outside domain boundaries\n", i);
            continue;
        }

        eddi_size_t idx = cell_x * nc_y * nc_z + cell_y * nc_z + cell_z;
        eddi_size_t occ = occupancy[idx]++;
        cells[idx].atoms_x[occ] = molecule->atoms_x[i];
        cells[idx].atoms_y[occ] = molecule->atoms_y[i];
        cells[idx].atoms_z[occ] = molecule->atoms_z[i];
        cells[idx].atomic_numbers[occ] = molecule->atomic_numbers[i];
    }

    // Loop over blocks
    #pragma omp parallel for
    for (eddi_size_t cell_i = 0; cell_i < nc_x; ++cell_i) {
        for (eddi_size_t cell_j = 0; cell_j < nc_y; ++cell_j) {
            for (eddi_size_t cell_k = 0; cell_k < nc_z; ++cell_k) {
                for (eddi_size_t point_i = 0; point_i < points_per_block_x; ++point_i) {
                    for (eddi_size_t point_j = 0; point_j < points_per_block_y; ++point_j) {
                        for (eddi_size_t point_k = 0; point_k < points_per_block_z; ++point_k) {

                            eddi_size_t gx = cell_i * points_per_block_x + point_i;
                            eddi_size_t gy = cell_j * points_per_block_y + point_j;
                            eddi_size_t gz = cell_k * points_per_block_z + point_k;

                            if (gx >= density_field->x_size || gy >= density_field->y_size || gz >= density_field->z_size)
                                continue;

                            eddi_size_t p_idx = gx * density_field->y_size * density_field->z_size +
                                                gy * density_field->z_size + gz;

                            eddi_real_t px = density_field->origin.x + gx * density_field->dx;
                            eddi_real_t py = density_field->origin.y + gy * density_field->dy;
                            eddi_real_t pz = density_field->origin.z + gz * density_field->dz;

                            eddi_real_t partial = 0.0;

                            for (int ni = -1; ni <= 1; ++ni) {
                                for (int nj = -1; nj <= 1; ++nj) {
                                    for (int nk = -1; nk <= 1; ++nk) {
                                        eddi_size_t ni_idx = cell_i + ni;
                                        eddi_size_t nj_idx = cell_j + nj;
                                        eddi_size_t nk_idx = cell_k + nk;

                                        if (ni_idx >= nc_x || nj_idx >= nc_y || nk_idx >= nc_z)
                                            continue;

                                        eddi_size_t n_idx = ni_idx * nc_y * nc_z + nj_idx * nc_z + nk_idx;

                                        for (eddi_size_t a = 0; a < occupancy[n_idx]; ++a) {
                                            eddi_real_t dx = px - cells[n_idx].atoms_x[a];
                                            eddi_real_t dy = py - cells[n_idx].atoms_y[a];
                                            eddi_real_t dz = pz - cells[n_idx].atoms_z[a];
                                            eddi_real_t r2 = SQR(dx) + SQR(dy) + SQR(dz);

                                            if (r2 <= max_radius_squared) {
                                                eddi_real_t r = sqrt(r2);
                                                partial += slater(cells[n_idx].atomic_numbers[a], r, 0.0, 0.0);
                                            }
                                        }
                                    }
                                }
                            }
                            density_field->field[p_idx] = partial;
                        }
                    }
                }
            }
        }
    }

    free(x_values);
    free(y_values);
    free(z_values);
    free(atomic_numbers_values);
}



void eddi_compute_density_field_atom(eddi_density_field_t* density_field, eddi_molecule_t* molecule)
{
    const eddi_real_t cutoff_radius = 12.0;
    const eddi_real_t max_radius_squared = cutoff_radius * cutoff_radius;
    const size_t n_radii = (size_t) (max_radius_squared);
    
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
        const eddi_real_t starting_x = ox + floor((molecule->atoms_x[atom_id] - cutoff_radius - ox) / dx) * dx;
        const eddi_real_t starting_y = oy + floor((molecule->atoms_y[atom_id] - cutoff_radius - oy) / dy) * dy;
        const eddi_real_t starting_z = oz + floor((molecule->atoms_z[atom_id] - cutoff_radius - oz) / dz) * dz;

        const eddi_real_t end_x = molecule->atoms_x[atom_id] + cutoff_radius;
        const eddi_real_t end_y = molecule->atoms_y[atom_id] + cutoff_radius;
        const eddi_real_t end_z = molecule->atoms_z[atom_id] + cutoff_radius;

        for (eddi_real_t x = starting_x;  x <= end_x; x += dx)
        {
            for (eddi_real_t y = starting_y; y <= end_y; y += dy)
            {
                for (eddi_real_t z = starting_z; z <= end_z; z += dz)
                {
                    const eddi_size_t idx_x = (eddi_size_t)floor(((x - ox) / dx + 1e-6));
                    const eddi_size_t idx_y = (eddi_size_t)floor(((y - oy) / dy + 1e-6));
                    const eddi_size_t idx_z = (eddi_size_t)floor(((z - oz) / dz + 1e-6));

                    if (idx_x < n_x && idx_y < n_y && idx_z < n_z)
                    {
                        const eddi_size_t idx = idx_x * n_y * n_z + idx_y * n_z + idx_z;
                        #pragma omp atomic
                        density_field->field[idx] += slater(molecule->atomic_numbers[atom_id], sqrt((x - molecule->atoms_x[atom_id]) * (x - molecule->atoms_x[atom_id]) +
                                                                                    (y - molecule->atoms_y[atom_id]) * (y - molecule->atoms_y[atom_id]) +
                                                                                    (z - molecule->atoms_z[atom_id]) * (z - molecule->atoms_z[atom_id])), 0.0, 0.0);
                    }
                }
            }
        }
        
    }
}

double eddi_compute_volume(eddi_density_field_t* density_field, eddi_real_t isodensity)
{
    eddi_size_t nx = density_field->x_size;
    eddi_size_t ny = density_field->y_size;
    eddi_size_t nz = density_field->z_size;

    int count = 0;
    for (int i = 0; i < nx-1; ++i)
    {
        for (int j = 0; j < ny-1; ++j)
        {
            for (int k = 0; k < nz-1; ++k)
            {
                // Check if the cube with origin in i, j, k has all vertices inside the isodensity
                bool all_inside = true;
                for (int di = 0; di <= 1; ++di)
                {
                    for (int dj = 0; dj <= 1; ++dj)
                    {
                        for (int dk = 0; dk <= 1; ++dk)
                        {
                            eddi_size_t vertex_idx = (i + di) * ny * nz + (j + dj) * nz + (k + dk);
                            if (density_field->field[vertex_idx] <= isodensity)
                            {
                                all_inside = false;
                                break;
                            }
                        }
                        if (!all_inside) break;
                    }
                    if (!all_inside) break;
                }

                if (all_inside)
                {
                    count ++;
                }
            }
        }
    }
    return count * density_field->dx * density_field->dy * density_field->dz;
}

void eddi_free_density_field(eddi_density_field_t* density_field)
{
    // Only need to deallocate the field
    free(density_field->field);
}