#include "eddi_molecule_readers.h"
#include "eddi_density_writers.h"

#include "time.h"
int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t density_field;

    const eddi_real_t resolution = 0.5;

    eddi_read_pdb(argv[1], &molecule);
    eddi_init_field_from_molecule(&density_field, &molecule, 4.0, resolution, resolution, resolution);

    
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    eddi_compute_density_field_cl(&density_field, &molecule);
    end = clock();

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_cl execution time: %f seconds\n", cpu_time_used);
    eddi_write_gaussian_cube(argv[2], &density_field, &molecule);

    memset(density_field.field, 0, density_field.x_size * density_field.y_size * density_field.z_size * sizeof(eddi_real_t));

    start = clock();
    eddi_compute_density_field(&density_field, &molecule);
    end = clock();

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field execution time: %f seconds\n", cpu_time_used);
    eddi_write_gaussian_cube(argv[3], &density_field, &molecule);

}