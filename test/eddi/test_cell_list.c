
#include "eddi_molecule_readers.h"
#include "eddi_density_writers.h"

#include "time.h"
int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t density_field;

    printf("Sizeof eddi_real_t %d\n (float: %d double: %d)\n", sizeof(eddi_real_t), sizeof(float), sizeof(double));
    const eddi_real_t resolution = 0.8;

    printf("Starting\n");
    eddi_read_pdb(argv[1], &molecule);
    printf("Init\n");
    eddi_init_field_from_molecule(&density_field, &molecule, 40.0, resolution, resolution, resolution);

    
    clock_t start, end;
    double cpu_time_used;
    printf("Data in\n");
    printf("Atom\n");
    start = clock();
    eddi_compute_density_field_atom(&density_field, &molecule);
    end = clock();

    eddi_cl_info_t info = {.cx = 10.0, .cy = 10.0, .cz = 10.0};

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_atom execution time: %f seconds\n", cpu_time_used);
    eddi_write_gaussian_cube(argv[2], &density_field, &molecule);

    memset(density_field.field, 0, density_field.x_size * density_field.y_size * density_field.z_size);

    start = clock();
    eddi_compute_density_field_cl(&density_field, &molecule, &info);
    end = clock();

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_cl execution time: %f seconds\n", cpu_time_used);
    eddi_write_gaussian_cube(argv[3], &density_field, &molecule);
    
    memset(density_field.field, 0, density_field.x_size * density_field.y_size * density_field.z_size);


    start = clock();
    eddi_compute_density_field_cl_opt(&density_field, &molecule, &info);
    end = clock();

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_cl_opt execution time: %f seconds\n", cpu_time_used);

    eddi_write_gaussian_cube(argv[4], &density_field, &molecule);

    printf("Saved execution of eddi_compute_density_field_atom in %s\n", argv[2]);
    printf("Saved execution of eddi_compute_density_field_cl in %s\n", argv[3]);
    printf("Saved execution of eddi_compute_density_field_cl_opt in %s\n", argv[4]);


    eddi_write_binary(argv[5], &density_field);

}