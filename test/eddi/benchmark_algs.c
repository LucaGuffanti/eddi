#include "eddi_molecule_readers.h"
#include "eddi_density_writers.h"

struct timespec start_w, finish_w;

int main(int argc, char** argv)
{
        
    if (argc != 5)
    {
        fprintf(stderr, "Usage %s <path> <resolution> <cell-list-res> <padding>\n", argv[0]);
        exit(1);
    }

    char* file_path = argv[1];

    eddi_real_t res;
    eddi_real_t cl_res;
    eddi_real_t padding;

#ifdef EDDI_HIGH_PRECISION
    res = strtod(argv[2], NULL);
    cl_res = strtod(argv[3], NULL);
    padding = strtod(argv[4], NULL);
#else
    res = strtof(argv[2], NULL);
    cl_res = strtof(argv[3], NULL);
    padding = strtof(argv[4], NULL);
#endif



    eddi_molecule_t molecule;
    eddi_density_field_t field;

    eddi_read_pdb(file_path, &molecule);
    eddi_init_field_from_molecule(&field, &molecule, padding, res, res, res);
    eddi_cl_info_t info = {.cx = cl_res, .cy = cl_res, .cz = cl_res};

    // Cell-list    
    clock_t start = clock();
    clock_gettime(CLOCK_MONOTONIC, &start_w);
    eddi_compute_density_field_cl(&field, &molecule, &info);
    clock_gettime(CLOCK_MONOTONIC, &finish_w);
    clock_t end = clock();

    memset(field.field, 0, field.x_size * field.y_size * field.z_size * sizeof(eddi_real_t));
    
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    double walltime = (finish_w.tv_sec - start_w.tv_sec);
    walltime += (finish_w.tv_nsec - start_w.tv_nsec) / 1000000000.0;
    printf("cl %f %f\n", time_taken, walltime);



    // Cell-list optimized
    start = clock();
    clock_gettime(CLOCK_MONOTONIC, &start_w);
    eddi_compute_density_field_cl_opt(&field, &molecule, &info);
    clock_gettime(CLOCK_MONOTONIC, &finish_w);
    end = clock();
        
    memset(field.field, 0, field.x_size * field.y_size * field.z_size * sizeof(eddi_real_t));

    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    walltime = (finish_w.tv_sec - start_w.tv_sec);
    walltime += (finish_w.tv_nsec - start_w.tv_nsec) / 1000000000.0;
    printf("opt %f %f\n", time_taken, walltime);



    // Atom
    start = clock();
    clock_gettime(CLOCK_MONOTONIC, &start_w);
    eddi_compute_density_field_atom(&field, &molecule);
    clock_gettime(CLOCK_MONOTONIC, &finish_w);
    end = clock();
    
    memset(field.field, 0, field.x_size * field.y_size * field.z_size * sizeof(eddi_real_t));


    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    walltime = (finish_w.tv_sec - start_w.tv_sec);
    walltime += (finish_w.tv_nsec - start_w.tv_nsec) / 1000000000.0;
    printf("atom %f %f\n", time_taken, walltime);



    eddi_free_density_field(&field);
    eddi_free_molecule(&molecule);
    eddi_free_atomic_data();
}