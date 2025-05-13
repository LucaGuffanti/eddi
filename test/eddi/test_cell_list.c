
#include "eddi_molecule_readers.h"
#include "eddi_density_writers.h"

#include "time.h"
int main(int argc, char** argv)
{
    eddi_molecule_t molecule;
    eddi_density_field_t density_field;

    printf("Sizeof eddi_real_t %d\n (float: %d double: %d)\n", sizeof(eddi_real_t), sizeof(float), sizeof(double));
    const eddi_real_t resolution = 0.4;

    printf("Starting\n");
    eddi_read_pdb(argv[1], &molecule);
    printf("Init\n");
    eddi_init_field_from_molecule(&density_field, &molecule, 10.0, resolution, resolution, resolution);

    
    clock_t start, end;
    double cpu_time_used;
    printf("Data in\n");
    printf("Atom\n");
    start = clock();
    eddi_compute_density_field_atom(&density_field, &molecule);
    end = clock();
    double v_atom = eddi_compute_volume(&density_field, 0.004);

    eddi_cl_info_t info = {.cx = 12.0, .cy = 12.0, .cz = 12.0};

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_atom execution time: %f seconds\n", cpu_time_used);
    eddi_write_gaussian_cube(argv[2], &density_field, &molecule);
    
    
    memset(density_field.field, 0, density_field.x_size * density_field.y_size * density_field.z_size);
    
    start = clock();
    eddi_compute_density_field_cl(&density_field, &molecule, &info);
    end = clock();
    double v_cl = eddi_compute_volume(&density_field, 0.004);

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_cl execution time: %f seconds\n", cpu_time_used);
    eddi_write_gaussian_cube(argv[3], &density_field, &molecule);
    
    memset(density_field.field, 0, density_field.x_size * density_field.y_size * density_field.z_size);


    start = clock();
    eddi_compute_density_field_cl_opt(&density_field, &molecule, &info);
    end = clock();
    double v_opt = eddi_compute_volume(&density_field, 0.004);
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_cl_opt execution time: %f seconds\n", cpu_time_used);
    

    memset(density_field.field, 0, density_field.x_size * density_field.y_size * density_field.z_size);


    start = clock();
    eddi_compute_density_field_cl_opt_v2(&density_field, &molecule, &info);
    end = clock();
    double v_opt2 = eddi_compute_volume(&density_field, 0.004);
    eddi_write_gaussian_cube(argv[4], &density_field, &molecule);
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("eddi_compute_density_field_cl_opt_v2 execution time: %f seconds\n", cpu_time_used);




    eddi_write_gaussian_cube(argv[4], &density_field, &molecule);

    printf("Saved execution of eddi_compute_density_field_atom in %s\n", argv[2]);
    printf("Saved execution of eddi_compute_density_field_cl in %s\n", argv[3]);
    printf("Saved execution of eddi_compute_density_field_cl_opt_v2 in %s\n", argv[4]);


    printf("Atom volume : %lf\n", v_atom);

    printf("Cl volume   : %lf\n", v_cl);
    printf("Rel. err    : %lf (%lf \%)\n ", fabs(v_atom - v_cl)/v_atom, fabs(v_atom - v_cl)/v_atom * 100);
    printf("Opt volume  : %lf\n", v_opt);
    printf("Rel. err    : %lf (%lf \%)\n", fabs(v_atom - v_opt)/v_atom, fabs(v_atom - v_opt)/v_atom * 100);
    printf("Opt volume  : %lf\n", v_opt2);
    printf("Rel. err    : %lf (%lf \%)\n", fabs(v_atom - v_opt2)/v_atom, fabs(v_atom - v_opt2)/v_atom * 100);



    eddi_write_binary(argv[5], &density_field);

}