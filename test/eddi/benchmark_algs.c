#include "eddi_molecule_readers.h"
#include "eddi_density_writers.h"

struct timespec start_w, finish_w;

int main(int argc, char** argv)
{
    const char* file_paths[] = {
        "/home/lucaguf/thesis/code/data/benchmark/134d_793.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/101m_1413.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a3y_3385.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/10gs_3521.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a00_4770.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a3q_5748.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a6z_6099.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a0l_7872.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1bj1_8694.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1al0_9521.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a0t_10077.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1avo_11361.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1ahj_12832.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/13pk_13271.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a0c_15116.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1be3_16222.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1afr_17326.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1aa1_18956.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/8spx_22000.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1a2v_33726.pdb",
        "/home/lucaguf/thesis/code/data/benchmark/1e7p_37072.pdb",
        NULL
    };

    const int atom_counts[] = {
        793,
        1413,
        3385,
        3521,
        4770,
        5748,
        6099,
        7872,
        8694,
        9521,
        10077,
        11361,
        12832,
        13271,
        15116,
        16222,
        17326,
        18956,
        22000,
        33726,
        37072,
    };
    const int file_count = 21;
    const int iterations_per_run = 1;
    int id = 0;

    const eddi_real_t resolutions[] = {0.8, 0.7, 0.6, 0.5, 0.4};
    const eddi_real_t padding = 15.0;
    const eddi_cl_info_t infos[] = {
        {.cx = 8.0 , .cy = 8.0 , .cz = 8.0 },
        {.cx = 10.0, .cy = 10.0, .cz = 10.0},
        {.cx = 12.0, .cy = 12.0, .cz = 12.0},
    };

    FILE* fp = fopen("benchmark.csv", "w");

    int ompthreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp master
        {
            ompthreads = omp_get_num_threads();
        }
    }
#endif
    
    fprintf(fp, "id,atom_count,mode,cpu_time,wall_time,res,cx\n");
    fflush(fp);
    
    
    printf("Benchmarking of CPU algorithms\n");
    printf("Using %d omp threads.\n", ompthreads);

    for (int file_id = 0; file_id < 1; ++file_id)
    {
        printf("-> File: %s\n", file_paths[file_id]);
        eddi_molecule_t molecule;
        eddi_density_field_t field;

        eddi_read_pdb(file_paths[file_id], &molecule);
        for (int resolution_id = 0; resolution_id < 1; ++resolution_id)
        {
            for(int info_id = 0; info_id < 2; ++info_id)
            {
                const eddi_real_t res = resolutions[resolution_id];
                eddi_cl_info_t info = infos[info_id];

                printf("resolution: %f\n", res);
                eddi_init_field_from_molecule(&field, &molecule, padding, res, res, res);

                printf("    Cell-list\n");
                for (int it = 0; it < iterations_per_run; ++it)
                {
                    // Time the cl iteration
                    printf("    it %d/%d\r", it+1, iterations_per_run);
                    fflush(stdout);
                    

                    clock_t start = clock();
                    clock_gettime(CLOCK_MONOTONIC, &start_w);
                    eddi_compute_density_field_cl(&field, &molecule, &info);
                    clock_gettime(CLOCK_MONOTONIC, &finish_w);
                    clock_t end = clock();

                    memset(field.field, 0, field.x_size * field.y_size * field.z_size * sizeof(eddi_real_t));

                    
                    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
                    double walltime = (finish_w.tv_sec - start_w.tv_sec);
                    walltime += (finish_w.tv_nsec - start_w.tv_nsec) / 1000000000.0;
                    fprintf(fp, "%d,%d,cl,%f,%f,%f,%f\n", id, atom_counts[file_id], time_taken, walltime, res, info.cx);

                    fflush(fp);

                    id++;
                }
                
                // printf("    Optimized Cell-list\n");
                // for (int it = 0; it < iterations_per_run; ++it)
                // {
                //     // Time the cl_opt iteration
                //     printf("    it %d/%d\r", it+1, iterations_per_run);
                //     fflush(stdout);
                    

                //     clock_t start = clock();
                //     clock_gettime(CLOCK_MONOTONIC, &start_w);
                //     eddi_compute_density_field_cl_opt(&field, &molecule, &info);
                //     clock_gettime(CLOCK_MONOTONIC, &finish_w);
                //     clock_t end = clock();
                    
                //     memset(field.field, 0, field.x_size * field.y_size * field.z_size * sizeof(eddi_real_t));

                //     double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
                //     double walltime = (finish_w.tv_sec - start_w.tv_sec);
                //     walltime += (finish_w.tv_nsec - start_w.tv_nsec) / 1000000000.0;
                //     fprintf(fp, "%d,%d,cl,%f,%f,%f,%f\n", id, atom_counts[file_id], time_taken, walltime, res, info.cx);

                //     fflush(fp);


                //     id++;
                // }
                
                printf("    Atom-based\n");
                for (int it = 0; it < iterations_per_run; ++it)
                {
                    // Time the atom-based iteration
                    printf("    it %d/%d\r", it+1, iterations_per_run);
                    fflush(stdout);

                    clock_t start = clock();
                    clock_gettime(CLOCK_MONOTONIC, &start_w);
                    eddi_compute_density_field_atom(&field, &molecule);
                    clock_gettime(CLOCK_MONOTONIC, &finish_w);
                    clock_t end = clock();
                    
                    memset(field.field, 0, field.x_size * field.y_size * field.z_size * sizeof(eddi_real_t));


                    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
                    double walltime = (finish_w.tv_sec - start_w.tv_sec);
                    walltime += (finish_w.tv_nsec - start_w.tv_nsec) / 1000000000.0;
                    fprintf(fp, "%d,%d,atom,%f,%f,%f,0.0\n", id, atom_counts[file_id], time_taken, walltime, res);
                    fflush(fp);

                    id++;
                }
                eddi_free_density_field(&field);
            }
        }
        eddi_free_molecule(&molecule);
    }
}