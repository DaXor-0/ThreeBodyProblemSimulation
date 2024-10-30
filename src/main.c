#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "tools.h"

int main(int argc, char** argv){
  MPI_Init(NULL, NULL);
  int rank, comm_sz;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &rank);
  
  body_system *system_history = NULL, system_status;
  if (argc < 3){
    if( rank == 0 ) fprintf(stderr, "Error: using '''%s''' as <n_of_bodies> <n_of_iter>\n", argv[0]);
    goto cleanup;
  }

  size_t n_of_bodies, n_of_iter;
  double delta_t;
  char *endptr;
  n_of_bodies = (size_t) strtol(argv[1], &endptr, 10);
  if (*endptr != '\0' || n_of_bodies < 3) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number of bodies, it must be >= 3. Aborting...\n");
    goto cleanup;
  }

  n_of_iter = (size_t) strtol(argv[2], &endptr, 10);
  if (*endptr != '\0' || n_of_iter < 100) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number of iter, it must be >= 1000. Aborting...\n");
    goto cleanup;
  } 

  // NOTE: delta_t should be (?) something like 
  // grid_size/(6*num of bodies*average mass), average mass is 1
  // just for safety we set it up to be one tenth smaller
  delta_t  = (GRID_MAX - GRID_MIN) / (80 * (double) n_of_bodies);

  int split_rank;
  size_t large_body_count, small_body_count;
  COMPUTE_BODY_COUNT(n_of_bodies, comm_sz, split_rank, large_body_count, small_body_count);
  int *count, *disp;
  count = (int*) malloc(comm_sz * sizeof(int));
  disp = (int*) malloc(comm_sz * sizeof(int));
  for (int i = 0; i < comm_sz; i++){
    count[i] = (i < split_rank) ? (int)large_body_count * 3 : (int)small_body_count * 3;
    disp[i] = (i == 0) ? 0 : disp[i - 1] + count [i - 1];
  }

  // if (rank == 0){
  //   body_system *system_history = (body_system *)malloc(SAVE_HISTORY * sizeof(body_system));
  //   if (system_history == NULL) goto cleanup;
  // }

  system_status.mass = (double *) malloc(n_of_bodies * sizeof(double));
  system_status.x_data = (double *) malloc(3 * n_of_bodies * sizeof(double));
  system_status.y_data = (double *) malloc(3 * n_of_bodies * sizeof(double));
  if (system_status.mass == NULL || system_status.x_data == NULL || system_status.y_data == NULL){
    goto cleanup;
  }
  
  if (rank == 0){
    set_initial_conditions(&system_status, n_of_bodies);
      printf("iter_number,body_id, mass, x_pos, x_vel, x_acc, y_pos, y_vel, y_acc\n"); 
      for (int i = 0; i < n_of_bodies; i++)
      printf("0,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", i, system_status.mass[i],
          system_status.x_data[i*3], system_status.x_data[i*3+1], system_status.x_data[i*3+2],
          system_status.y_data[i*3], system_status.y_data[i*3+1], system_status.y_data[i*3+2]);
  }

  MPI_Bcast(system_status.mass, n_of_bodies, MPI_DOUBLE, 0, comm);
  MPI_Bcast(system_status.x_data, 3 * n_of_bodies, MPI_DOUBLE, 0, comm);
  MPI_Bcast(system_status.y_data, 3 * n_of_bodies, MPI_DOUBLE, 0, comm);

  // const char *filename = "simulation_output.csv";

  for (int iter = 0; iter < n_of_iter; iter++){
    acceleration_update(system_status.x_data, system_status.y_data, system_status.mass, n_of_bodies, count[rank], disp[rank]);
    
    time_step_update(system_status.x_data, n_of_bodies, delta_t, count[rank], disp[rank]);
    time_step_update(system_status.y_data, n_of_bodies, delta_t, count[rank], disp[rank]);

    MPI_Allgatherv(MPI_IN_PLACE, count[rank], MPI_DOUBLE,
                  system_status.x_data, count, disp, MPI_DOUBLE, comm);
    MPI_Allgatherv(MPI_IN_PLACE, count[rank], MPI_DOUBLE,
                  system_status.y_data, count, disp, MPI_DOUBLE, comm);
    
    if (rank == 0) {
      for (int i = 0; i < n_of_bodies; i++)
        printf("%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", iter+1, i, system_status.mass[i],
          system_status.x_data[i*3], system_status.x_data[i*3+1], system_status.x_data[i*3+2],
          system_status.y_data[i*3], system_status.y_data[i*3+1], system_status.y_data[i*3+2]);
    }
    // if(iter % SAVE_HISTORY == SAVE_HISTORY - 1 && rank == 0){
    //   print_data(filename, &system_history, n_of_bodies, iter, (iter < SAVE_HISTORY) ? 1 : 0);
    // }
  }
  
  

  // if (rank == 0) free(system_history);
  free(system_status.mass);
  free(system_status.x_data);
  free(system_status.y_data);
  free(count);
  free(disp);

  MPI_Finalize();

  return 0;

cleanup:
  // if ( NULL != system_history) free(system_history);
  if ( NULL != system_status.mass) free(system_status.mass);
  if ( NULL != system_status.x_data) free(system_status.x_data);
  if ( NULL != system_status.y_data) free(system_status.y_data);
  MPI_Abort(comm, 1);
}

