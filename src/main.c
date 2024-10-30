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
  
  body_system *system_history = NULL;
  body_system system_status;
  if (argc < 4){
    if( rank == 0 ) fprintf(stderr, "Error: using '''%s''' as <n_of_bodies> <n_of_iter> <delta_t>\n", argv[0]);
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
  if (*endptr != '\0' || n_of_iter < 1000) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number of iter, it must be >= 1000. Aborting...\n");
    goto cleanup;
  }
  
  delta_t = (double) strtod(argv[3], &endptr);
  if (*endptr != '\0' || delta_t <= 0.0) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number delta_t, it must be > 0. Aborting...\n");
    goto cleanup;
  }
  // int split_rank;
  // size_t large_body_count, small_body_count;
  // COMPUTE_BODY_COUNT(n_of_bodies, comm_sz, split_rank, large_body_count, small_body_count);
  if (rank == 0){
    body_system *system_history = (body_system *)malloc(SAVE_HISTORY * sizeof(body_system));
    if (system_history == NULL) goto cleanup;
  }

  system_status.mass = (double *) malloc(n_of_bodies * sizeof(double));
  system_status.x_data = (double *) malloc(3 * n_of_bodies * sizeof(double));
  system_status.y_data = (double *) malloc(3 * n_of_bodies * sizeof(double));
  if (system_status.mass == NULL || system_status.mass == NULL || system_status.mass == NULL){
    goto cleanup;
  }

  for (int iter = 0; iter < n_of_iter; iter++){
    acceleration_update(system_status.x_data, system_status.mass, n_of_bodies);
    acceleration_update(system_status.y_data, system_status.mass, n_of_bodies);
    time_step_update(system_status.x_data, n_of_bodies, delta_t);
    time_step_update(system_status.y_data, n_of_bodies, delta_t);

  }
  
  

  if (rank == 0) free(system_history);
  free(system_status.mass);
  free(system_status.x_data);
  free(system_status.y_data);

  MPI_Finalize();

  return 0;

cleanup:
  if ( NULL != system_history) free(system_history);
  if ( NULL != system_status.mass) free(system_status.mass);
  if ( NULL != system_status.x_data) free(system_status.x_data);
  if ( NULL != system_status.y_data) free(system_status.y_data);
  MPI_Abort(comm, 1);
}

