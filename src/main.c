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
  if (*endptr != '\0' || n_of_iter < 1000) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number of iter, it must be >= 1000. Aborting...\n");
    goto cleanup;
  }
  
  // delta_t = (double) strtod(argv[3], &endptr);
  // if (*endptr != '\0' || delta_t <= 0.0) {
  //   if( rank == 0 ) fprintf(stderr, "Error: Invalid number delta_t, it must be > 0. Aborting...\n");
  //   goto cleanup;
  // }

  // WARN: need to set up later, now is like this to remove compiler warnings
  // delta_t should be (?) something like grid_size/(6*num of bodies*average mass), average mass is 1
  // just for safety we set it up to be one tenth smaller
  delta_t  = (GRID_MAX-GRID_MIN) / (60 * n_of_bodies);

  int split_rank;
  size_t large_body_count, small_body_count, my_count;
  COMPUTE_BODY_COUNT(n_of_bodies, comm_sz, split_rank, large_body_count, small_body_count);
  my_count = (rank < split_rank) ? large_body_count: small_body_count;
  
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
  
  if (rank == 0){
    set_initial_conditions(&system_history[0], n_of_bodies);

    MPI_Bcast(system_history[0].mass, n_of_bodies, MPI_DOUBLE, 0, comm);
    MPI_Bcast(system_history[0].x_data, 3 * n_of_bodies, MPI_DOUBLE, 0, comm);
    MPI_Bcast(system_history[0].y_data, 3 * n_of_bodies, MPI_DOUBLE, 0, comm);
  } 
  else {
    MPI_Bcast(system_status.mass, n_of_bodies, MPI_DOUBLE, 0, comm);
    MPI_Bcast(system_status.x_data, 3 * n_of_bodies, MPI_DOUBLE, 0, comm);
    MPI_Bcast(system_status.y_data, 3 * n_of_bodies, MPI_DOUBLE, 0, comm);
  }

  const char *filename = "simulation_output.csv";

  for (int iter = 0; iter < n_of_iter; iter++){
    if (rank == 0){
      acceleration_update(system_history[iter % SAVE_HISTORY].x_data,
                          system_history[iter % SAVE_HISTORY].mass, n_of_bodies);
      acceleration_update(system_history[iter % SAVE_HISTORY].y_data,
                          system_history[iter % SAVE_HISTORY].mass, n_of_bodies);
      
      time_step_update(system_history[iter % SAVE_HISTORY].x_data, n_of_bodies, delta_t);
      time_step_update(system_history[iter % SAVE_HISTORY].y_data, n_of_bodies, delta_t);

      MPI_Allgather(system_history[iter % SAVE_HISTORY].x_data, 3 * my_count, MPI_FLOAT,
                    system_history[iter % SAVE_HISTORY].x_data, 3 * n_of_bodies, MPI_FLOAT, comm);
      MPI_Allgather(system_history[iter % SAVE_HISTORY].y_data, 3 * my_count, MPI_FLOAT,
                    system_history[iter % SAVE_HISTORY].y_data, 3 * n_of_bodies, MPI_FLOAT, comm);
    }
    else{
      acceleration_update(system_status.x_data, system_status.mass, n_of_bodies);
      acceleration_update(system_status.y_data, system_status.mass, n_of_bodies);
      
      time_step_update(system_status.x_data, n_of_bodies, delta_t);
      time_step_update(system_status.y_data, n_of_bodies, delta_t);

      MPI_Allgather(system_status.x_data, 3 * my_count, MPI_FLOAT,
                    system_status.x_data, 3 * n_of_bodies, MPI_FLOAT, comm);
      MPI_Allgather(system_status.y_data, 3 * my_count, MPI_FLOAT,
                    system_status.y_data, 3 * n_of_bodies, MPI_FLOAT, comm);
    }
    
    if(iter % SAVE_HISTORY == SAVE_HISTORY - 1 && rank == 0){
      print_data(filename, &system_history, n_of_bodies, iter, (iter < SAVE_HISTORY) ? 1 : 0);
    }
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

