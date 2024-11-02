#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "tools_simulation.h"
#include "utils.h"

int main(int argc, char** argv){
  MPI_Init(NULL, NULL);
  int rank, comm_sz;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &rank);
  
  body_system system_status;
  double* store_buffer = NULL;
  int ret;
  size_t n_of_bodies, n_of_iter;
  char* filename;
  ret = set_inputs(argc, argv, &n_of_bodies, &n_of_iter, &filename);
  if ( ret == -1 ) goto cleanup;

  int *count, *disp;
  count = (int*) malloc(comm_sz * sizeof(int));
  disp  = (int*) malloc(comm_sz * sizeof(int));
  compute_body_count(n_of_bodies, comm_sz, count, disp);

  if (rank == 0 ){
    store_buffer = (double*) malloc(2 * TMP_BUF_SIZE * n_of_bodies * sizeof(double));
    if ( store_buffer == NULL ) goto cleanup;
  }

  system_status.mass = (double *) malloc(n_of_bodies * sizeof(double));
  system_status.pos = (double *) malloc(2 * n_of_bodies * sizeof(double));
  system_status.vel = (double *) malloc(2 * n_of_bodies * sizeof(double));
  system_status.acc = (double *) calloc(2 * n_of_bodies, sizeof(double));
  if (system_status.mass == NULL || system_status.pos == NULL ||
    system_status.vel == NULL || system_status.acc == NULL){
    goto cleanup;
  }
  
  srand(time(0));
  get_init_ranges(n_of_bodies);
  if (rank == 0){
    set_initial_conditions(&system_status, n_of_bodies);
    print_status(&system_status, n_of_bodies);
  }

  ret = MPI_Bcast(system_status.mass, n_of_bodies, MPI_DOUBLE, 0, comm);
  if (ret != MPI_SUCCESS) goto cleanup;
  ret = MPI_Bcast(system_status.pos, 2 * n_of_bodies, MPI_DOUBLE, 0, comm);
  if (ret != MPI_SUCCESS) goto cleanup;
  ret = MPI_Bcast(system_status.vel, 2 * n_of_bodies, MPI_DOUBLE, 0, comm);
  if (ret != MPI_SUCCESS) goto cleanup;

  int print_iter = 0;
  double delta_t, elapsed_time = 0;
  for (int iter = 0; iter < n_of_iter; iter++){
    // delta_t is calculated at each iteration based on highest velocity body
    delta_t = compute_partial_delta_t(system_status.vel, count[rank], disp[rank]);
    ret = MPI_Allreduce(MPI_IN_PLACE, &delta_t, 1, MPI_DOUBLE, MPI_MIN, comm);
    if (ret != MPI_SUCCESS) goto cleanup;
    elapsed_time += delta_t;

    // rank 0 accumulates data every PRINT_INTERVAL elapsed time
    if(rank == 0 && elapsed_time > print_iter * PRINT_INTERVAL) { 
      accumulate_data(store_buffer, print_iter % TMP_BUF_SIZE, n_of_bodies, &system_status);
      print_iter++;

      // write the data to disk every when the buffer is full
      if( print_iter % TMP_BUF_SIZE == 0) {
        ret = write_data_to_disk(store_buffer, system_status.mass, n_of_bodies, print_iter, filename);
        if (ret == -1) goto cleanup;
      }
    }

    compute_new_accelerations(system_status.mass, system_status.pos, system_status.acc, n_of_bodies,
                        count[rank], disp[rank], NEWTON);

    time_step_update(system_status.pos, system_status.vel, system_status.acc,
                       n_of_bodies, delta_t, count[rank], disp[rank]);
    
    // Synchronize only the positions, acc ad velocity remains private data of the process
    ret = MPI_Allgatherv(MPI_IN_PLACE, count[rank], MPI_DOUBLE,
                  system_status.pos, count, disp, MPI_DOUBLE, comm);
    if (ret != MPI_SUCCESS) goto cleanup;
  }
  
  if(rank == 0) printf("\nLast delta_t is: %.5f\n TOTAL ELAPSED (simulation) TIME:%.3f\n", delta_t, elapsed_time);
  
  if (rank == 0) free(store_buffer);
  free(system_status.mass);
  free(system_status.pos);
  free(system_status.vel);
  free(system_status.acc);
  free(count);
  free(disp);

  MPI_Finalize();

  return 0;

cleanup:
  if ( NULL != store_buffer && rank == 0) free(store_buffer);
  if ( NULL != system_status.mass) free(system_status.mass);
  if ( NULL != system_status.pos) free(system_status.pos);
  if ( NULL != system_status.vel) free(system_status.vel);
  if ( NULL != system_status.acc) free(system_status.acc);
  MPI_Abort(comm, 1);
}

