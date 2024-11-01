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
  
  body_system store_buffer, system_status;
  int ret;
  size_t n_of_bodies, n_of_iter;
  char* filename;
  ret = set_inputs(argc, argv, &n_of_bodies, &n_of_iter, &filename);
  if ( ret == -1 ) goto cleanup;

  int split_rank, *count, *disp;
  size_t large_body_count, small_body_count;
  COMPUTE_BODY_COUNT(n_of_bodies, comm_sz, split_rank, large_body_count, small_body_count);
  count = (int*) malloc(comm_sz * sizeof(int));
  disp  = (int*) malloc(comm_sz * sizeof(int));
  for (int i = 0; i < comm_sz; i++){
    count[i] = (i < split_rank) ? (int) large_body_count * 6 : (int) small_body_count * 6;
    disp[i]  = (i == 0) ? 0 : disp[i - 1] + count [i - 1];
  }

  if (rank == 0 ){
    ret = allocate_store_buffer(&store_buffer, n_of_bodies);
    if ( ret == -1 ) goto cleanup;
  }

  system_status.mass = (double *) malloc(n_of_bodies * sizeof(double));
  system_status.data = (double *) malloc(6 * n_of_bodies * sizeof(double));
  if (system_status.mass == NULL || system_status.data == NULL){
    goto cleanup;
  }
  
  double vel_range;
  srand(time(0));
  if (rank == 0){
    set_initial_conditions(&system_status, n_of_bodies, &vel_range);
    print_status(&system_status, n_of_bodies);
  }

  MPI_Bcast(system_status.mass, n_of_bodies, MPI_DOUBLE, 0, comm);
  MPI_Bcast(system_status.data, 6 * n_of_bodies, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&vel_range, 1, MPI_DOUBLE, 0, comm);

  int print_iter = 0;
  double delta_t, elapsed_time = 0;
  for (int iter = 0; iter < n_of_iter; iter++){
    delta_t = compute_new_delta_t(system_status.data, n_of_bodies);
    elapsed_time += delta_t;
    if(rank == 0){
      // accumulate data every STORE_VAR iterations
      if(elapsed_time >  print_iter * PRINT_INTERVAL) { 
        accumulate_data(&store_buffer, print_iter % SAVE_HISTORY, n_of_bodies, &system_status);
        print_iter++;
        // write the data to disk every when the buffer is full
        if( print_iter % SAVE_HISTORY == 0) {
          ret = write_data_to_disk(&store_buffer, n_of_bodies, print_iter, filename);
          if (ret == -1) goto cleanup;
        }
      }
    }
    
    compute_new_accelerations(system_status.data, system_status.mass, n_of_bodies,
                        count[rank], disp[rank], NEWTON);
    

    time_step_update(system_status.data, n_of_bodies, delta_t, count[rank], disp[rank], vel_range);

    MPI_Allgatherv(MPI_IN_PLACE, count[rank], MPI_DOUBLE,
                  system_status.data, count, disp, MPI_DOUBLE, comm);
  }
  
  if(rank == 0) printf("%.3f %.3f\n", delta_t, elapsed_time);
  
  if (rank == 0) free_store_buffer(&store_buffer);
  free(system_status.mass);
  free(system_status.data);
  free(count);
  free(disp);

  MPI_Finalize();

  return 0;

cleanup:
  if ( rank == 0 && ( NULL != store_buffer.mass || NULL != store_buffer.data) )
    free_store_buffer(&store_buffer);
  if ( NULL != system_status.mass) free(system_status.mass);
  if ( NULL != system_status.data) free(system_status.data);
  MPI_Abort(comm, 1);
}

