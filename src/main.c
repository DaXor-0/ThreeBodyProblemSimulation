#include <stdio.h>
#include <time.h>

#include "tools_simulation.h"
#include "utils.h"

int main(int argc, char *argv[]){
  body_system store_buffer, system_status;
  int ret;
  size_t n_of_bodies, n_of_iter;
  char* filename;

  ret = set_inputs(argc, argv, &n_of_bodies, &n_of_iter, &filename);
  if(ret == -1) goto cleanup;

  ret = allocate_store_buffer(&store_buffer, n_of_bodies);
  if(ret == -1) goto cleanup;

  system_status.mass = (double *) malloc(n_of_bodies * sizeof(double));
  system_status.data = (double *) malloc(6 * n_of_bodies * sizeof(double));
  if (system_status.mass == NULL || system_status.data == NULL){
    goto cleanup;
  }

  srand(time(0));
  get_init_ranges(n_of_bodies);
  set_initial_conditions(&system_status, n_of_bodies);
  print_status(&system_status, n_of_bodies);

  int print_iter = 0;
  double delta_t, elapsed_time = 0;

  for (int iter = 0; iter < n_of_iter; iter++){
    delta_t = compute_new_delta_t(system_status.data, n_of_bodies);
    elapsed_time += delta_t;
      
    // accumulate data every STORE_VAR iterations
    if(elapsed_time >  print_iter * PRINT_INTERVAL) { 
      accumulate_data(&store_buffer, print_iter % TMP_BUF_SIZE, n_of_bodies, &system_status);
      print_iter++;
      // write the data to disk every when the buffer is full
      if( print_iter % TMP_BUF_SIZE == 0) {
        ret = write_data_to_disk(&store_buffer, n_of_bodies, print_iter, filename);
        if (ret == -1) goto cleanup;
      }
    }

    compute_new_accelerations(system_status.data, system_status.mass, n_of_bodies, NEWTON);

    time_step_update(system_status.data, n_of_bodies, delta_t);
  }
<<<<<<< HEAD

  printf("%.3f %.3f\n", delta_t, elapsed_time);
  free_store_buffer(&store_buffer);
=======
  
  if(rank == 0) printf("\nLast delta_t is: %.5f\n TOTAL ELAPSED TIME:%.3f\n", delta_t, elapsed_time);
  
  if (rank == 0) free_store_buffer(&store_buffer);
>>>>>>> main
  free(system_status.mass);
  free(system_status.data);

  return 0;

  cleanup:
    if ( NULL != store_buffer.mass || NULL != store_buffer.data )
    free_store_buffer(&store_buffer);
  if ( NULL != system_status.mass) free(system_status.mass);
  if ( NULL != system_status.data) free(system_status.data);
}
