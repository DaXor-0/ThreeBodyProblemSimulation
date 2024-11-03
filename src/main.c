#include <stdio.h>
#include <time.h>

#include "tools_simulation.h"
#include "utils.h"

int main(int argc, char *argv[]){
  body_system system_status;
  double* store_buffer = NULL;
  int ret;
  size_t n_of_bodies, n_of_iter;
  char* filename;

  ret = set_inputs(argc, argv, &n_of_bodies, &n_of_iter, &filename);
  if(ret == -1) goto cleanup;

  store_buffer = (double*) malloc(2 * TMP_BUF_SIZE * n_of_bodies * sizeof(double));
  if ( store_buffer == NULL ) goto cleanup;

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
  set_initial_conditions(&system_status, n_of_bodies);
  print_status(&system_status, n_of_bodies);

  int print_iter = 0;
  double delta_t, elapsed_time = 0;
  for (int iter = 0; iter < n_of_iter; iter++){
    delta_t = compute_new_delta_t(system_status.vel, n_of_bodies);
    elapsed_time += delta_t;
      
    // accumulate data every PRINT_INTERVAL time elapsed
    if(elapsed_time >  print_iter * PRINT_INTERVAL) { 
      accumulate_data(store_buffer, print_iter % TMP_BUF_SIZE, n_of_bodies, &system_status);
      print_iter++;

      // write the data to disk every when the buffer is full
      if( print_iter % TMP_BUF_SIZE == 0) {
        ret = write_data_to_disk(store_buffer, system_status.mass, n_of_bodies, print_iter, filename);
        if (ret == -1) goto cleanup;
      }
    }

    compute_new_accelerations(system_status.mass, system_status.pos, system_status.acc,
                              n_of_bodies, NEWTON);

    time_step_update(system_status.pos, system_status.vel, system_status.acc,
                     n_of_bodies, delta_t);
  }

  printf("\nLast delta_t is: %.5f\n TOTAL ELAPSED (simulation) TIME:%.3f\n", delta_t, elapsed_time);
  
  free(store_buffer);
  free(system_status.mass);
  free(system_status.pos);
  free(system_status.vel);
  free(system_status.acc);

  return 0;

cleanup:
  if ( NULL != store_buffer)       free(store_buffer);
  if ( NULL != system_status.mass) free(system_status.mass);
  if ( NULL != system_status.pos)  free(system_status.pos);
  if ( NULL != system_status.vel)  free(system_status.vel);
  if ( NULL != system_status.acc)  free(system_status.acc);

  return -1;
}
