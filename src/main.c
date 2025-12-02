#include <stdio.h>

#include "tools_simulation.h"
#include "utils.h"


int main(int argc, char **argv){
  Body* system = NULL;
  
  // Check if correct number of command line input arguments are given
  if (argc != 3){
    fprintf(stderr, "Error: using '%s' as <input_file> <output_file>\n", argv[0]);
    goto cleanup;
  }

  // Count number of bodies of the system
  int ret;
  size_t num_bodies;
  ret = count_bodies(argv[1], &num_bodies);
  if(!ret) goto cleanup;

  // Allocate memory for system
  system = (Body*)malloc(num_bodies * sizeof(Body));
  if (!system) {
    perror("Error allocating memory");
    goto cleanup;
  }

  // Parse input file to set initial conditions
  ret = parse_input(argv[1], system, num_bodies);
  if(ret == -1) goto cleanup;
  
  // Main computation loop
  double delta_t, simulation_time = 0;
  for (int iter = 0; iter < NUM_ITER; iter++){
    // Calculate delta_t of iteration based on the fastest body in the simulation
    delta_t = compute_new_delta_t(system, num_bodies);
    simulation_time += delta_t;
    
    // For each element iterate through the other elements to compute new accelerations
    compute_new_accelerations(system, num_bodies);
    
    // Update new positions and velocities with the calculated accelerations
    update_pos_and_vel(system, num_bodies, delta_t);
  }
  
  // Allocate memory for x and y arrays
  double *x = malloc(num_bodies * sizeof(double));
  double *y = malloc(num_bodies * sizeof(double));
  
  ret = write_results(argv[2], x, y, num_bodies, simulation_time, real_time);

  free(system);

  return 0;

cleanup:
  if ( NULL != system )   free(system);
  return -1;
}
