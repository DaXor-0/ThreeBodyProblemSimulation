#ifndef TOOLS_SIMULATION_H
#define TOOLS_SIMULATION_H

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

#define MAX_LINE_LENGTH 1024

#define GRID_MIN 0
#define GRID_MAX 200
#define NUM_ITER 20000

#define GLOBAL_CONSTANT_G 100
#define DELTA_T_MOD 0.0003

/**
 * @brief Struct representing ranges for initializing body system properties based on the number of bodies.
 */
typedef struct {
    double x;
    double y;
    double v_x;
    double v_y;
    double a_x;
    double a_y;
    double mass;
} Body;

double compute_new_delta_t(Body* system, size_t n_of_bodies);

void compute_new_accelerations(Body* system, size_t n_of_bodies);

void update_pos_and_vel(Body* system, size_t n_of_bodies, double delta_t);

/**
 * @brief Read input file (argv[1]) to initialize the number of bodies.
 *
 * @param input_file The input file containing initial conditions.
 * @param num_bodies Pointer to the variable for storing the number of bodies.
 * 
 * @return 0 on success, -1 on error.
 */
static inline int count_bodies(const char* input_file, size_t* num_bodies) {
  FILE* file = fopen(input_file, "r");
  if (!file) {
    perror("Error opening file");
    return -1;
  }

  char line[MAX_LINE_LENGTH];
  size_t count = 0;

  fgets(line, MAX_LINE_LENGTH, file);                 // Skip the header line
  while (fgets(line, MAX_LINE_LENGTH, file)) {        // Count the remaining lines
    count++;
  }
  *num_bodies = count;
  
  fclose(file);
  return 0;
}


/**
 * @brief Read input file (argv[1]) to initialize the system.
 *
 * @param input_file The input file containing initial conditions.
 * @param system Pointer to the Body array of structures for storing conditions of the system.
 * @param num_bodies Number of bodies in the system.
 * 
 * @return 0 on success, -1 on error.
 */
static inline int parse_input(const char *input_file, Body* system, size_t num_bodies){
  FILE* file = fopen(input_file, "r");
  if (!file) {
    perror("Error opening file");
    return -1;
  }

  char line[MAX_LINE_LENGTH];
  fgets(line, MAX_LINE_LENGTH, file); // Skip header

  // Parse each line and populate the Body structure
  int i = 0;
  while (fgets(line, MAX_LINE_LENGTH, file) && i < num_bodies) {
    if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &system[i].x, &system[i].y, &system[i].v_x,
               &system[i].v_y, &system[i].a_x, &system[i].a_y, &system[i].mass) != 7) {
      fprintf(stderr, "Error parsing line %d\n", i + 1);
      fclose(file);
      return -1;
    }
    i++;
  }

  fclose(file);
  return 0;
}


/**
 * @brief Prints to screen the current status of each body in the system.
 *
 * @param system_status Pointer to the body_system structure containing current system status.
 * @param num_bodies The number of bodies in the system.
 */
static inline void print_bodies(Body* bodies, int num_bodies) {
  for (size_t i = 0; i < num_bodies; i++) {
    printf("Body %ld: x=%.2f, y=%.2f, v_x=%.2f, v_y=%.2f, a_x=%.2f, a_y=%.2f, mass=%.2f\n",
           i, bodies[i].x, bodies[i].y, bodies[i].v_x, bodies[i].v_y, bodies[i].a_x, bodies[i].a_y, bodies[i].mass);
  }
}

#endif
