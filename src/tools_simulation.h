#ifndef TOOLS_SIMULATION_H
#define TOOLS_SIMULATION_H

#include <stdlib.h>

// #define GRID_MAX 400 <- works for 4 bodies
// #define GRID_MAX 800 <- works for 16 bodies
// #define GRID_MAX 1600
#define GRID_MIN 0

#define GLOBAL_CONSTANT_G   100
#define FORTH_POW_THRESHOLD 0.1

/**
 * @brief Enum defining different types of acceleration laws.
 * 
 * - `NEWTON`: Inverse-square law of gravitation.
 * - `YUKAWA`: Modified Newtonian gravity with an exponential decay factor.
 * - `GAUSSIAN`: Gaussian function based on distance.
 * - `FORTH_POW`: Inverse-fourth power law with a minimum threshold (FORTH_POW_THRESHOLD).
 */
typedef enum{
  NEWTON,
  YUKAWA,
  GAUSSIAN,
  FORTH_POW,
}accel_t;

/**
 * @brief Structure representing a system of bodies, storing mass and associated positional and velocity data.
 * 
 * - `mass`: Pointer to an array of masses for each body.
 * - `data`: Pointer to an array containing position, velocity, and acceleration for each body.
 *           For each body, the data format is: [x_pos, x_vel, x_acc, y_pos, y_vel, y_acc].
 */
typedef struct{
  double* mass;
  double* data;
} body_system;

/**
 * @brief Struct representing ranges for initializing body system properties based on the number of bodies.
 * 
 * - `bodies`: Number of bodies in the system.
 * - `mass_range`: Maximum mass for the bodies in this range.
 * - `velocity_range`: Maximum initial velocity for bodies in this range.
 */
typedef struct{
  size_t bodies;
  const double mass_range;
  const double velocity_range;
  const double grid_size;
} ranges;


/**
 * @brief Predefined ranges for initializing body properties. Adjusts mass and velocity ranges based on system size.
 */
static const ranges init_ranges[]={
  {4,     60.0,   15     , 200},
  {16,    60.0,   5      , 400},
  {64,    60.0,   5      , 800},
  {256,   2.5,    1.25   , 1600},
  {1024,  0.05,   0.625  , 3200},
  {4096,  0.01,   0.3125 , 6400},
  {16384, 0.002,  0.15625, 12800}
};

/**
 * @brief Global variables to set up ranges of mass, velocity and grid size
 */
double mass_range, vel_range, grid_max;

/**
 * @brief Retrieves mass, velocity, and position ranges for initializing a system of a given size.
 *
 * @param n_of_bodies Number of bodies to initialize.
 */
static inline void get_init_ranges(size_t n_of_bodies){
  int num_ranges = sizeof(init_ranges) / sizeof(ranges);
  for (int idx = 0; idx < num_ranges - 1; idx++){
    if (n_of_bodies <= init_ranges[idx].bodies){
      mass_range = init_ranges[idx].mass_range;
      vel_range  = init_ranges[idx].velocity_range;
      grid_max   = init_ranges[idx].grid_size;
      return;
    }
  }
  mass_range = init_ranges[num_ranges - 1].mass_range;
  vel_range  = init_ranges[num_ranges - 1].velocity_range;
  grid_max   = init_ranges[num_ranges - 1].grid_size;
}

/**
 * @brief Distributes bodies across processes for load balancing.
 *
 * @param n_of_bodies [in] Total number of bodies.
 * @param comm_sz [in] Total number of processes.
 * @param split_index [out] Index of first late_body_count rank.
 * @param early_body_count [out] Bodies assigned to the first `split_index` processes.
 * @param late_body_count [out] Bodies assigned to remaining processes.
 */
static inline void compute_body_count(size_t n_of_bodies, int comm_sz, int *split_index, 
                        size_t *early_body_count, size_t *late_body_count) {
  *early_body_count = *late_body_count = n_of_bodies / (size_t)comm_sz;

  *split_index = n_of_bodies % comm_sz;

  // If there are any extra bodies, assign one to each of the first `split_index` processes
  if (*split_index != 0) {
    (*early_body_count)++;
  }
}

void set_initial_conditions(body_system *system, size_t n_of_bodies);

void time_step_update(double *data, size_t n_of_bodies, double delta_t, size_t my_count, size_t my_first);

int compute_new_accelerations(double* data, double* mass, size_t n_of_bodies, size_t my_count, size_t my_first, accel_t type);

double compute_new_delta_t(double* data, size_t n_of_bodies);

#endif
