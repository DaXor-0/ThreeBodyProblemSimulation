#ifndef TOOLS_SIMULATION_H
#define TOOLS_SIMULATION_H

#include <stdlib.h>

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
 * - `pos`: Pointer to an array containing position. For each body the data format is: [x_pos, y_pos].
 * - `vel`: Pointer to an array containing velocity. For each body the data format is: [x_vel, y_vel].
 * - `acc`: Pointer to an array containing acceleration. For each body the data format is: [x_acc, y_acc].
 */
typedef struct{
  double* mass;
  double* pos;
  double* vel;
  double* acc;
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

extern const ranges init_ranges[];

extern double mass_range, vel_range, grid_max;

void compute_body_count(size_t n_of_bodies, int comm_sz, int *counts, int *disp);

void get_init_ranges(size_t n_of_bodies);

void set_initial_conditions(body_system *system, size_t n_of_bodies);

void time_step_update(double *pos, double *vel, double *acc, size_t n_of_bodies,
                        double delta_t, size_t my_count, size_t my_first);

int compute_new_accelerations(double* mass, double* pos, double* acc, size_t n_of_bodies,
                                size_t my_count, size_t my_first, accel_t type);

double compute_partial_delta_t(double* vel, int my_count, int my_first);

#endif
