#ifndef TOOLS_SIMULATION_H
#define TOOLS_SIMULATION_H

#include <stdlib.h>

// #define GRID_MAX 400 <- works for 4 bodies
// #define GRID_MAX 800 <- works for 16 bodies
#define GRID_MAX 1600
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

extern double mass_range, vel_range, grid_max;

extern const ranges init_ranges[];

/**
 * @brief Distributes bodies across processes for load balancing.
 *
 * This macro divides `TOTAL_BODY_COUNT` bodies among `COMM_SZ` processes, assigning the
 * first `SPLIT_INDEX` processes `EARLY_BODY_COUNT` bodies, and the remaining processes 
 * `LATE_BODY_COUNT` bodies. If `TOTAL_BODY_COUNT` is not divisible by `COMM_SZ`, the first 
 * `SPLIT_INDEX` processes handle one extra body.
 *
 * @param TOTAL_BODY_COUNT [in] Total number of bodies.
 * @param COMM_SZ [in] Total number of processes.
 * @param SPLIT_INDEX [out] Number of processes with one extra body.
 * @param EARLY_BODY_COUNT [out] Bodies assigned to the first `SPLIT_INDEX` processes.
 * @param LATE_BODY_COUNT [out] Bodies assigned to remaining processes.
 */
#define COMPUTE_BODY_COUNT( TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX,   \
                                  EARLY_BODY_COUNT, LATE_BODY_COUNT ) \
  EARLY_BODY_COUNT = LATE_BODY_COUNT = TOTAL_BODY_COUNT / COMM_SZ;    \
  SPLIT_INDEX = TOTAL_BODY_COUNT % COMM_SZ;                           \
  if (0 != SPLIT_INDEX) {                                             \
    EARLY_BODY_COUNT = EARLY_BODY_COUNT + 1;                          \
  }                                                                   \

void set_initial_conditions(body_system *system, size_t n_of_bodies);

void time_step_update(double *data, size_t n_of_bodies, double delta_t, size_t my_count, size_t my_first);

int compute_new_accelerations(double* data, double* mass, size_t n_of_bodies, size_t my_count, size_t my_first, accel_t type);

double compute_new_delta_t(double* data, size_t n_of_bodies);

#endif
