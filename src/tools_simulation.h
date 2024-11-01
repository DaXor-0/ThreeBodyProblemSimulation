#ifndef TOOLS_SIMULATION_H
#define TOOLS_SIMULATION_H

#include <stdlib.h>

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
 *
 * @example
 * ```
 * size_t TOTAL_BODY_COUNT = 10;
 * int COMM_SZ = 3;
 * size_t SPLIT_INDEX, EARLY_BODY_COUNT, LATE_BODY_COUNT;
 * COMPUTE_BODY_COUNT(TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX, EARLY_BODY_COUNT, LATE_BODY_COUNT);
 * ```
 */
#define COMPUTE_BODY_COUNT( TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX,   \
                                  EARLY_BODY_COUNT, LATE_BODY_COUNT ) \
  EARLY_BODY_COUNT = LATE_BODY_COUNT = TOTAL_BODY_COUNT / COMM_SZ;    \
  SPLIT_INDEX = TOTAL_BODY_COUNT % COMM_SZ;                           \
  if (0 != SPLIT_INDEX) {                                             \
    EARLY_BODY_COUNT = EARLY_BODY_COUNT + 1;                          \
  }                                                                   \

#define GRID_MAX 100
#define GRID_MIN 0

#define GLOBAL_CONSTANT_G   2
#define FORTH_POW_THRESHOLD 0.1

typedef enum{
  NEWTON,
  YUKAWA,
  GAUSSIAN,
  FORTH_POW,
}accel_t;

typedef struct{
  double* mass;
  double* data;
} body_system;

void set_initial_conditions(body_system *system, size_t n_of_bodies);

void time_step_update(double *data, size_t n_of_bodies ,double delta_t, size_t my_count, size_t my_first);

int compute_new_accelerations(double* data, double* mass, size_t n_of_bodies, size_t my_count, size_t my_first, accel_t type);

#endif
