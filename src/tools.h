#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>

#define SAVE_HISTORY 50

#define GLOBAL_CONSTANT_G 1
#define GRID_MAX 100
#define GRID_MIN 0

#define COMPUTE_BODY_COUNT( TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX,        \
                                       EARLY_BODY_COUNT, LATE_BODY_COUNT ) \
    EARLY_BODY_COUNT = LATE_BODY_COUNT = TOTAL_BODY_COUNT / COMM_SZ;       \
    SPLIT_INDEX = TOTAL_BODY_COUNT % COMM_SZ;                              \
    if (0 != SPLIT_INDEX) {                                                \
        EARLY_BODY_COUNT = EARLY_BODY_COUNT + 1;                           \
    }                                                                      \

typedef struct{
  double* mass;
  double* x_data; // [0]-> x_pos [1]->x_vel [2]->x_acc
  double* y_data; // [0]-> y_pos [1]->y_vel [2]->y_acc
} body_system;

void set_initial_conditions(body_system *system, size_t n_of_bodies);

void acceleration_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t my_count, size_t my_first);

void time_step_update(double *data, size_t n_of_bodies ,double delta_t, size_t my_count, size_t my_first);

int print_data(const char *filename, body_system **system, size_t n_of_bodies, int iter, int write_header);

#endif
