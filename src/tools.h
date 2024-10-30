#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>

#define GLOBAL_CONSTANT_G 9.81
#define X_GRID_MAX 100
#define X_GRID_MIN 0
#define Y_GRID_MAX 100
#define Y_GRID_MIN 0

#define COMPUTE_BODY_COUNT( TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX,        \
                                       EARLY_BODY_COUNT, LATE_BODY_COUNT ) \
    EARLY_BODY_COUNT = LATE_BODY_COUNT = TOTAL_BODY_COUNT / COMM_SZ;       \
    SPLIT_INDEX = TOTAL_BODY_COUNT % COMM_SZ;                              \
    if (0 != SPLIT_INDEX) {                                                \
        EARLY_BODY_COUNT = EARLY_BODY_COUNT + 1;                           \
    }                                                                      \

typedef struct{
  size_t n_of_bodies;
  double mass;
  double* x_data[3]; // [0]-> x_pos [1]->x_vel [2]->x_acc
  double* y_data[3]; // [0]-> y_pos [1]->y_vel [2]->y_acc
} body_system;


int set_initial_conditions(body_system **system);

#endif
