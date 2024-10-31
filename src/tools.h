#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>

#define SAVE_HISTORY 50
#define STORE_VAR 3
#define GLOBAL_CONSTANT_G 2
#define GRID_MIN 0
#define GRID_MAX 100

#define COMPUTE_BODY_COUNT( TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX,        \
                                       EARLY_BODY_COUNT, LATE_BODY_COUNT ) \
    EARLY_BODY_COUNT = LATE_BODY_COUNT = TOTAL_BODY_COUNT / COMM_SZ;       \
    SPLIT_INDEX = TOTAL_BODY_COUNT % COMM_SZ;                              \
    if (0 != SPLIT_INDEX) {                                                \
        EARLY_BODY_COUNT = EARLY_BODY_COUNT + 1;                           \
    }                                                                      \

// problema:
// facciamo 2 Allgather per iterazione
// soluzione:
// possiamo rendere x_data e y_data solo data, e poi ricontrollare per bene
// come aggiorniamo le cose
typedef struct{
  double* mass;
  double* data;
} body_system;

void allocate_buffer(body_system* buffer, size_t n_of_bodies);

void free_buffer(body_system* buffer);

void accumulate_data(body_system* buffer, int buffer_index, size_t n_of_bodies, body_system* system_status);

void write_data_to_disk(body_system* buffer, size_t n_of_bodies, int true_iter, const char * filename);

void set_initial_conditions(body_system *system, size_t n_of_bodies);

void acceleration_newton_update(double* data, double* mass, size_t n_of_bodies, size_t my_count, size_t my_first);

void time_step_update(double *data, size_t n_of_bodies ,double delta_t, size_t my_count, size_t my_first);

#endif
