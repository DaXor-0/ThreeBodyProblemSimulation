#ifndef TOOLS_H
#define TOOLS_H

#define COMPUTE_BODY_COUNT( TOTAL_BODY_COUNT, COMM_SZ, SPLIT_INDEX,        \
                                       EARLY_BODY_COUNT, LATE_BODY_COUNT ) \
    EARLY_BODY_COUNT = LATE_BODY_COUNT = TOTAL_BODY_COUNT / COMM_SZ;       \
    SPLIT_INDEX = TOTAL_BODY_COUNT % COMM_SZ;                              \
    if (0 != SPLIT_INDEX) {                                                \
        EARLY_BODY_COUNT = EARLY_BODY_COUNT + 1;                           \
    }                                                                      \

typedef struct{
  int mass;
  //double radius;
  double pos[2];
  double vel[2];
  double acc[2];
} planet;

// 10.24 10.40

int set_initial_conditions(planet **target, int n_of_bodies);

#endif
