#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>

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
