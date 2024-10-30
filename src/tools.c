#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "tools.h"

int set_initial_conditions(planet **target, int n_of_bodies){
  unsigned int seed = time(NULL);
  
  for (int body = 0; body < n_of_bodies; body++){
    target[body]->mass = (int)rand_r(&seed) / RAND_MAX * 100;
    //target[body]->radius = (double)rand_r(&seed) / RAND_MAX * 100.0;
    target[body]->pos[0] = (double)rand_r(&seed) / RAND_MAX * 100.0;
    target[body]->pos[1] = (double)rand_r(&seed) / RAND_MAX * 100.0;
    target[body]->vel[0] = (double)rand_r(&seed) / RAND_MAX * 100.0;
    target[body]->vel[1] = (double)rand_r(&seed) / RAND_MAX * 100.0;
    target[body]->acc[0] = 0;
    target[body]->acc[1] = 0;
  }

  return 0;  // Success
}



