#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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
  }

  return 0;  // Success
}

//update acceleration with actual conditions
int acceleration_update(planet **target, int n_of_bodies){

  double d_x, d_y, cubed_d_x, cubed_d_y;
  
  for (int body_i = 0; body_i < n_of_bodies; body_i++){
    target[body_i]->acc[0] = 0;
    target[body_i]->acc[1] = 0;
    for(int body_k = 0; body_k < n_of_bodies; body_k++){
      //distances between the two masses 
      d_x = target[body_k]->pos[0]-target[body_i]->pos[0];
      d_y = target[body_k]->pos[1]-target[body_i]->pos[1];
      cubed_d_x = abs(d_x*d_x*d_x);
      cubed_d_y = fabs(d_y*d_y*d_y);
      
      target[body]->acc[0] += target[body_2]->mass*d_x/cubed_d_x;
      target[body]->acc[1] += target[body_2]->mass*d_y/cubed_d_y;
    }
    target[body]->acc[0] = target[body]->acc[0]*GLOBAL_CONSTANT_G;
    target[body]->acc[1] = target[body]->acc[1]*GLOBAL_CONSTANT_G;
  }

  return 0;  // Success
}

//Time-step update
int time_step_update(planet **target, int n_of_bodies, double delta_t){
  
  double delta_x, delta_y, delta_vx, delta_vy;
  
  acceleration_update(target, n_of_bodies);
  
  for (int body = 0; body < n_of_bodies; body++){
    //delta for positions
    delta_x = target[body]->vel[0]*delta_t + 0.5*target[body]->acc[0]*delta_t*delta_t;
    delta_y = target[body]->vel[1]*delta_t + 0.5*target[body]->acc[1]*delta_t*delta_t;
    //delta for velocities
    delta_vx = target[body]->acc[0]*delta_t;
    delta_vy = target[body]->acc[1]*delta_t;
    
    //update of positions
    target[body]->pos[0] += delta_x;
    target[body]->pos[1] += delta_y;
    //update of velocities
    target[body]->vel[0] += delta_vx;
    target[body]->vel[1] += delta_vy;
  }

  return 0;  // Success
}
