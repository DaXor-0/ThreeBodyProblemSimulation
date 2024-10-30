#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "tools.h"

int set_initial_conditions(body_system *system){
  unsigned int seed = time(NULL);
  
  int idx;
  for (int body = 0; body < system->n_of_bodies; body++){
    system->mass[body] = (double)rand_r(&seed) / RAND_MAX * 100;
    
    idx = 3 * body;

    for(int i = 0; i < 2; i++){
      system->x_data[idx + i] = (double)rand_r(&seed) / RAND_MAX * 100.0;
      system->y_data[idx + i] = (double)rand_r(&seed) / RAND_MAX * 100.0;
    }
    system->x_data[idx + 2] = 0.0;
    system->y_data[idx + 2] = 0.0;
  }

  return 0;  // Success
}

  // void acceleration_update(float* data, float* mass, size_t n_of_bodies);

//update acceleration with actual conditions
int acceleration_update(body_system *system, size_t n_of_bodies){

  double dist, cubed_dist;
  
  int t_idx, s_idx;
  for (int target_body = 0; target_body < n_of_bodies; target_body++){
    
    t_idx = target_body * 3;
    double new_x_acc = 0.0;
    double x_pos = system->x_data[t_idx]; // data[t_idx];

    // To update the acc, we need to check every other body distances 
    for(int source_body = 0; source_body < n_of_bodies; source_body++){
      //distances between the two masses
      s_idx = source_body * 3;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist = system->x_data[s_idx] - x_pos;
      cubed_dist = fabs(dist * dist * dist);
      
      new_x_acc += system->mass[source_body] * dist / cubed_dist;
    }

    system->x_data[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
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
