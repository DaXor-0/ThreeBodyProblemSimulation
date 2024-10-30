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
void acceleration_update(double* data, double* mass, size_t n_of_bodies){

  double dist, cubed_dist;
  
  int t_idx, s_idx;
  for (int target_body = 0; target_body < n_of_bodies; target_body++){
    
    t_idx = target_body * 3;
    double new_x_acc = 0.0;
    double x_pos = data[t_idx]; // data[t_idx];

    // To update the acc, we need to check every other body distances 
    for(int source_body = 0; source_body < n_of_bodies; source_body++){
      //distances between the two masses
      s_idx = source_body * 3;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist = data[s_idx] - x_pos;
      cubed_dist = fabs(dist * dist * dist);
      
      new_x_acc += mass[source_body] * dist / cubed_dist;
    }

    data[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
  }
  
}

//Time-step update, considering acceleration already updated
void time_step_update(double *data, size_t n_of_bodies ,double delta_t){
  
  double delta_p, delta_v, position, velocity;
  int t_idx;
  
  for (int target_body = 0; target_body < n_of_bodies; target_body++){
    t_idx = target_body * 3;
    
    //delta for position
    delta_p = data[t_idx+1]*delta_t + 0.5*data[t_idx+2]*delta_t*delta_t;
    //delta for velocity
    delta_v = data[t_idx+2]*delta_t;
    
    //evaluate postition and velocity and control if in boundary
    position = data[t_idx] + delta_p;
    velocity = data[t_idx+1] + delta_v;
    if(position < GRID_MIN){
      position = 2*GRID_MIN - position;
      velocity = -velocity;
    }
    else if(position > GRID_MAX){
      position = 2*GRID_MAX - position;
      velocity = -velocity;
    }
    
    //update position and velocity
    data[t_idx] = position;
    data[t_idx+1] = velocity;
  }

}
