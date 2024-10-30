#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "tools.h"

void set_initial_conditions(body_system *system, size_t n_of_bodies){
  unsigned int seed = time(NULL);

  double range_pos = (double) RAND_MAX / (GRID_MAX - GRID_MIN);
  double vel_size = 5 * (double) n_of_bodies / (GRID_MAX - GRID_MIN);
  
  int idx;
  for (int body = 0; body < n_of_bodies; body++){
    system->mass[body] = (double)rand_r(&seed) / RAND_MAX;
    
    idx = 3 * body;

    system->x_data[idx] = (double)rand_r(&seed) / range_pos + GRID_MIN;
    system->y_data[idx] = (double)rand_r(&seed) / range_pos + GRID_MIN;
 
    system->x_data[idx+1] = (2*(double)rand_r(&seed) / RAND_MAX - 1) * vel_size;
    system->y_data[idx+1] = (2*(double)rand_r(&seed) / RAND_MAX - 1) * vel_size;
 
 
    system->x_data[idx + 2] = 0.0;
    system->y_data[idx + 2] = 0.0;
  }

}


//update acceleration with actual conditions
// approximation of newton potential
void acceleration_newton_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){

  double dist_x, dist_y, radius, cubed_radius;
  
  int t_idx, s_idx;
  for (size_t target_body = 0; target_body < (count / 3); target_body++){
    
    t_idx = (int)(target_body * 3) + (int)first;
    double new_x_acc = 0.0;
    double new_y_acc = 0.0;
    double x_pos = data_x[t_idx]; // data[t_idx];
    double y_pos = data_y[t_idx];
    // To update the acc, we need to check every other body distances 
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      //distances between the two masses
      s_idx = (int) source_body * 3;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist_x = data_x[s_idx] - x_pos;
      dist_y = data_y[s_idx] - y_pos;

      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      cubed_radius = radius * radius * radius;
      
      new_x_acc += mass[source_body] * dist_x / cubed_radius;
      new_y_acc += mass[source_body] * dist_y / cubed_radius;
    }

    data_x[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
    data_y[t_idx + 2] = new_y_acc * GLOBAL_CONSTANT_G;
  }
  
}

//update acceleration with actual conditions
// approximation with yukawa potential
void acceleration_yukawa_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){

  double dist_x, dist_y, radius, cubed_radius;
  
  int t_idx, s_idx;
  for (size_t target_body = 0; target_body < (count / 3); target_body++){
    
    t_idx = (int)(target_body * 3) + (int)first;
    double new_x_acc = 0.0;
    double new_y_acc = 0.0;
    double x_pos = data_x[t_idx]; // data[t_idx];
    double y_pos = data_y[t_idx];
    // To update the acc, we need to check every other body distances 
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      //distances between the two masses
      s_idx = (int) source_body * 3;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist_x = data_x[s_idx] - x_pos;
      dist_y = data_y[s_idx] - y_pos;

      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      cubed_radius = radius * radius * radius;
      
      new_x_acc += mass[source_body] * dist_x * exp(- radius) * (radius + 1) / cubed_radius;
      new_y_acc += mass[source_body] * dist_x * exp(- radius) * (radius + 1) / cubed_radius;
    }

    data_x[t_idx + 2] = new_x_acc / mass[t_idx];
    data_y[t_idx + 2] = new_y_acc / mass[t_idx];
  }
  
}

//update acceleration with actual conditions
// approximation with gaussian potential
void acceleration_gaussian_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){

  double dist_x, dist_y, radius, cubed_radius;
  
  int t_idx, s_idx;
  for (size_t target_body = 0; target_body < (count / 3); target_body++){
    
    t_idx = (int)(target_body * 3) + (int)first;
    double new_x_acc = 0.0;
    double new_y_acc = 0.0;
    double x_pos = data_x[t_idx]; // data[t_idx];
    double y_pos = data_y[t_idx];
    // To update the acc, we need to check every other body distances 
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      //distances between the two masses
      s_idx = (int) source_body * 3;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist_x = data_x[s_idx] - x_pos;
      dist_y = data_y[s_idx] - y_pos;

      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      
      new_x_acc += mass[source_body] * 2 * dist_x * exp(- a*radius*radius);
      new_y_acc += mass[source_body] * 2 * dist_y * exp(- a*radius*radius);
    }

    data_x[t_idx + 2] = new_x_acc / mass[t_idx];
    data_y[t_idx + 2] = new_y_acc / mass[t_idx];
  }
  
}

//update acceleration with actual conditions
// approximation of 1/r^4 force + 0 under certain distance
void acceleration_forthpow_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){

  double dist_x, dist_y, radius, cubed_radius;
  
  int t_idx, s_idx;
  for (size_t target_body = 0; target_body < (count / 3); target_body++){
    
    t_idx = (int)(target_body * 3) + (int)first;
    double new_x_acc = 0.0;
    double new_y_acc = 0.0;
    double x_pos = data_x[t_idx]; // data[t_idx];
    double y_pos = data_y[t_idx];
    // To update the acc, we need to check every other body distances 
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      //distances between the two masses
      s_idx = (int) source_body * 3;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist_x = data_x[s_idx] - x_pos;
      dist_y = data_y[s_idx] - y_pos;
      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      
      if (radius > 0.1){
        forth_pow_radius = radius * radius * radius *radius;
        new_x_acc += mass[source_body] * dist_x / forth_pow_radius;
        new_y_acc += mass[source_body] * dist_y / forth_pow_radius;
      }
    }

    data_x[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
    data_y[t_idx + 2] = new_y_acc * GLOBAL_CONSTANT_G;
  }
  
}

//Time-step update, considering acceleration already updated
void time_step_update(double *data, size_t n_of_bodies, double delta_t, size_t count, size_t first){

  double delta_p, delta_v, position, velocity;
  int t_idx;

  
  for (size_t target_body = 0; target_body < (count / 3); target_body++){
    t_idx = (int)(target_body * 3) + (int) first;
    
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


int print_data(const char *filename, body_system **system, size_t n_of_bodies, int iter, int write_header) {
  FILE *file = fopen(filename, "a");
  if (file == NULL) {
    fprintf(stderr, "Error opening file for appending\n");
    return -1;
  }

  if (write_header) {
    fprintf(file, "iteration,body_index,mass,x_pos,x_vel,x_acc,y_pos,y_vel,y_acc\n");
  }
  
  for (int idx = 0; idx < SAVE_HISTORY; idx++){
    for (size_t i = 0; i < n_of_bodies; i++) {
      fprintf(file, "%d, %ld, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
        iter + idx, i, system[idx]->mass[i],
        system[idx]->x_data[3*i], system[idx]->x_data[3*i + 1], system[idx]->x_data[3*i + 2],
        system[idx]->y_data[3*i], system[idx]->y_data[3*i + 1], system[idx]->y_data[3*i + 2]);
    }
  }

  fclose(file); // Close file after each append
  
  return 0;
}
