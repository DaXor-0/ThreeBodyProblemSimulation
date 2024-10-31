#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "tools.h"

void allocate_buffer(body_system* buffer, size_t n_of_bodies) {
  buffer->mass = (double*) malloc(SAVE_HISTORY * n_of_bodies * sizeof(double));
  buffer->data = (double*) malloc(6 * SAVE_HISTORY * n_of_bodies * sizeof(double));
}

void free_buffer(body_system* buffer) {
  free(buffer->mass);
  free(buffer->data);
}

void accumulate_data(body_system* buffer, int buffer_index, size_t n_of_bodies, body_system* system_status) {
  double* mass_offset = buffer->mass + (ptrdiff_t)(buffer_index * n_of_bodies);
  double* data_offset = buffer->data + (ptrdiff_t)(buffer_index * n_of_bodies * 6);
  
  for (int idx = 0; idx < n_of_bodies; ++idx) {
    mass_offset[idx]         = system_status->mass[idx];
    data_offset[6 * idx]     = system_status->data[6 * idx];     // x_pos
    data_offset[6 * idx + 1] = system_status->data[6 * idx + 1]; // x_vel
    data_offset[6 * idx + 2] = system_status->data[6 * idx + 2]; // x_acc
    data_offset[6 * idx + 3] = system_status->data[6 * idx + 3]; // y_pos
    data_offset[6 * idx + 4] = system_status->data[6 * idx + 4]; // y_vel
    data_offset[6 * idx + 5] = system_status->data[6 * idx + 5]; // y_acc
  }
}

void write_data_to_disk(body_system* buffer, size_t n_of_bodies, int true_iter, const char* filename) {
  FILE* file = fopen(filename, (true_iter + 1 == SAVE_HISTORY) ? "w" : "a");
  if (!file) {
    perror("Failed to open file for writing");
    return;
  }

  if (true_iter + 1 == SAVE_HISTORY){
    // printf("CISNSJANIOHGSIOANGIOASNGIONSAIONGISOANSIONGAKDSJNGJIASB");
    fprintf(file, "iter_number,body_id, mass, x_pos, x_vel, x_acc, y_pos, y_vel, y_acc\n");
  }

  for (int step = 0; step < SAVE_HISTORY; ++step) {
    int iter_idx;
    for (int i = 0; i < n_of_bodies; ++i) {
      iter_idx = step * n_of_bodies + i;
      fprintf(file, "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
        true_iter - SAVE_HISTORY + step + 1, i,
        buffer->mass[ iter_idx ],
        buffer->data[ iter_idx * 6],      // x_pos
        buffer->data[ iter_idx * 6 + 1],  // x_vel
        buffer->data[ iter_idx * 6 + 2],  // x_acc
        buffer->data[ iter_idx * 6 + 3],  // y_pos
        buffer->data[ iter_idx * 6 + 4],  // y_vel
        buffer->data[ iter_idx * 6 + 5]); // y_acc
      }
    }

    fclose(file);
}


void set_initial_conditions(body_system *system, size_t n_of_bodies){
  unsigned int seed = time(NULL);
  
  double range_pos = (double) RAND_MAX / (GRID_MAX - GRID_MIN);
  double vel_size = 5 * (double) n_of_bodies / (GRID_MAX - GRID_MIN);
  
  int idx;
  for (int body = 0; body < n_of_bodies; body++){
    system->mass[body] = (double)rand_r(&seed) / RAND_MAX;
    
    idx = 6 * body;

    system->data[idx]     = (double)rand_r(&seed) / range_pos + GRID_MIN;
    system->data[idx + 1] = (2*(double)rand_r(&seed) / RAND_MAX - 1) * vel_size;
    system->data[idx + 2] = 0.0;
    system->data[idx + 3] = (double)rand_r(&seed) / range_pos + GRID_MIN;
    system->data[idx + 4] = (2*(double)rand_r(&seed) / RAND_MAX - 1) * vel_size;
    system->data[idx + 5] = 0.0;
 
  }

}


//update acceleration with actual conditions
// approximation of newton potential
void acceleration_newton_update(double* data, double* mass, size_t n_of_bodies, size_t count, size_t first){

  double dist_x, dist_y, radius, cubed_radius;
  
  for (size_t target_body = 0; target_body < (count / 6); target_body++){
    int t_idx = (int)(target_body * 6) + (int)first;
    
    double t_x_pos = data[t_idx];
    double t_y_pos = data[t_idx + 3];
    double new_x_acc = 0.0;
    double new_y_acc = 0.0;
    
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      int s_idx = (int) source_body * 6;
      // skip same body to avoid division by 0
      if (s_idx == t_idx) continue;

      dist_x = data[s_idx] - t_x_pos;
      dist_y = data[s_idx + 3] - t_y_pos;

      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      cubed_radius = radius * radius * radius;
      
      new_x_acc += mass[source_body] * dist_x / cubed_radius;
      new_y_acc += mass[source_body] * dist_y / cubed_radius;
    }

    data[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
    data[t_idx + 5] = new_y_acc * GLOBAL_CONSTANT_G;
  }
  
}

//update acceleration with actual conditions
// approximation with yukawa potential
// void acceleration_yukawa_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){
//
//   double dist_x, dist_y, radius, cubed_radius;
//   
//   int t_idx, s_idx;
//   for (size_t target_body = 0; target_body < (count / 3); target_body++){
//     
//     t_idx = (int)(target_body * 3) + (int)first;
//     double new_x_acc = 0.0;
//     double new_y_acc = 0.0;
//     double x_pos = data_x[t_idx]; // data[t_idx];
//     double y_pos = data_y[t_idx];
//     // To update the acc, we need to check every other body distances 
//     for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
//       //distances between the two masses
//       s_idx = (int) source_body * 3;
//       // skip same body to avoid division by 0
//       if (s_idx == t_idx) continue;
//
//       dist_x = data_x[s_idx] - x_pos;
//       dist_y = data_y[s_idx] - y_pos;
//
//       radius = sqrt(dist_y * dist_y + dist_x * dist_x);
//       cubed_radius = radius * radius * radius;
//       
//       new_x_acc += mass[source_body] * dist_x * exp(- radius) * (radius + 1) / cubed_radius;
//       new_y_acc += mass[source_body] * dist_x * exp(- radius) * (radius + 1) / cubed_radius;
//     }
//
//     data_x[t_idx + 2] = new_x_acc / mass[t_idx];
//     data_y[t_idx + 2] = new_y_acc / mass[t_idx];
//   }
//   
// }
//
// //update acceleration with actual conditions
// // approximation with gaussian potential
// void acceleration_gaussian_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){
//
//   double dist_x, dist_y, radius;
//   
//   int t_idx, s_idx;
//   for (size_t target_body = 0; target_body < (count / 3); target_body++){
//     
//     t_idx = (int)(target_body * 3) + (int)first;
//     double new_x_acc = 0.0;
//     double new_y_acc = 0.0;
//     double x_pos = data_x[t_idx]; // data[t_idx];
//     double y_pos = data_y[t_idx];
//     // To update the acc, we need to check every other body distances 
//     for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
//       //distances between the two masses
//       s_idx = (int) source_body * 3;
//       // skip same body to avoid division by 0
//       if (s_idx == t_idx) continue;
//
//       dist_x = data_x[s_idx] - x_pos;
//       dist_y = data_y[s_idx] - y_pos;
//
//       radius = sqrt(dist_y * dist_y + dist_x * dist_x);
//       
//       new_x_acc += mass[source_body] * 2 * dist_x * exp(-radius*radius);
//       new_y_acc += mass[source_body] * 2 * dist_y * exp(-radius*radius);
//     }
//
//     data_x[t_idx + 2] = new_x_acc / mass[t_idx];
//     data_y[t_idx + 2] = new_y_acc / mass[t_idx];
//   }
//   
// }
//
// //update acceleration with actual conditions
// // approximation of 1/r^4 force + 0 under certain distance
// void acceleration_forthpow_update(double* data_x, double* data_y, double* mass, size_t n_of_bodies, size_t count, size_t first){
//
//   double dist_x, dist_y, radius, forth_pow_radius;
//   
//   int t_idx, s_idx;
//   for (size_t target_body = 0; target_body < (count / 3); target_body++){
//     
//     t_idx = (int)(target_body * 3) + (int)first;
//     double new_x_acc = 0.0;
//     double new_y_acc = 0.0;
//     double x_pos = data_x[t_idx]; // data[t_idx];
//     double y_pos = data_y[t_idx];
//     // To update the acc, we need to check every other body distances 
//     for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
//       //distances between the two masses
//       s_idx = (int) source_body * 3;
//       // skip same body to avoid division by 0
//       if (s_idx == t_idx) continue;
//
//       dist_x = data_x[s_idx] - x_pos;
//       dist_y = data_y[s_idx] - y_pos;
//       radius = sqrt(dist_y * dist_y + dist_x * dist_x);
//       
//       if (radius > 0.1){
//         forth_pow_radius = radius * radius * radius *radius;
//         new_x_acc += mass[source_body] * dist_x / forth_pow_radius;
//         new_y_acc += mass[source_body] * dist_y / forth_pow_radius;
//       }
//     }
//
//     data_x[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
//     data_y[t_idx + 2] = new_y_acc * GLOBAL_CONSTANT_G;
//   }
//   
// }

//Time-step update, considering acceleration already updated
void time_step_update(double *data, size_t n_of_bodies, double delta_t, size_t count, size_t first){

  double new_x_pos, new_x_vel, new_y_pos, new_y_vel;
  int t_idx;
  
  for (size_t target_body = 0; target_body < (count / 6); target_body++){
    t_idx = (int)(target_body * 6) + (int) first;
    
    //evaluate postition and velocity and control if in boundary
    new_x_pos = data[t_idx]     + data[t_idx+1]*delta_t + 0.5*data[t_idx+2]*delta_t*delta_t;
    new_x_vel = data[t_idx + 1] + data[t_idx+2]*delta_t;
    new_y_pos = data[t_idx + 3] + data[t_idx+4]*delta_t + 0.5*data[t_idx+5]*delta_t*delta_t;
    new_y_vel = data[t_idx + 4] + data[t_idx+5]*delta_t;

    if(new_x_pos < GRID_MIN){
      new_x_pos = 2 * GRID_MIN - new_x_pos;
      new_x_vel = - new_x_vel;
    }
    else if(new_x_pos > GRID_MAX){
      new_x_pos = 2 * GRID_MAX - new_x_pos;
      new_x_vel = - new_x_vel;
    }
    
    if(new_y_pos < GRID_MIN){
      new_y_pos = 2 * GRID_MIN - new_y_pos;
      new_y_vel = - new_y_vel;
    }
    else if(new_y_pos > GRID_MAX){
      new_y_pos = 2 * GRID_MAX - new_y_pos;
      new_y_vel = - new_y_vel;
    }
    //update position and velocity
    data[t_idx]   = new_x_pos;
    data[t_idx+1] = new_x_vel;
    data[t_idx+3] = new_x_pos;
    data[t_idx+4] = new_x_vel;
  }

}


