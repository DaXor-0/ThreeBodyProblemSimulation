#include <stdio.h>
#include <time.h>
#include <math.h>

#include "tools_simulation.h"

double mass_range, vel_range, grid_max;

/**
 * @brief Predefined ranges for initializing body properties. Adjusts mass and velocity ranges based on system size.
 */
const ranges init_ranges[]={
  {4,     60.0,   15     , 200},
  {16,    60.0,   5      , 400},
  {64,    60.0,   5      , 800},
  {256,   60.0,   4      , 1600},
  {1024,  60.0,   2.5    , 3200},
  {4096,  0.01,   0.3125 , 6400},
  {16384, 0.002,  0.15625, 12800}
};

/**
 * @brief Retrieves mass, velocity, and position ranges for initializing a system of a given size.
 *
 * @param n_of_bodies Number of bodies to initialize.
 */
void get_init_ranges(size_t n_of_bodies){
  int num_ranges = sizeof(init_ranges) / sizeof(ranges);
  for (int idx = 0; idx < num_ranges - 1; idx++){
    if (n_of_bodies <= init_ranges[idx].bodies){
      mass_range = init_ranges[idx].mass_range;
      vel_range  = init_ranges[idx].velocity_range;
      grid_max   = init_ranges[idx].grid_size;
      return;
    }
  }
  mass_range = init_ranges[num_ranges - 1].mass_range;
  vel_range  = init_ranges[num_ranges - 1].velocity_range;
  grid_max   = init_ranges[num_ranges - 1].grid_size;
}


/**
 * @brief Randomly initializes the system with masses, positions, and velocities for each body.
 * 
 * @param system Pointer to the body system to initialize.
 * @param n_of_bodies Number of bodies in the system.
 */
void set_initial_conditions(body_system *system, size_t n_of_bodies){
  get_init_ranges(n_of_bodies);

  double pos_range = (double) RAND_MAX / (grid_max - GRID_MIN);
  
  int idx;
  for (int body = 0; body < n_of_bodies; body++){
    idx = 6 * body;
    system->mass[body]    = (double) rand() / RAND_MAX * mass_range;
    system->data[idx]     = (double) rand() / pos_range + GRID_MIN;
    system->data[idx + 1] = (2 * (double)rand() / RAND_MAX - 1) * vel_range;
    system->data[idx + 2] = 0.0;
    system->data[idx + 3] = (double) rand() / pos_range + GRID_MIN;
    system->data[idx + 4] = (2 * (double)rand() / RAND_MAX - 1) * vel_range;
    system->data[idx + 5] = 0.0;
  }
}


/**
 * @brief Adjusts position if a body moves out of grid bounds.
 *
 * @note
 * Not inteded to be used outside of the scope of this file
 */
static inline int check_out_of_bound(double *position){
  if(*position < GRID_MIN){
    *position = 2 * GRID_MIN - *position;
    return 1;
  } else if(*position > grid_max){
    *position = 2 * grid_max - *position;
    return -1;
  }
  return 0;
}


/**
 * @brief Updates positions and velocities of bodies based on current acceleration.
 * 
 * Is called by each rank on roughly n/comm_sz bodies, and updates
 * the position and velocity of the bodies
 * 
 * @param data Array of body data (positions, velocities, accelerations).
 * @param n_of_bodies Total number of bodies.
 * @param delta_t Time step for the update.
 *
 * @note
 * Acceleration must be calculated beforehand
 */
void time_step_update(double *data, size_t n_of_bodies, double delta_t){
  double new_x_pos, new_x_vel, new_y_pos, new_y_vel;
  int t_idx, out[2];
  for (size_t target_body = 0; target_body < n_of_bodies; target_body++){
    t_idx = (int)(target_body * 6);
    
    //evaluate postition and velocity and control if in boundary
    new_x_pos = data[t_idx]     + data[t_idx+1]*delta_t + 0.5*data[t_idx+2]*delta_t*delta_t;
    new_x_vel = data[t_idx + 1] + data[t_idx+2]*delta_t;
    new_y_pos = data[t_idx + 3] + data[t_idx+4]*delta_t + 0.5*data[t_idx+5]*delta_t*delta_t;
    new_y_vel = data[t_idx + 4] + data[t_idx+5]*delta_t;
    
    // Check if is in bound, if not warp it back and randomize new velocity
    out[0] = check_out_of_bound(&new_x_pos);
    if (out[0] != 0){
      new_x_vel = (double)rand() / RAND_MAX * vel_range * (double)out[0];
    }

    out[1] = check_out_of_bound(&new_y_pos);
    if (out[1] != 0){
      new_y_vel = (double)rand() / RAND_MAX * vel_range * (double)out[1];
    }

    //update position and velocity
    data[t_idx]   = new_x_pos;
    data[t_idx+1] = new_x_vel;
    data[t_idx+3] = new_y_pos;
    data[t_idx+4] = new_y_vel;
  }
}

/**
 * @brief Helper function to calculate acceleration based on the distance and mass for different interaction models.
 *
 * @note
 * Not inteded to be used outside of the scope of this file
 */
static inline void acceleration_update(double* acc, double mass, double dist, double radius, accel_t type){
  double cubed_radius;
  
  switch (type) {
    case NEWTON:
      cubed_radius = radius * radius * radius;
      if (cubed_radius < 10) cubed_radius = 10;
      *acc += mass * dist / cubed_radius;
      break;
    case YUKAWA:
      cubed_radius = radius * radius * radius;
      *acc += mass * dist * exp(- radius) * (radius + 1) / cubed_radius;
      break;
    case GAUSSIAN:
      *acc += mass * 2 * dist * exp(-radius*radius);
      break;
    case FORTH_POW:
      if (radius > FORTH_POW_THRESHOLD){
        radius = radius * radius * radius *radius;
        *acc += mass * dist / radius;
      }
      break;
    default:
      fprintf(stderr, "ERROR: <accel_t> not valid. Aborting...\n");
      break;
  }
}

/**
 * @brief Computes new accelerations for each body by summing contributions from all other bodies.
 * 
 * Is called by each rank on roughly n/comm_sz bodies, and calculates the accelleration
 * of each body by computing the distances of that body w.r.t. every other body on the system.
 * So each rank must loop through n/comm_sz bodies and each one of those iterations
 * goes through all the other n bodies.
 *
 * The acceleration is updated in the helper function which considers also what type of
 * model we're simulating.
 *
 * @param data Array of body data (positions, velocities, accelerations).
 * @param mass Array of masses for each body.
 * @param n_of_bodies Total number of bodies.
 * @param type Acceleration model type.
 * 
 * @return int Status code (0 if successful, -1 if error).
 */
int compute_new_accelerations(double* data, double* mass, size_t n_of_bodies, accel_t type){
  double dist_x, dist_y, radius;
  double t_x_pos, t_y_pos, new_x_acc, new_y_acc;

  for (size_t target_body = 0; target_body < n_of_bodies; target_body++){
    int t_idx = (int)(target_body * 6);
    
    t_x_pos = data[t_idx];
    t_y_pos = data[t_idx + 3];
    new_x_acc = 0.0;
    new_y_acc = 0.0;
    
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      int s_idx = (int) source_body * 6;

      if (s_idx == t_idx) continue;   // skip same body to avoid division by 0
      
      dist_x = data[s_idx] - t_x_pos;
      dist_y = data[s_idx + 3] - t_y_pos;

      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      
      acceleration_update(&new_x_acc, mass[source_body], dist_x, radius, type);
      acceleration_update(&new_y_acc, mass[source_body], dist_y, radius, type);
    }

    switch (type) {
      case FORTH_POW:
      case NEWTON:
        data[t_idx + 2] = new_x_acc * GLOBAL_CONSTANT_G;
        data[t_idx + 5] = new_y_acc * GLOBAL_CONSTANT_G;
        break;
      case GAUSSIAN:
      case YUKAWA:
        data[t_idx + 2] = new_x_acc / mass[target_body];
        data[t_idx + 5] = new_y_acc / mass[target_body];
        break;
      default:
        fprintf(stderr, "ERROR: <accel_t> not valid. Aborting...\n");
        return -1;
    }
  }
  return 0;
}

/**
 * @brief Calculates an appropriate time step for the simulation based on the maximum velocity.
 *
 * The function computes the largest velocity component across all bodies and derives a stable time step.
 *
 * @param data Array of body data containing positions, velocities, and accelerations.
 * @param n_of_bodies Total number of bodies in the system.
 * @return double Suggested time step for the simulation.
 */
double compute_new_delta_t(double* data, size_t n_of_bodies){
  size_t idx = 0;
  double this_velocity, max_velocity = 0.0;
  
  for(size_t body = 0; body < n_of_bodies; body++){
    idx = body * 6;
    this_velocity = data[idx + 1] * data[idx + 1] + data[idx + 4] * data[idx + 4];
    if (this_velocity > max_velocity) max_velocity = this_velocity;
  }

  return 0.0003 * (grid_max - GRID_MIN) / sqrt(max_velocity);
}
