#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "tools_simulation.h"

double mass_range, vel_range, grid_max;

/**
 * @brief Predefined ranges for initializing body properties. Adjusts mass and velocity ranges based on system size.
 */
const ranges init_ranges[]={
  {4,     60.0,  15   , 200},
  {16,    60.0,  5    , 400},
  {64,    60.0,  5    , 800},
  {256,   60.0,  4    , 1600},  // don't print gif after this, it might require a loooot of time
  {1024,  30.0,  2.5  , 3200},
  {4096,  30.0,  2.5  , 6400},
  {16384, 30.0,  1.0  , 12800}
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

  double pos_range = (double) RAND_MAX / (grid_max - GRID_MIN);
  #pragma omp parallel
  {
    unsigned int seed = (unsigned int)(omp_get_thread_num() * 6 + 23);
    
    int idx;
    #pragma omp for private(idx)
    for (int body = 0; body < n_of_bodies; body++){
      idx = body * 2;
      system->mass[body]   = (double) rand_r(&seed) / RAND_MAX * mass_range;
      system->pos[idx]     = (double) rand_r(&seed) / pos_range + GRID_MIN;
      system->pos[idx + 1] = (double) rand_r(&seed) / pos_range + GRID_MIN;
      system->vel[idx]     = (2 * (double)rand_r(&seed) / RAND_MAX - 1) * vel_range;
      system->vel[idx + 1] = (2 * (double)rand_r(&seed) / RAND_MAX - 1) * vel_range;
    }
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
 * @param pos Array of body positions.
 * @param vel Array of body velocities.
 * @param acc Array of body accelerations.
 * @param n_of_bodies Total number of bodies.
 * @param delta_t Time step for the update.
 *
 * @note
 * Acceleration must be calculated beforehand
 */
void time_step_update(double *pos, double *vel, double *acc,
                      size_t n_of_bodies, double delta_t){
  double new_x_pos, new_x_vel, new_y_pos, new_y_vel;
  int t_idx, is_out[2];
  #pragma omp parallel
  {
    unsigned int seed = (unsigned int) (omp_get_thread_num() * 4 + 11);

    #pragma omp for private(new_x_pos, new_x_vel, new_y_pos, new_y_vel, t_idx, is_out)
    for (size_t target_body = 0; target_body < n_of_bodies; target_body++){
      t_idx = (int)(target_body * 2);
      
    //evaluate postition and velocity and control if in boundary
      new_x_pos = pos[t_idx]     + vel[t_idx]     * delta_t + 0.5 * acc[t_idx] * delta_t * delta_t;
      new_x_vel = vel[t_idx]     + acc[t_idx]     * delta_t;
      new_y_pos = pos[t_idx + 1] + vel[t_idx + 1] * delta_t + 0.5 * acc[t_idx + 1] * delta_t * delta_t;
      new_y_vel = vel[t_idx + 1] + acc[t_idx + 1] * delta_t;
      
      // Check if is in bound, if not warp it back and randomize new velocity
      is_out[0] = check_out_of_bound(&new_x_pos);
      if (is_out[0] != 0){
        new_x_vel = (double)rand_r(&seed) / RAND_MAX * vel_range * (double)is_out[0];
      }

      is_out[1] = check_out_of_bound(&new_y_pos);
      if (is_out[1] != 0){
        new_y_vel = (double)rand_r(&seed) / RAND_MAX * vel_range * (double)is_out[1];
      }

      //update position and velocity
      pos[t_idx]     = new_x_pos;
      pos[t_idx + 1] = new_y_pos;
      vel[t_idx]     = new_x_vel;
      vel[t_idx + 1] = new_y_vel;
    }
  }
}

/**
 * @brief Helper function to calculate acceleration based on the distance and mass for different interaction models.
 *
 * @note
 * Not inteded to be used outside of the scope of this file
 */
static inline void acceleration_update(double* acc, double mass, double dist, double radius,
                                       double cubed_radius, accel_t type){
  switch (type) {
    case NEWTON:
      if (cubed_radius < 10) cubed_radius = 10;
      *acc += mass * dist / cubed_radius;
      break;
    case YUKAWA:
      *acc += mass * dist * exp(- radius) * (radius + 1) / cubed_radius;
      break;
    case GAUSSIAN:
      *acc += mass * 2 * dist * exp(-radius*radius);
      break;
    case FORTH_POW:
      if (radius > FORTH_POW_THRESHOLD){
        radius = cubed_radius *radius;
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
 * @param mass Array of masses for each body.
 * @param pos Array of body positions.
 * @param acc Array of body accelerations.
 * @param n_of_bodies Total number of bodies.
 * @param type Acceleration model type.
 * 
 * @return int Status code (0 if successful, -1 if error).
 */
int compute_new_accelerations(double* mass, double* pos, double* acc,
                              size_t n_of_bodies, accel_t type){
  double dist_x, dist_y, radius, cubed_radius;
  double t_x_pos, t_y_pos, new_x_acc, new_y_acc;
  volatile int error_flag = 0;

  #pragma omp parallel for private(dist_x, dist_y, radius, cubed_radius, \
                                   t_x_pos, t_y_pos, new_x_acc, new_y_acc) \
                           shared(error_flag) schedule(static)
  for (size_t target_body = 0; target_body < n_of_bodies; target_body++){
    int t_idx = (int)(target_body * 2);
    
    t_x_pos = pos[t_idx];
    t_y_pos = pos[t_idx + 1];
    new_x_acc = 0.0;
    new_y_acc = 0.0;
    
    for(size_t source_body = 0; source_body < n_of_bodies; source_body++){
      int s_idx = (int) source_body * 2;

      if (s_idx == t_idx) continue;   // skip same body to avoid division by 0

      dist_x = pos[s_idx]     - t_x_pos;
      dist_y = pos[s_idx + 1] - t_y_pos;

      radius = sqrt(dist_y * dist_y + dist_x * dist_x);
      cubed_radius = radius * radius * radius;

      acceleration_update(&new_x_acc, mass[source_body], dist_x, radius, cubed_radius, type);
      acceleration_update(&new_y_acc, mass[source_body], dist_y, radius, cubed_radius, type);
    }

    switch (type) {
      case FORTH_POW:
      case NEWTON:
        acc[t_idx]     = new_x_acc * GLOBAL_CONSTANT_G;
        acc[t_idx + 1] = new_y_acc * GLOBAL_CONSTANT_G;
        break;
      case GAUSSIAN:
      case YUKAWA:
        acc[t_idx]     = new_x_acc / mass[target_body];
        acc[t_idx + 1] = new_y_acc / mass[target_body];
        break;
      default:
        #pragma omp atomic write
        error_flag = 1;
    }
  }

  if (error_flag){
    fprintf(stderr, "ERROR: <accel_t> not valid. Aborting...\n");
    return -1;
  }

  return 0;
}

/**
 * @brief Calculates an appropriate time step for the simulation based on the maximum velocity.
 *
 * The function computes the largest velocity component across all bodies and derives a stable time step.
 *
 * @param vel Array of body velocities.
 * @param n_of_bodies Total number of bodies in the system.
 * @return double Suggested time step for the simulation.
 */
double compute_new_delta_t(double* vel, size_t n_of_bodies){
  double this_velocity, max_velocity = 0.0;

  #pragma omp parallel for private(this_velocity) reduction(max:max_velocity)
  for(size_t idx = 0; idx < n_of_bodies; idx++){
    this_velocity = vel[idx * 2] * vel[idx * 2] + vel[idx * 2 + 1] * vel[idx * 2 + 1];
    if (this_velocity > max_velocity) max_velocity = this_velocity;
  }

  return 0.0003 * (grid_max - GRID_MIN) / sqrt(max_velocity);
}
