#include <stdio.h>
#include <math.h>

#include "tools_simulation.h"


/**
 * @brief Calculates an appropriate time step for the simulation based on the maximum velocity.
 *
 * The function computes the largest velocity component across all bodies and derives a stable time step.
 *
 * @param system Array of struct of bodies representing the system.
 * @param num_bodies Total number of bodies in the system.
 *
 * @return double New time step for the simulation.
 */
double compute_new_delta_t(Body* system, size_t num_bodies){
  double this_velocity, max_velocity = 0.0;
  
  for(size_t i = 0; i < num_bodies; i++){
    this_velocity = sqrt(system[i].v_x * system[i].v_x + system[i].v_y * system[i].v_y);
    if (this_velocity > max_velocity) max_velocity = this_velocity;
  }

  return DELTA_T_MOD * (GRID_MAX - GRID_MIN) / max_velocity;
}


/**
 * @brief Computes new accelerations for each body by summing contributions from all other bodies.
 * 
 * @param system Array of struct of bodies representing the system.
 * @param num_bodies Total number of bodies.
 * 
 * @return int Status code (0 if successful, -1 if error).
 */
void compute_new_accelerations(Body* system, size_t num_bodies){
  double x_dist, y_dist, radius, cubed_radius;

  for (size_t i = 0; i < num_bodies; i++){
    double new_x_a = 0.0;
    double new_y_a = 0.0;
    
    for(size_t j = 0; j < num_bodies; j++){
      if (i == j) continue;   // skip same body to avoid division by 0

      x_dist = system[j].x - system[i].x;
      y_dist = system[j].y - system[i].y;

      radius = sqrt(y_dist * y_dist + x_dist * x_dist);
      
      cubed_radius = radius * radius * radius;
      if (cubed_radius < 10){
        cubed_radius = 10;
      }
      
      new_x_a += system[i].mass * x_dist / cubed_radius;
      new_y_a += system[i].mass * y_dist / cubed_radius;
    }
    system[i].a_x = new_x_a * GLOBAL_CONSTANT_G;
    system[i].a_y = new_y_a * GLOBAL_CONSTANT_G;
  }
}


/**
 * @brief Adjusts position if a body moves out of grid bounds
 * and return if out of bound or not.
 *
 * @param position Pointer to position to check.
 * 
 * @return 0 if out of bound, -1 if not.
 *
 * @note
 * Not inteded to be used outside of the scope of this file
 */
static inline int is_out_of_bound(double *position){
  if(*position < GRID_MIN){
    *position = 2 * GRID_MIN - *position;
    return 0;
  } else if(*position > GRID_MAX){
    *position = 2 * GRID_MAX - *position;
    return 0;
  }
  return -1;
}


/**
 * @brief Updates positions and velocities of bodies based on current acceleration.
 * 
 * @param system Array of struct of bodies representing the system.
 * @param num_bodies Total number of bodies.
 * @param delta_t Time step for the update.
 *
 * @note
 * Acceleration must be calculated beforehand
 */
void update_pos_and_vel(Body* system, size_t num_bodies, double delta_t){
  double new_x, new_v_x, new_y, new_v_y;

  for (size_t idx = 0; idx < num_bodies; idx++){

    //evaluate postition and velocity and control if in boundary
    new_x   = system[idx].x   + system[idx].v_x * delta_t + 0.5 * system[idx].a_x * delta_t * delta_t;
    new_v_x = system[idx].v_x + system[idx].a_x * delta_t;
    new_y   = system[idx].y   + system[idx].v_y * delta_t + 0.5 * system[idx].a_y * delta_t * delta_t;
    new_v_y = system[idx].v_y + system[idx].a_y * delta_t;
    
    // Check if new_x cooirdinate is in bound, if not warp it back and give new velocity
    if (!is_out_of_bound(&new_x)){
      new_v_x *= -1; // New velocity is opposite of previous velocity
    }

    // Do the same for the y position and velocity
    if (!is_out_of_bound(&new_y)){
      new_v_y *= -1;
    }

    //update position and velocity
    system[idx].x   = new_x;
    system[idx].y   = new_y;
    system[idx].v_x = new_v_x;
    system[idx].v_y = new_v_y;
  }
}
