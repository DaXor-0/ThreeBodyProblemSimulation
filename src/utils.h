#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stddef.h>

#include "tools_simulation.h"

#define TMP_BUF_SIZE     100
#define PRINT_INTERVAL 0.125


/**
 * @brief Parses command line arguments to initialize the number of bodies, iterations, and output filename.
 *
 * @param argc The argument count from the command line.
 * @param argv The argument values from the command line.
 * @param n_of_bodies Pointer to the variable for storing the number of bodies.
 * @param n_of_iter Pointer to the variable for storing the number of iterations.
 * @param filename Pointer to the filename string for output.
 * 
 * @return 0 on success, -1 on error (if inputs are invalid).
 */
static inline int set_inputs(int argc, char **argv, size_t *n_of_bodies, size_t *n_of_iter, char **filename){
  if (argc < 4){
    fprintf(stderr, "Error: using '%s' as <n_of_bodies> <n_of_iter> <filename>\n", argv[0]);
    return -1;
  }

  char *endptr;
  *n_of_bodies = (size_t) strtol(argv[1], &endptr, 10);
  if (*endptr != '\0' || *n_of_bodies < 2) {
    fprintf(stderr, "Error: Invalid number of bodies, it must be >= 3. Aborting...\n");
    return -1;
  }

  *n_of_iter = (size_t) strtol(argv[2], &endptr, 10);
  if (*endptr != '\0' || *n_of_iter < 100) {
    fprintf(stderr, "Error: Invalid number of iter, it must be >= 1000. Aborting...\n");
    return -1;
  } 
  
  *filename = argv[3];
  
  return 0;
}


/**
 * @brief Prints to screen the current status of each body in the system.
 *
 * @param system_status Pointer to the body_system structure containing current system status.
 * @param n_of_bodies The number of bodies in the system.
 */
static inline void print_status(body_system *system_status, size_t n_of_bodies){
  printf("iter_number,body_id, mass, x_pos, x_vel, x_acc, y_pos, y_vel, y_acc\n"); 
  for (int idx = 0; idx < n_of_bodies; idx++)
    printf("0,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", idx,
      system_status->mass[idx],
      system_status->pos[idx * 2],      // x_pos
      system_status->vel[idx * 2],      // x_vel
      system_status->acc[idx * 2],      // x_acc
      system_status->pos[idx * 2 + 1],  // y_pos
      system_status->vel[idx * 2 + 1],  // y_vel
      system_status->acc[idx * 2 + 1]); // y_acc
  fflush(stdout);
}


/**
 * @brief Copies the current system positions to a buffer at a specified index.
 *
 * @param store_buffer Pointer to the buffer storing system data history.
 * @param store_buffer_index The index in the buffer where the data will be stored.
 * @param n_of_bodies The number of bodies in the system.
 * @param system_status Pointer to the current system state.
 */
static inline void accumulate_data(double* store_buffer, int buffer_index,
                                     size_t n_of_bodies, body_system* system_status){
  double *offset = NULL;
  offset  = store_buffer + (ptrdiff_t)(buffer_index * n_of_bodies * 2);

  for (int idx = 0; idx < n_of_bodies; ++idx) {
    offset[2 * idx]     = system_status->pos[2 * idx];     // x_pos
    offset[2 * idx + 1] = system_status->pos[2 * idx + 1]; // y_pos
  }
}


/**
 * @brief Writes accumulated buffer data to disk.
 * Open and overwrites the file the first time it's called, then the
 * next times appends to the file.
 *
 * @param store_buffer Pointer to the buffer containing the stored positions.
 * @param mass Array of masses for each body.
 * @param n_of_bodies The number of bodies in the system.
 * @param print_iter The current iteration number.
 * @param filename The name of the file to write the data to.
 * 
 * @return 0 on success, -1 on file I/O error.
 */
static inline int write_data_to_disk(double* store_buffer, double* mass, size_t n_of_bodies,
                                       int print_iter, const char* filename) {
  FILE* file = fopen(filename, (print_iter == TMP_BUF_SIZE) ? "w" : "a");
  if (!file) {
    fprintf(stderr, "Failed to open file for writing");
    return -1;
  }
  
  // If file is opened for the first time, write a legenda of the values
  if (print_iter == TMP_BUF_SIZE){
    // printf("CISNSJANIOHGSIOANGIOASNGIONSAIONGISOANSIONGAKDSJNGJIASB");
    fprintf(file, "iter_number, body_id, mass, x_pos, y_pos\n");
  }

  for (int step = 0; step < TMP_BUF_SIZE; ++step) {
    int iter_idx;
    for (int i = 0; i < n_of_bodies; ++i) {
      iter_idx = step * n_of_bodies + i;
      fprintf(file, "%d,%d,%lf,%lf,%lf\n",
        print_iter - TMP_BUF_SIZE + step + 1, i,
        mass[i],
        store_buffer[iter_idx * 2],       // x_pos
        store_buffer[iter_idx * 2 + 1]);  // y_pos
    }
  }

  fclose(file);

  return 0;
}

#endif
