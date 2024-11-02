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
 * @brief Allocates memory for the buffer to store mass and positional data of a body system.
 *
 * @param store_buffer Pointer to the body_system structure.
 * @param n_of_bodies The number of bodies in the system.
 * 
 * @return 0 on success, -1 if memory allocation fails.
 */
static inline int allocate_store_buffer(body_system* store_buffer, size_t n_of_bodies) {
  store_buffer->mass = (double*) malloc(TMP_BUF_SIZE * n_of_bodies * sizeof(double));
  store_buffer->data = (double*) malloc(6 * TMP_BUF_SIZE * n_of_bodies * sizeof(double));
  
  if (store_buffer->mass == NULL || store_buffer->data == NULL){
    fprintf(stderr, "Error: failed memory allocation for store_buffer. Aborting\n");
    return -1;
  }

  return 0;
}


/**
 * @brief Frees memory allocated for the buffer storing body system data.
 *
 * @param store_buffer Pointer to the body_system structure whose memory needs to be freed.
 */
static inline void free_store_buffer(body_system* store_buffer) {
  free(store_buffer->mass);
  free(store_buffer->data);
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
      system_status->data[idx * 6],      // x_pos
      system_status->data[idx * 6 + 1],  // x_vel
      system_status->data[idx * 6 + 2],  // x_acc
      system_status->data[idx * 6 + 3],  // y_pos
      system_status->data[idx * 6 + 4],  // y_vel
      system_status->data[idx * 6 + 5]); // y_acc
  fflush(stdout);
}


/**
 * @brief Copies the current state of the system to a buffer at a specified index.
 *
 * @param store_buffer Pointer to the buffer storing system data history.
 * @param store_buffer_index The index in the buffer where the data will be stored.
 * @param n_of_bodies The number of bodies in the system.
 * @param system_status Pointer to the current system state.
 */
static inline void accumulate_data(body_system* store_buffer, int buffer_index, size_t n_of_bodies, body_system* system_status) {
  double* mass_offset = store_buffer->mass + (ptrdiff_t)(buffer_index * n_of_bodies);
  double* data_offset = store_buffer->data + (ptrdiff_t)(buffer_index * n_of_bodies * 6);
  
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

/**
 * @brief Writes accumulated buffer data to disk.
 * Open and overwrites the file the first time it's called, then the
 * next times appends to the file.
 *
 * @param store_buffer Pointer to the buffer containing the stored data.
 * @param n_of_bodies The number of bodies in the system.
 * @param print_iter The current iteration number.
 * @param filename The name of the file to write the data to.
 * 
 * @return 0 on success, -1 on file I/O error.
 */
static inline int write_data_to_disk(body_system* store_buffer, size_t n_of_bodies, int print_iter, const char* filename) {
  FILE* file = fopen(filename, (print_iter == TMP_BUF_SIZE) ? "w" : "a");
  if (!file) {
    fprintf(stderr, "Failed to open file for writing");
    return -1;
  }
  
  // If file is opened for the first time, write a legenda of the values
  if (print_iter == TMP_BUF_SIZE){
    // printf("CISNSJANIOHGSIOANGIOASNGIONSAIONGISOANSIONGAKDSJNGJIASB");
    fprintf(file, "iter_number,body_id, mass, x_pos, x_vel, x_acc, y_pos, y_vel, y_acc\n");
  }

  for (int step = 0; step < TMP_BUF_SIZE; ++step) {
    int iter_idx;
    for (int i = 0; i < n_of_bodies; ++i) {
      iter_idx = step * n_of_bodies + i;
      fprintf(file, "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
        print_iter - TMP_BUF_SIZE + step + 1, i,
        store_buffer->mass[ iter_idx ],
        store_buffer->data[ iter_idx * 6],      // x_pos
        store_buffer->data[ iter_idx * 6 + 1],  // x_vel
        store_buffer->data[ iter_idx * 6 + 2],  // x_acc
        store_buffer->data[ iter_idx * 6 + 3],  // y_pos
        store_buffer->data[ iter_idx * 6 + 4],  // y_vel
        store_buffer->data[ iter_idx * 6 + 5]); // y_acc
    }
  }

  fclose(file);

  return 0;
}

#endif
