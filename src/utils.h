#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

/**
 * @brief Writes the positions of bodies to an output file.
 *
 * This function saves the `x` and `y` coordinates of a set of bodies to the specified file. 
 * Each line in the output file represents the position of one body in the format: `<x> <y>`.
 *
 * @param output_file The name of the file where the results will be written.
 * @param x_pos Pointer to an array of `double` values representing the x-coordinates of the bodies.
 * @param y_pos Pointer to an array of `double` values representing the y-coordinates of the bodies.
 * @param num_bodies The number of bodies (length of the `x_pos` and `y_pos` arrays).
 * @return Returns `0` on success, or `-1` if the file could not be opened.
 * 
 * @note It is the caller's responsibility to ensure that `x_pos` and `y_pos` are valid arrays 
 *       of at least `num_bodies` elements. The function does not perform bounds checking.
 */
static inline int write_results(const char* output_file, const double *x_pos, const double *y_pos, double simulation_time, size_t num_bodies) {
  FILE *file = fopen(output_file, "w");
  if (!file) {
    perror("Failed to open file");
    return -1;
  }

  // Write each pair of coordinates (x, y) to the file
  for (size_t i = 0; i < num_bodies; ++i) {
    fprintf(file, "%lf %lf\n", x_pos[i], y_pos[i]);
  }
  fprintf(file, "%lf\n", simulation_time);

  fclose(file);
  return 0;
}

#endif
