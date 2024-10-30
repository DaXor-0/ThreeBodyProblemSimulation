#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "tools.h"

int main(int argc, char** argv){

  MPI_Init(NULL, NULL);
  int rank, comm_sz;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &rank);

  if (argc < 3){
    if( rank == 0 ) fprintf(stderr, "Usage: %s <n_of_bodies> <n_of_iter>\n", argv[0]);
    goto cleanup;
  }

  size_t n_of_bodies, n_of_iter;
  char *endptr;
  n_of_bodies = (size_t) strtol(argv[1], &endptr, 10);
  if (*endptr != '\0' || n_of_bodies < 3) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number of bodies, it must be >= 3. Aborting...\n");
    goto cleanup;
  }

  n_of_iter = (size_t) strtol(argv[2], &endptr, 10);
  if (*endptr != '\0' || n_of_iter < 1000) {
    if( rank == 0 ) fprintf(stderr, "Error: Invalid number of iter, it must be >= 1000. Aborting...\n");
    goto cleanup;
  }
  
  int split_rank;
  size_t large_body_count, small_body_count;
  COMPUTE_BODY_COUNT(n_of_bodies, comm_sz, split_rank, large_body_count, small_body_count);

  MPI_Finalize();

  return 0;

cleanup:
  MPI_Abort(comm, 1);
<<<<<<< HEAD
=======
=======
int main(){







>>>>>>> test-leo
>>>>>>> main
}
