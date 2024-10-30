#include <stdio.h>
#include <mpi.h>

#include "tools.h"

int main(int argc, char** argv){
  size_t n_of_bodies = 3;
  size_t n_of_iter = 1000;

  switch (argc) {
    case 2:
      n_of_bodies = (size_t) strtol(argv[1], NULL, 10);
      break;
    case 3:
      n_of_bodies = (size_t) strtol(argv[1], NULL, 10);
      n_of_iter = (size_t) strtol(argv[2], NULL, 10);
      break;
    default:
      break;
  }

  MPI_Init(NULL, NULL);

  int rank, comm_sz;
  MPI_Comm_size(&comm_sz, MPI_COMM_WORLD);
  MPI_Comm_rank(&rank, MPI_COMM_WORLD);

  MPI_Finalize();

cleanup:
  return -1;


}
