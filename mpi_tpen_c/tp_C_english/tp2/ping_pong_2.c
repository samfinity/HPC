/*
 * ping_pong_2.c   --- TP2 : point-to-point communications : ping-pong
 *
 * Author          : Dimitri LECAS (CNRS/IDRIS - France) <dimitri.lecas@idris.fr>
 *
*/

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char *argv[]) {
  int rank,iter;
  int nb_values=1000;
  int tag=99;
  double values[nb_values];
  MPI_Status status;
  double time_begin,time_end;

  /* ? */

  if (rank == 0) {
    for (iter = 0; iter<nb_values; iter++)
      values[iter] = rand() / (RAND_MAX + 1.);
    time_begin=MPI_Wtime();
    /* ? */
    time_end=MPI_Wtime();
    printf("Me, process 0, sent and received %d values"
           "(last = %g) from process 1 in %f seconds.\n",
           nb_values, values[nb_values-1], time_end-time_begin);
    /* ? */
  }

  /* ? */
  return 0;
}
