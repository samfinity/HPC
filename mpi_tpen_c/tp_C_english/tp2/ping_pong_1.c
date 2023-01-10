/*
 * ping_pong_1.c   --- TP2 : point to point communications :
 *                             the process 0 send a message to the process 1
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

  MPI_Init(&argc, &argv);
  /* ? */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  for(int i=0; i < nb_values;i++){
    values[i] = rand() / (RAND_MAX+1.);
  }

  if (rank == 0) {
    for (iter = 0; iter<nb_values; iter++)
      values[iter] = rand() / (RAND_MAX + 1.);
    /* ? */
    MPI_Send(values, nb_values,MPI_DOUBLE, 1,tag,MPI_COMM_WORLD);
  } else if(rank == 1) {
    /* ? */
    MPI_Recv(values, nb_values, MPI_DOUBLE,0,tag,MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
    printf("Me, process 1, received %d values (last = %g)"
           "from process 0.\n", nb_values, values[nb_values-1]);
  }

  /* ? */
  MPI_Finalize();
  return 0;
}
