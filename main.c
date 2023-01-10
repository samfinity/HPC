#include <stdio.h>
#include <mpi.h>


int main(int argc, char **argv){

   int nbp, rank;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nbp);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


   printf("I am the pro %d among %d\n", rank,nbp);

    MPI_Finalize();
}