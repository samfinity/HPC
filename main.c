#include <stdio.h>
#include <mpi.h>


int main(int argc, char **argv){

   int nbp, rank;
   int value;
   MPI_Status stat_msg;


   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nbp);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if(rank==4){
       value=4000;
       MPI_Send(&value, 1,MPI_INT,6,10, MPI_COMM_WORLD);
   }else if(rank==6){
       MPI_Recv(&value, 1,MPI_INT,4,10,MPI_COMM_WORLD,&stat_msg);
       printf("i proc %d recieved %d from the proc 2",rank, value);
   }


    MPI_Finalize();
}