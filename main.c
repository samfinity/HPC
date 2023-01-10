#include <stdio.h>
#include <mpi.h>


//int main(int argc, char **argv){

   //int nbp, rank;
   //int value;
   //MPI_Status stat_msg;


//
   //MPI_Init(&argc, &argv);
   //MPI_Comm_size(MPI_COMM_WORLD, &nbp);
   //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/*
    //point to point with combi de send recv
   if(rank==4){
       value=4000;
       MPI_Send(&value, 1,MPI_INT,6,10, MPI_COMM_WORLD);
   }else if(rank==6){
       MPI_Recv(&value, 1,MPI_INT,4,10,MPI_COMM_WORLD,&stat_msg);
       printf("i proc %d recieved %d from the proc 2",rank, value);
   }*/


    //simul send and recv avec MPI_Sendrecv
/*
    int msg,value;


 int num_proc=(rank+1)%2;
 msg = rank+1000;
 MPI_Sendrecv (&msg,1,MPI_INT,num_proc,10,&value,1,MPI_INT,
 num_proc,10,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

printf("I, process %d, I received %d from process %d.\n", rank,value,num_proc);

    MPI_Finalize();

    */

int main(int argc, char *argv[])
{
int my_rank, nb_procs;
MPI_Status status;
char s;
char r[20];

 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

 switch (my_rank){
 case 0 : s='h'; break ;
 case 1 : s='e'; break ;
 case 2 : s='l'; break ;
 case 3 : s='l'; break ;
 case 4 : s='o'; break ;
 }

 MPI_Send(&s, 1, MPI_CHAR, 0 /* destinataire */, 99 /* tag */, MPI_COMM_WORLD);

 if (my_rank == 0) {
 for (int i=0; i<nb_procs; i++){
 MPI_Recv(r+i, 1, MPI_CHAR, i /* emetteur */,
 99 /* tag */, MPI_COMM_WORLD, &status);
 }
 r[nb_procs] = 0;
 printf("%d a recu : %s\n", my_rank, r);
 }

 MPI_Finalize();
 return 0;
 }



//}