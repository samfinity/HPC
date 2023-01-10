/*
 * produit_matrices.c --- TP5 : matrix products
 *
 * Auteur          : Dimitri LECAS (CNRS/IDRIS - France) <dimitri.lecas@idris.fr>
 *
 *
 *
 * Remarks :
 * --------- 
 *
 *   * We would like to do the matrix product  C = A * B in parallel.
 *
 *   * We suppose the matrices are square and the size N is divide by 
 *     the number of processes Nprocs.
 *
 *   * The process 0 initializes the matrices A and B and distributes
 *     to other processes.
 *
 *   * The A distribution is done by horizontal slices.
 *     The B distribution is done by vertical slices.
 *
 *
 *   * Each process has a slice for A and B.
 *
 *   * Each process calculates the diagonal block of the matrix C that he have.
 *     The blocks calculations of the non-diagonal blocks de C needs 
 *     communicati!ons between other processes.
 *
 *   * In fact, the operation in this case is a matrix block product.
 */

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef USE_MKL
#include <mkl.h>
#endif

int tag=1000;

/* Allocate matrix and init to zero */
void matrix_allocate(double **mat, int linenbr, int colnbr) {
  int iter;
  (*mat) = (double *) malloc(linenbr*colnbr*sizeof(double));
  for (iter=0; iter<linenbr*colnbr; iter++) {
    (*mat)[iter] = 0.; }
}

/* Init matrix with random number */
void random_number(double *mat,int n) {
  int iterl,iterc;
  for(iterl=0; iterl<n; iterl++)
    for(iterc=0;iterc<n; iterc++)
      mat[iterl*n+iterc] =  rand() / (RAND_MAX + 1.);
}

/* Matrix product C = A*B */
void matmul(double *A, double *B, double *C,int nl, int nc, int nk) {
  int iterl,iterc,iterk;
  double sum;

#ifdef USE_MKL
  double alpha,beta;
  alpha = 1.0;
  beta = 0.0;
  mkl_set_num_threads(1);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              nl, nc, nk, alpha, A, nk, B, nc, beta, C, nc);
#else
  for(iterl=0;iterl<nl;iterl++) {
    for(iterc=0;iterc<nc;iterc++) {
      sum = 0;
      for(iterk=0;iterk<nk;iterk++) {
        sum += A[iterl*nk+iterk]*B[iterk*nc+iterc]; }
      C[iterl*nc+iterc] = sum;
    }
  }
#endif
}

int main(int argc, char *argv[]) {
  int rank,Nprocs,N,NL;
  double *A,*B,*C,*CC;
  double *AL,*BL,*CL,*TEMP;
  MPI_Datatype type_temp,type_slice;
  int double_size;
  MPI_Aint lbound,displacement;
  int iter;
  MPI_Status status;
  int previous_rank,following_rank;

  /* Init MPI */

  if (rank == 0) {
    FILE * fichier = fopen("matrix_products.data","r");
    fscanf(fichier, "%d", &N);
    fclose(fichier);
  }

  /* Process 0 broadcast N */

  /* N need to be divisible by Nprocs */
  if ( (N % Nprocs) == 0)
    NL = N/Nprocs;
  else {
    fprintf(stderr, "N is not divisible by Nprocs\n");
    /* We stop the program */

  }

  /* Process 0 init matrices A and B */
  if (rank ==0) {
    /* Dynamic allocation of matrices A, B & C */
    matrix_allocate(&A,N,N);matrix_allocate(&B,N,N);
    matrix_allocate(&C,N,N);matrix_allocate(&CC,N,N);

    /* Init A & B */
    random_number(A,N);
    random_number(B,N);

    /* Sequential computation of A*B */
    matmul(A,B,CC,N,N,N);
  }

  /* Dynamoc allocation of local matrices */
  matrix_allocate(&AL,NL,N);matrix_allocate(&BL,N,NL);
  matrix_allocate(&CL,N,NL);matrix_allocate(&TEMP,NL,N);

  /* Build datatype for 1 chunck of N lignes & NL colonnes */

  /* Processus 0 distribute in AL the horizontal slices of A */

  /* Processus 0 distribute in BL the vertical slices of B */

  /* Computation of the diagonal blocks */
  matmul(AL,BL,&(CL[rank*NL*NL]),NL,NL,N);

  /* Computation for none diagonal blocks */

  /* First algorithm */
  for (iter=0; iter<Nprocs; iter++) {
    // Each process sned his AL slice to process k
    // and receives in TEMP the AL slice of process k
    if (rank != iter) {
    

    // Each process calculates his blocks above or below the diagonal block
    matmul(TEMP,BL,&(CL[iter*NL*NL]),NL,NL,N);
    }
  }

  /* Second algorithm */
  //previous_rank = (Nprocs+rank-1)%Nprocs;
  //following_rank = (rank+1)%Nprocs;
  //for (iter=1; iter<Nprocs; iter++) {
    // Each process sned his AL slice to his previous process
    // and receives the AL slice from the next process ( the values of AL changed)

    // Each process calculates his block above or below the diagonal block
    //matmul(AL,BL,&(CL[( (rank+iter)%Nprocs )*NL*NL]),NL,NL,N);
  //}

  /* The process 0 gather all CL slices from each processes to form the C matrix */

  /* Deallocate locals arrays */
  free(AL); free(BL); free(CL); free(TEMP);

  /* Verification of the results */
  if (rank == 0) {
    double Emax=0;
    for(iter=0; iter<N*N; iter++) {
      if (Emax < fabs(C[iter]-CC[iter])) {
        Emax = fabs(C[iter]-CC[iter]); }}
    free(A); free(B); free(C); free(CC);

    if (Emax < 1e-10)
      printf("Super !\nMatrix product A*B in parallel\n"
             "equal the sequential one\n");
    else
      printf("False result !\nMatrix product A*B in parallel\n"
             "different from the sequential one\n");
  }


  return 0;
}
