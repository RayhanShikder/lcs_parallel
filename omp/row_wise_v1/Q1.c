
/*
 * Name: Rayhan Shikder
 * email: shikderr@cs.umanitoba.ca
 * student ID: 007833032
 * Question 1 (part 1)
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"

#define M 1024
#define N 1024
#define K 1024

#define PRINT_VECS 1  // flag so we can turn off printing when N is large
#define MAX_RAND 100  // max value of elements generated for array

void init_matrix(int **mat,int P,int Q,int zero_flag)
{
    int i,j;
    for(i=0;i<P;i++)
    {
        for(j=0;j<Q;j++)
        {
            mat[i][j] = zero_flag==0?(i+j):0;
        }
    }
}

void print_matrix(int **mat, int P, int Q)
{
    for(int i=0;i<P;i++)
    {
        for(int j=0;j<Q;j++)
        {
            printf("%d ",mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[])
{
    // Declare matrixes for data and results
    // Note: all threads share one copy of this
    int **A; // input matrix
    int **B; // input matrix
    int **C; //for storing the result

    double start_time; // use these for timing
    double stop_time;

    // Fill the vector and start the timer
    printf("N: %d\n", N);

    //allocate size for matrix A
    A = (int **) malloc(M * sizeof(int*));
    for(int i=0;i<M;i++)
    {
        A[i] = (int *) malloc(N * sizeof(int));
    }

    //allocate size for matrix B
    B = (int **) malloc(N * sizeof(int*));
    for(int i=0;i<N;i++)
    {
        B[i] = (int *) malloc(K * sizeof(int));
    }

    //allocate size for matrix C
    C = (int **) malloc(M * sizeof(int*));
    for(int i=0;i<M;i++)
    {
        C[i] = (int *) malloc(K * sizeof(int));
    }

    init_matrix(A, M,N,0);
    init_matrix(B, N,K,0);
    init_matrix(C, M,K,1); //initialize with zero values


   
    start_time = omp_get_wtime(); // can use this function to grab a
                                  // timestamp (in seconds)

    // matrix-matrix multiplication code
   
    #pragma omp parallel for
    for(int i=0;i<M;i++)
    {
       // #pragma omp parallel for
        for(int j=0;j<K;j++)
        {
            //#pragma omp parallel for
            for(int k=0;k<N;k++)
            {
                C[i][j]+=(A[i][k]*B[k][j]);
            }
        }
    }

    
    // Print result & timing info
    stop_time = omp_get_wtime();
    printf("Total time (sec): %f\n", stop_time - start_time);

    return EXIT_SUCCESS;;
}
