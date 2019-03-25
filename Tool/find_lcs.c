#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include<time.h>
#include "omp.h"

//macros
#define ALPHABET_LENGTH 4
#define max(x,y) ((x)>(y)?(x):(y))


//global variables
char *string_A;
char *string_B;
char *unique_chars_C; //unique alphabets
int c_len;
int **P_Matrix;
int *DP_Results; //to store the DP values of current row
int *dp_prev_row; //to store the DP values of previous row

//function prototypes
int get_index_of_character(char *str,char x, int len);
void print_matrix(int **x, int row, int col);
void calc_P_matrix_v2(int **P, char *b, int len_b, char *c, int len_c);
int lcs_yang_v2(int *DP, int *prev_row, int **P, char *A, char *B, char *C, int m, int n, int u);
int lcs(int *DP, int *prev_row, char *A, char *B, int m, int n);


int get_index_of_character(char *str,char x, int len)
{
    for(int i=0;i<len;i++)
    {
        if(str[i]== x)
        {
            return i;
        }
    }
    return -1;//not found the character x in str
}


void calc_P_matrix_v2(int **P, char *b, int len_b, char *c, int len_c)
{
    #pragma omp parallel for
    for(int i=0;i<len_c;i++)
    {
        for(int j=1;j<len_b+1;j++)
        {
            if(b[j-1]==c[i])
            {
                P[i][j] = j;
            }
            else
            {
                P[i][j] = P[i][j-1];
            }
        }
    }
}

int lcs_yang_v2(int *DP, int *prev_row, int **P, char *A, char *B, char *C, int m, int n, int u)
{
    for(int i=1;i<m+1;i++)
    {
        int c_i = get_index_of_character(C,A[i-1],u);
        int t,s;
	
	   #pragma omp parallel for private(t,s) schedule(static)
        for(int j=0;j<n+1;j++)
        {
            t= (0-P[c_i][j])<0;
            s= (0 - (prev_row[j] - (t*prev_row[P[c_i][j]-1]) ));
            DP[j] = ((t^1)||(s^0))*(prev_row[j]) + (!((t^1)||(s^0)))*(prev_row[P[c_i][j]-1] + 1);
        }

        #pragma omp parallel for schedule(static)
        for(int j=0;j<n+1;j++){
            prev_row[j] = DP[j];
        }
    }
    return DP[n];
}



int lcs(int *DP, int *prev_row, char *A, char *B, int m, int n)
{
    for(int i=1;i<(m+1);i++)
    {
        for(int j=1;j<(n+1);j++)
        {
            if(A[i-1] == B[j-1])
            {
                DP[j] = prev_row[j-1] + 1;
            }
            else
            {
                DP[j] = max(prev_row[j],DP[j-1]);
            }
        }

        for(int j=0;j<n+1;j++){
            prev_row[j] = DP[j];
        }
    }

    return DP[n];
}


int main(int argc, char *argv[])
{

   printf("\nYour input file: %s \n",argv[1]);

    FILE *fp;
    int len_a,len_b;
    double start_time, stop_time;

    fp = fopen(argv[1], "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    
    printf("Length of sequence 1: %d bp\n", len_a);
    printf("Length of sequence 2: %d bp\n", len_b);


    // printf("\n##################################\n");
    printf("\n######## Parallel Results ########\n");
    // printf("##################################\n");
//looking at the number of available threads
    #pragma omp parallel 
  {          
      #pragma omp single
      {
          printf("Number of threads used: %d\n", omp_get_num_threads() );
      }
  }


    string_A = (char *)malloc((len_a+1) * sizeof(char *));
    string_B = (char *)malloc((len_b+1) * sizeof(char *));
    unique_chars_C = (char *)malloc((c_len+1) * sizeof(char *));

    fscanf(fp, "%s %s %s", string_A,string_B,unique_chars_C);

    //allocate memory for DP Results
    DP_Results = (int *)malloc((len_b+1) * sizeof(int));
    dp_prev_row = (int *)malloc((len_b+1) * sizeof(int));
    // for(int k=0;k<len_a+1;k++)
    // {
    //     DP_Results[k] = (int *)calloc((len_b+1), sizeof(int));
    // }


    //allocate memory for P_Matrix array
    P_Matrix = (int **)malloc(c_len * sizeof(int *));
    for(int k=0;k<c_len;k++)
    {
        P_Matrix[k] = (int *)calloc((len_b+1), sizeof(int));
    }

    calc_P_matrix_v2(P_Matrix,string_B,len_b,unique_chars_C,c_len);

    start_time = omp_get_wtime();
    calc_P_matrix_v2(P_Matrix,string_B,len_b,unique_chars_C,c_len);
    int res = lcs_yang_v2(DP_Results, dp_prev_row, P_Matrix,string_A,string_B,unique_chars_C,len_a,len_b,c_len);
    stop_time = omp_get_wtime();
    printf("Length of the LCS is: %d\n",res);
    printf("Total time taken: %lf seconds\n",stop_time - start_time);


    for(int l=0;l<len_b+1;l++)
    {
        DP_Results[l]=0;
        dp_prev_row[l]=0;
    }
    //resetting DP to zero values
    // for(int k=0;k<len_a+1;k++)
    // {
    //     for(int l=0;l<len_b+1;l++)
    //     {
    //         DP_Results[k][l]=0;
    //     }
    // }

    // printf("\n##################################\n");
    printf("\n######## Sequential Results ########\n");
    // printf("##################################\n");
    start_time = omp_get_wtime();
    int n_res = lcs(DP_Results, dp_prev_row, string_A, string_B, len_a, len_b);
    stop_time = omp_get_wtime();
    printf("Length of the LCS is: %d\n",n_res);
    printf("Total time taken: %lf seconds\n\n",stop_time - start_time);



    free(P_Matrix);
    free(DP_Results);
    return 0;
}
