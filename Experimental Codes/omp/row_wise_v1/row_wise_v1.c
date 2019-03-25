
/****
    Author: Rayhan Shikder,
    email: shikderr@myumanitoba.ca
    MSc Student,
    Department of Computer Science,
    University of Manitoba, Winnipeg, MB, Canada
****/


#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"

//macros
#define ALPHABET_LENGTH 4
#define max(x,y) ((x)>(y)?(x):(y))


//global variables
char *string_A;
char *string_B;
char *unique_chars_C; //unique alphabets
int c_len;
short **P_Matrix;
short **DP_Results; //to store the DP values

//function prototypes
int get_index_of_character(char *str,char x, int len);
void print_matrix(short **x, int row, int col);
void calc_P_matrix_v1(short **P, char *b, int len_b, char *c, int len_c);
short lcs_yang_v1(short **DP, short **P, char *A, char *B, char *C, int m, int n, int u);
short lcs(short **DP, char *A, char *B, int m, int n);


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

void print_matrix(short **x, int row, int col)
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            printf("%d ",x[i][j]);
        }
        printf("\n");
    }
}

void calc_P_matrix_v1(short **P, char *b, int len_b, char *c, int len_c)
{
    #pragma omp parallel for
    for(int i=0;i<len_c;i++)
    {
        for(int j=2;j<len_b+1;j++)
        {
            if(b[j-2]==c[i]) //j-2 as b we assume here that b has a empty character in the beginning
            {
                P[i][j] = j-1;
            }
            else
            {
                P[i][j] = P[i][j-1];
            }
        }
    }
}


short lcs_yang_v1(short **DP, short **P, char *A, char *B, char *C, int m, int n, int u)
{
    for(int i=1;i<m+1;i++)
    {
        int c_i = get_index_of_character(C,A[i-1],u);
        #pragma omp parallel for schedule(static)
	for(int j=0;j<n+1;j++)
        {
            if(A[i-1]==B[j-1])
            {
                DP[i][j] = DP[i-1][j-1] + 1;
            }
            else if(P[c_i][j]==0)
            {
                DP[i][j] = max(DP[i-1][j], 0);
            }
            else
            {
                DP[i][j] = max(DP[i-1][j], DP[i-1][P[c_i][j]-1] + 1);
            }
        }
    }
    return DP[m][n];
}

short lcs(short **DP, char *A, char *B, int m, int n)
{
   // printf("%s %d \n%s %d\n",A,m,B,n );

    for(int i=1;i<(m+1);i++)
    {
        for(int j=1;j<(n+1);j++)
        {
            if(A[i-1] == B[j-1])
            {
                DP[i][j] = DP[i-1][j-1] + 1;
            }
            else
            {
                DP[i][j] = max(DP[i-1][j],DP[i][j-1]);
            }
        }
    }

    return DP[m][n];
}

int main(int argc, char *argv[])
{
    if(argc <= 1){
        printf("Error: No input file specified! Please specify the input file, and run again!\n");
        return 0;
    }
    printf("\nYour input file: %s \n",argv[1]);
    
    FILE *fp;
    int len_a,len_b;
    double start_time,stop_time,start_time_yang,stop_time_yang;

    fp = fopen(argv[1], "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    //printf("1 : %d %d %d\n", len_a, len_b, c_len );

    string_A = (char *)malloc((len_a+1) * sizeof(char *));
    string_B = (char *)malloc((len_b+1) * sizeof(char *));
    unique_chars_C = (char *)malloc((c_len+1) * sizeof(char *));

    fscanf(fp, "%s %s %s", string_A,string_B,unique_chars_C);
   // printf("Strings : %s\n %s\n %s\n", string_A, string_B, unique_chars_C );


    //allocate memory for DP Results
    DP_Results = (short **)malloc((len_a+1) * sizeof(short *));
    for(int k=0;k<len_a+1;k++)
    {
        DP_Results[k] = (short *)calloc((len_b+1), sizeof(short));
    }


    //allocate memory for P_Matrix array
    P_Matrix = (short **)malloc(c_len * sizeof(short *));
    for(int k=0;k<c_len;k++)
    {
        P_Matrix[k] = (short *)calloc((len_b+1), sizeof(short));
    }

    start_time = omp_get_wtime();
    printf("lcs is: %d\n",lcs(DP_Results,string_A,string_B,len_a,len_b));
    stop_time = omp_get_wtime();
    printf("time taken by normal algorithm is: %lf",stop_time-start_time);



    //resetting DP to zero values
    for(int k=0;k<len_a+1;k++)
    {
        //memset(DP_Results[k],0,len_b+1);
        for(int l=0;l<len_b+1;l++)
        {
            DP_Results[k][l]=0;
        }
    }
    printf("\n");

    //deallocate pointers
    free(P_Matrix);
    free(DP_Results);
    return 0;
}
