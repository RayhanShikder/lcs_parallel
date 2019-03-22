/*
 *  * Name: Rayhan Shikder
 *  * email: shikderr@cs.umanitoba.ca
 *  * student ID: 007833032
 *  * LCS_YANG Version 1
 **/

#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"
//macros
#define max(x,y) ((x)>(y)?(x):(y))


//global variables
char *string_A;
char *string_B;
char *unique_chars_C; //unique alphabets
int c_len;
short **DP_Results; //to store the DP values

//function prototypes
void print_matrix(int **x, int row, int col);
short lcs(short **DP, char *A, char *B, int m, int n);



void print_matrix(int **x, int row, int col)
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
    
    FILE *fp;
    int len_a,len_b;
    double start_time,stop_time;

    fp = fopen("/home/cs/grad/shikderr/lcs/data/32768.txt", "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    printf("1 : %d %d %d\n", len_a, len_b, c_len );

    string_A = (char *)malloc((len_a+1) * sizeof(char *));
    string_B = (char *)malloc((len_b+1) * sizeof(char *));
    unique_chars_C = (char *)malloc((c_len+1) * sizeof(char *));

    fscanf(fp, "%s %s %s", string_A,string_B,unique_chars_C);


    //allocate memory for DP Results
    DP_Results = (short **)malloc((len_a+1) * sizeof(short *));
    for(int k=0;k<len_a+1;k++)
    {
        DP_Results[k] = (short *)calloc((len_b+1), sizeof(short));
    }


    start_time = omp_get_wtime();
    printf("lcs is: %d\n",lcs(DP_Results,string_A,string_B,len_a,len_b));
    stop_time = omp_get_wtime();
    printf("time taken by normal algorithm is: %lf\n",stop_time-start_time);




    //deallocate pointers
    free(DP_Results);
    return 0;
}
