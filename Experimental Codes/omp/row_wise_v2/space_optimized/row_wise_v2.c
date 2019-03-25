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
int *DP_Results; //to store the DP values
int *dp_prev_row; // to store the previous row of DP results
//function prototypes
int get_index_of_character(char *str,char x, int len);
void print_matrix(int **x, int row, int col);
void calc_P_matrix_v2(int **P, char *b, int len_b, char *c, int len_c);
int lcs_yang_v2(int *DP, int *prev_dp, int **P, char *A, char *B, char *C, int m, int n, int u);
int lcs(int **DP, char *A, char *B, int m, int n);


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

int lcs_yang_v2(int *DP, int *prev_dp, int **P, char *A, char *B, char *C, int m, int n, int u)
{

    for(int i=1;i<m+1;i++)
    {
        int c_i = get_index_of_character(C,A[i-1],u);
        int t,s;

        #pragma omp parallel for private(t,s) schedule(static)
        for(int j=0;j<n+1;j++)
        {
            t= (0-P[c_i][j])<0;
            s= (0 - (prev_dp[j] - (t*prev_dp[P[c_i][j]-1]) ));

            DP[j] = ((t^1)||(s^0))*(prev_dp[j]) + (!((t^1)||(s^0)))*(prev_dp[P[c_i][j]-1] + 1);
        }

        #pragma omp parallel for schedule(static)
        for(int j=0;j<n+1;j++){
                prev_dp[j] = DP[j];
        }
    }
    return DP[n];
}


int lcs(int **DP, char *A, char *B, int m, int n)
{
//    printf("%s %d \n%s %d\n",A,m,B,n );

    //print_matrix(DP,m+1,n+1);

    for(int i=1;i<(m+1);i++)
    {
        for(int j=1;j<(n+1);j++)
        {
//            if(i==0 || j==0)
//            {
//                DP[i][j]=0;
//            }
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
    double start_time, stop_time;

    fp = fopen(argv[1], "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    printf("1 : %d %d %d\n", len_a, len_b, c_len );

    string_A = (char *)malloc((len_a+1) * sizeof(char *));
    string_B = (char *)malloc((len_b+1) * sizeof(char *));
    unique_chars_C = (char *)malloc((c_len+1) * sizeof(char *));

    fscanf(fp, "%s %s %s", string_A,string_B,unique_chars_C);
    //printf("Strings :  %s\n",  string_B );
    printf("length of string B: %zu \n C: %zu\n", strlen(string_B), strlen(unique_chars_C));
    printf("string C is: %s\n",unique_chars_C);
//allocate memory for DP Results

    DP_Results = (int *)malloc((len_b+1) * sizeof(int));
    dp_prev_row = (int *)malloc((len_b+1) * sizeof(int));
/*
 for(int k=0;k<len_a+1;k++)
    {
        DP_Results[k] = (int *)calloc((len_b+1), sizeof(int));
        if(DP_Results[k]==NULL)
        {
                printf("Can not malloc anymore on %d \n",k);
                break;
        }
    }
*/
        //allocate memory for P_Matrix array
    P_Matrix = (int **)malloc(c_len * sizeof(int *));
    for(int k=0;k<c_len;k++)
    {
        P_Matrix[k] = (int *)calloc((len_b+1), sizeof(int));
    }


//    printf("initial DP_Results: \n");
//    print_matrix(DP_Results,len_a+1,len_b+1);

//    printf("lcs is: %d\n",lcs(DP_Results,string_A,string_B,len_a,len_b));
//    printf("after normal lcs, DP_Results: \n");
//    print_matrix(DP_Results,len_a+1,len_b+1);

   // calc_P_matrix_v2(P_Matrix,string_B,len_b,unique_chars_C,c_len);
//    printf("\n P_Matrix is: \n");
//    print_matrix(P_Matrix,len_c,len_b+1);

//resetting dp array
    for(int k=0;k<len_b+1;k++){
        DP_Results[k] = 0;
        dp_prev_row[k] = 0;
    }
  //  printf("\n");
//    print_matrix(DP_Results,len_a+1,len_b+1);
    start_time = omp_get_wtime();
    calc_P_matrix_v2(P_Matrix,string_B,len_b,unique_chars_C,c_len);
    int res = lcs_yang_v2(DP_Results, dp_prev_row, P_Matrix,string_A,string_B,unique_chars_C,len_a,len_b,c_len);
    //printf("lcs_yang_v2 is: %d\n",res);
    stop_time = omp_get_wtime();
    printf("lcs_yang_v2 is: %d\n",res);
    printf("total time taken: %lf\n",stop_time-start_time);
//    printf("final DP_Results: \n");
//    print_matrix(DP_Results,len_a+1,len_b+1);
    //deallocate pointers
    free(P_Matrix);
    free(DP_Results);
    return 0;
}
