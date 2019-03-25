
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
#include "mpi.h"
#include<time.h>
//macros
#define ALPHABET_LENGTH 4
#define max(x,y) ((x)>(y)?(x):(y))


//global variables
char *string_A;
char *string_B;
char *unique_chars_C; //unique alphabets
int c_len;
short *P_Matrix;
short **DP_Results; //to store the DP values

//function prototypes
int get_index_of_character(char *str,char x, int len);
void print_matrix(short **x, int row, int col);
void calc_P_matrix_v2(short *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size);
short lcs_yang_v2(short **DP, short *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size);
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



void calc_P_matrix_v2(short *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size)
{
	char receive_array_for_scatter_c[chunk_size];
   short receive_array_for_scatter_p[chunk_size*(len_b+1)];
//Scatter the char array chunks by sending each process a particular chunk
MPI_Scatter(c, chunk_size, MPI_CHAR,&receive_array_for_scatter_c,chunk_size,MPI_CHAR, 0, MPI_COMM_WORLD);
//Scatter the char array chunks by sending each process a particular chunk
MPI_Scatter(P, chunk_size*(len_b+1), MPI_SHORT,&receive_array_for_scatter_p,chunk_size*(len_b+1),MPI_SHORT, 0, MPI_COMM_WORLD);
// Broadcast the whole b  array to everybody
MPI_Bcast(b, len_b, MPI_CHAR, 0, MPI_COMM_WORLD);

    
for(int i=0;i<chunk_size;i++)
    {
        for(int j=1;j<len_b+1;j++)
        {
            if(b[j-1]==receive_array_for_scatter_c[i])
            {
		receive_array_for_scatter_p[(i*(len_b+1))+j] = j;
            }
            else
            {
		receive_array_for_scatter_p[(i*(len_b+1))+j] = receive_array_for_scatter_p[(i*(len_b+1))+j-1];
            }
        }
    }
//now gather all the calculated values of P matrix in process 0
MPI_Gather(receive_array_for_scatter_p, chunk_size*(len_b+1), MPI_SHORT, P, chunk_size*(len_b+1), MPI_SHORT, 0, MPI_COMM_WORLD);
}

short lcs_yang_v2(short **DP, short *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size)
{
    MPI_Bcast(P, (u*(n+1)), MPI_SHORT, 0, MPI_COMM_WORLD);
    for(int i=1;i<m+1;i++)
    {
        int c_i = get_index_of_character(C,A[i-1],u);
	short dp_i_receive[chunk_size];
	MPI_Scatter(DP[i], chunk_size, MPI_SHORT,&dp_i_receive,chunk_size,MPI_SHORT, 0, MPI_COMM_WORLD);
	int start_id = (myrank * chunk_size);
        int end_id = (myrank * chunk_size) + chunk_size;

        int t,s;

        for(int j=start_id;j<end_id;j++)
        {
 	    if(j==start_id && myrank==0)j=j+1;
            t= (0-P[(c_i*(n+1))+j])<0;
            s= (0 - (DP[i-1][j] - (t*DP[i-1][P[(c_i*(n+1))+j]-1]) ));
            dp_i_receive[j-start_id] = ((t^1)||(s^0))*(DP[i-1][j]) + (!((t^1)||(s^0)))*(DP[i-1][P[(c_i*(n+1))+j]-1] + 1);
        }
	//now gather all the calculated values of P matrix in process 0
	MPI_Allgather(dp_i_receive, chunk_size, MPI_SHORT,DP[i], chunk_size, MPI_SHORT, MPI_COMM_WORLD);
    }
    return DP[m][n];
}


short lcs(short **DP, char *A, char *B, int m, int n)
{

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
    
    
	int my_rank;
    int num_procs;
    int chunk_size_p,chunk_size_dp;//chunk_size for P matrix and DP matrix
    int res;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //grab this process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); //grab the total num of processes

    FILE *fp;
    int len_a,len_b;
    double start_time,stop_time,start_time_yang,stop_time_yang;
    
    if(my_rank == 0)printf("\nYour input file: %s \n",argv[1]);
    fp = fopen(argv[1], "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    string_A = (char *)malloc((len_a+1) * sizeof(char *));
    string_B = (char *)malloc((len_b+1) * sizeof(char *));
    unique_chars_C = (char *)malloc((c_len+1) * sizeof(char *));

    fscanf(fp, "%s %s %s", string_A,string_B,unique_chars_C);

    chunk_size_p = (c_len/num_procs);
    chunk_size_dp = ((len_b+1)/num_procs);

    if(my_rank==0)
{
        printf("chunk_p: %d chunk_dp: %d procs: %d\n",chunk_size_p,chunk_size_dp,num_procs);
}
     
     DP_Results = (short **)malloc((len_a+1) * sizeof(short *));
    for(int k=0;k<len_a+1;k++)
    {
        DP_Results[k] = (short *)calloc((len_b+1), sizeof(short));
    }

    P_Matrix = (short *)malloc((c_len*(len_b+1)) * sizeof(short));

   
    for(int k=0;k<len_a+1;k++)
    {
        for(int l=0;l<len_b+1;l++)
        {
            DP_Results[k][l]=0;
        }
    }
    start_time_yang = MPI_Wtime();

    calc_P_matrix_v2(P_Matrix,string_B,len_b,unique_chars_C,c_len, my_rank, chunk_size_p);
    
    res = lcs_yang_v2(DP_Results,P_Matrix,string_A,string_B,unique_chars_C,len_a,len_b,c_len,my_rank, chunk_size_dp);

    stop_time_yang = MPI_Wtime();
    
    if(my_rank == 0)
        {
                printf("lcs_yang_v2 is: %d\n",res);
                printf("time taken for lcs_yang_v1 is: %lf\n",stop_time_yang-start_time_yang);
       
        }
    //deallocate pointers
    free(P_Matrix);
    free(DP_Results);
	
    // Shutdown MPI (important - don't forget!)
    MPI_Finalize();
    return 0;
}
