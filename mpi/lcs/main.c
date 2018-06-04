#include<stdio.h>
#include<string.h>
#include<stdlib.h>

char *C;
int c_len;

int max(int x,int y)
{
    return x>y?x:y;
}


int lcs(char *A, char *B, int m, int n)
{
    printf("%s %d \n%s %d\n",A,m,B,n );
    int DP[m+1][n+1];

    for(int i=0;i<(m+1);i++)
    {
        for(int j=0;j<(n+1);j++)
        {
            if(i==0 || j==0)
            {
                DP[i][j]=0;
            }
            else if(A[i-1] == B[j-1])
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
    char *A;// = "abcdef";
    char *B;// = "ace";
//    printf("arg 1: %s argv2: %s\n",argv[1],argv[2]);

    FILE *fp;
    int len_a,len_b;

    fp = fopen("/home/cs/grad/shikderr/lcs/data/test.txt", "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    printf("1 : %d %d %d\n", len_a, len_b, c_len );

    A = (char *)malloc((len_a+1) * sizeof(char *));
    B = (char *)malloc((len_b+1) * sizeof(char *));
    C = (char *)malloc((c_len+1) * sizeof(char *));

    fscanf(fp, "%s %s %s", A,B,C);
    printf("Strings are:\n A: %s\nB: %s\nUnique: %s\n", A, B, C );




    printf("lcs is: %d\n",lcs(A,B,strlen(A),strlen(B)));
    fclose(fp);
    return 0;
}
