#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
mat_t Cholesky(mat_t A)
{
    int i,j,k,h;
    int m=A.m,n=A.n;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(A.mat[i][j]!=A.mat[j][i])
            {
                printf("矩阵非对称！\n");
                exit(-1);
            }                                  //疑问：正定与否如何判定？
        }
    }
    mat_t L=new_mat(m,n);
    L.mat[0][0]=sqrt(A.mat[0][0]);
        for(j=0;j<n;j++)
        {
            double sum1=0;
            for(k=0;k<=j-1;k++)
            {   
                sum1+=pow(L.mat[j][k],2);
            }
            L.mat[j][j]=sqrt(A.mat[j][j]-sum1);
            for(i=j+1;i<n;i++)
           {
            double sum2=0;
            for(h=0;h<=j-1;h++)
            {
                sum2+=(L.mat[i][h]*L.mat[j][h]);
            }
            L.mat[i][j]=(A.mat[i][j]-sum2)/L.mat[j][j]; 
        }
        }
        printf("对称正定矩阵的Cholesky分解L为:\n");
    mat_print(L);
    free_mat(L);
}
int main()
{
    int m,n;
    printf("行数m列数n：\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_set_all(A);
    Cholesky(A);
    free_mat(A);
    return 0;
}
