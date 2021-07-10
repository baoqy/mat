#include<stdio.h>  
#include<stdlib.h>
#include"mat.c"
mat_t Jacobi(mat_t A,mat_t b)
{
    int m=A.m,n=A.n;
    int i,j,k;
    mat_t X=new_mat_vec(m);
    for(i=0;i<m;i++)
    {
        X.mat[i][0]=b.mat[i][0]/A.mat[i][i];
    }
    mat_t Xk=new_mat_vec(m);
    mat_t v=new_mat_row(m);
    for(k=0;k<=7;k++)
    {
        for(i=0;i<n;i++)
        {
           double sum1=0;
          for(j=0;j<=i-1;j++)
          {
              sum1+=A.mat[i][j]*X.mat[j][0];
          }

          double sum2=0;
          for(j=i+1;j<n;j++)
          {
              sum2+=A.mat[i][j]*X.mat[j][0];
          }

          Xk.mat[i][0]=-(1/A.mat[i][i])*(sum1+sum2-b.mat[i][0]);
         }
        mat_copy(X,Xk);
        mat_transpose(v,X);
        printf("第%d次运算结果为:\n",k+1);
        mat_print(v);
    }
}
int main()
{
    int m,n;
    printf("矩阵行数与列数\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_t b=new_mat_vec(m);
    mat_set_all(A);
    mat_set_all(b);
    Jacobi(A,b);
    free_mat(A);
    free_mat(b);
    return 0;
}
