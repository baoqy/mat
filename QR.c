#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
void mat_QR(mat_t A,mat_t b)//n阶方阵A
{
    int n=A.n,m=A.m;
    int i,j;
    mat_t v=new_mat_vec(m);
    mat_t T=new_mat(n,n);
    mat_t bt=new_mat_vec(m);
    mat_t Ax=new_mat(n,n);
     for(j=0;j<n-1;j++)
     {
         for(i=j;i<n-1;i++)
         {
           mat_t v=new_vec_get(A,j);
           mat_t T=Givens(v,j,i+1);
           mat_mul(Ax,T,A);
          mat_mul(bt,T,b);
           mat_copy(b,bt);
           mat_copy(A,Ax);
         }
     }
     printf("Givens矩阵作用后得上三角矩阵A=:\n");
     mat_print(A);
     printf("Ax=b的解向量为x=:\n");
    mat_print(mat_U_solve(A,bt));
}
int main()
{
    int m,n;
    printf("矩阵行数与列数:\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_set_all(A);
    printf("输入矩阵为:\n");
    mat_print(A);
   mat_t b=new_mat_vec(n);
   printf("输入向量:\n");
    mat_set_all(b);
    mat_QR(A,b);
    free_mat(A);
    free_mat(b);
}
