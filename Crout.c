#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
mat_t Crout(mat_t A,mat_t b)//n阶方阵A的Crout分解
{
    int i,j,k;
    int n=A.m;
    double sum;
    mat_t L=new_mat(n,n);
    mat_t U=new_mat(n,n);
    for(k=0;k<n;k++)
    {
        L.mat[k][0]=A.mat[k][0];
        U.mat[0][k]=(A.mat[0][k])/(L.mat[0][0]);//先求的L的第一列和U的第一行
    }
    for(i=1;i<n;i++)
    {
        sum=0;
        for(j=0;j<i;j++)
        {sum+=(L.mat[i][j])*(U.mat[j][i]);}
         L.mat[i][i]=A.mat[i][i]-sum;     //计算L.mat[k][k]用于计算U时作为除数    
      for(k=i;k<n;k++)
      {
         U.mat[i][k]=(A.mat[i][k]-sum)/(L.mat[i][i]);
         //计算L的下一列   
         sum=0;
         for(j=0;j<i;j++)
         {sum+=(L.mat[k][j])*(U.mat[j][i]);}
          L.mat[k][i]=A.mat[k][i]-sum;
      }
    }
  mat_t y= mat_L_solve(L,b);
  mat_t x=mat_U_solve(U,y);
  return x;
}
int main()
{
  int m,n;
  printf("矩阵行数与列数:\n");
  scanf("%d%d",&m,&n);
  mat_t A=new_mat(m,n);
  mat_set_all(A);
  mat_t b=new_mat_vec(m);
  mat_set_all(b);
  printf("解向量为:\n");
  mat_t x=Crout(A,b);
  mat_print(x);
  free_mat(A);
  free_mat(b);
  return 0;
  }
