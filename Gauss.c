#include "mat.c"
#include <stdio.h>
#include <stdlib.h>
void Gauss(mat_t A, mat_t b) {
  if (A.m !=A.n) 
  {
    printf("输入矩阵非方阵！\n");
    exit(-1);
  } 
  else 
  {
    int m=A.m,n=A.n;
    int i, j;
    mat_t M=new_mat(m,n+1);
    mat_t G=new_mat(m,n);
    mat_t I=mat_scalar(n,1);
    mat_t Mx=new_mat(m,n+1);
    mat_t Gx=new_mat(m,n);
    mat_t L=new_mat(m,n);
    mat_t Lx=new_mat(m,n);
    mat_t Ly=new_mat(m,n);
    mat_t U=new_mat(m,n);
    for (i = 0; i < m; i++) //计算矩阵U的增广矩阵M
    {
      for (j = 0; j < n; j++) 
      {M.mat[i][j] = A.mat[i][j];}
      M.mat[i][n] = b.mat[i][0];
    }
    printf("方程的增广矩阵为:\n");
    mat_print(M);
    int k;
    for(k=0;k<n-1;k++)
    {//取单位阵I的第K列
       mat_t r=new_vec_get(I,k);
       mat_t rr=new_mat_row(n);
       mat_transpose(rr,r);
      //取增广矩阵第k列的下半部分v
      mat_t v=new_vec_get(M,k);
      for (i = 0; i <n; i++) 
      {
        if (i <=k) 
        {v.mat[i][0] = 0;}
      }
      mat_mul(G,v,rr);
      mat_sub(Gx,I,G);
      mat_add(Lx,I,G);
      if (k > 0) 
      {mat_mul(L, Ly, Lx);}
      if (k <(n -1)) 
      {
       mat_mul(Mx, Gx, M);
       mat_copy(U,Mx);//计算LU分解的U
      }
      mat_copy(M,Mx); //将计算中间过程的增广矩阵和算子迭代
      mat_copy(G,Gx);
      mat_copy(Ly,Lx);
    }
    printf("所用LU分解为L=:\n");
    mat_print(L);
    printf("所用LU分解U=:\n");
    mat_print(U);
    printf("消元后上三角增广矩阵为:\n");
    mat_print(Mx);
    printf("方程的解向量为:\n");
    b=new_vec_get(Mx,n);
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {A.mat[i][j]=Mx.mat[i][j];}
    }
    mat_print(mat_U_solve(A,b));
  }
}
int main() 
{
  int m, n;
  printf("输入矩阵行数m与列数n：\n");
  scanf("%d%d", &m, &n);
  mat_t U=new_mat(m,n);
  mat_set_all(U);
  mat_t b=new_mat_vec(m);
  mat_set_all(b);
  Gauss(U, b);
  free_mat(U);
  free_mat(b);
  return 0;
}
