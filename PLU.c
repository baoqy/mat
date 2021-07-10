#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
mat_t PLU(mat_t A,mat_t b)
{
   int m=A.m,n=A.n;
   int i,j,k,r;
   mat_t G=new_mat(m,n);
   mat_t Gk=new_mat(m,n);
   mat_t Ak=new_mat(m,n);
   mat_t R=new_mat(m,n);
   mat_t I=mat_scalar(m,1);
   mat_t bk=new_mat_vec(m);
       for(k=0;k<n-1;k++)
   {
      // 对第k列做主元选取
      mat_t Ac=new_mat_clone(A);
       for(i=k+1;i<m;i++)
       {
           if(Ac.mat[i][k]<0)              
           {Ac.mat[i][k]=0-Ac.mat[i][k];}
           if(Ac.mat[i][k]>Ac.mat[k][k])
           {
               Ac.mat[k][k]=Ac.mat[i][k];
               r=i;
           }
       }
      // 交换主元行
      mat_t A1=new_mat_clone(A);
      mat_t b1=new_mat_clone(b);
       for(j=0;j<n;j++)
       {
           A.mat[k][j]=A.mat[r][j];
           b.mat[k][0]=b.mat[r][0];
       }
       
       for(j=0;j<n;j++)
       {
           A.mat[r][j]=A1.mat[k][j];
           b.mat[r][0]=b1.mat[k][0];
       }
       //消元计算
     mat_t r=new_mat_row(n);
     r.mat[0][k]=1;
     mat_t v=new_vec_get(A,k);
     for(i=0;i<=k;i++)
     {
         v.mat[i][0]=0;
     }
     mat_mul(R,v,r);
     mat_scaler(G,R,1/A.mat[k][k]);
     mat_sub(Gk,I,G);
     mat_mul(Ak,Gk,A);mat_mul(bk,Gk,b);
     mat_copy(A,Ak);mat_copy(b,bk);
}
mat_print(mat_U_solve(A,b));
}
int main()
{
    int m,n;
    printf("矩阵行数与列数:\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_t b=new_mat_vec(m);
    mat_set_all(A);
    mat_set_all(b);
    PLU(A,b);
    free_mat(A);
    free_mat(b);
    return 0;
}
