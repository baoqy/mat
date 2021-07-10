#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
mat_t S_D(mat_t A,mat_t b)//最速下降法求解方程Ax=b(要求A为正定矩阵！）;
{
   int m=A.m,n=A.n;
   int i,j,k;
   double a,a1,a2;
   mat_t R=new_mat_vec(m);
   mat_t R1=new_mat_vec(m);
   mat_t r=new_mat_vec(m);
   mat_t rk=new_mat_vec(m);
   mat_t x=new_mat_vec(m);
   mat_t xk=new_mat_clone(x);
for(k=0;k<=100;k++)           //进行100次迭代
{
   mat_mul(R,A,x);
   mat_sub(r,b,R);
   a1=dot_product(r,r);       //dot_product为计算两向量的内积函数
   mat_mul(R,A,r);
   a2=dot_product(R,r);
   a=a1/a2;
   mat_scaler(R,r,a);
   mat_add(xk,x,R);
   mat_copy(x,xk);
   mat_mul(R,A,r);
   mat_scaler(R1,R,a);
   mat_sub(rk,r,R1);
   mat_copy(r,rk);
}
  a=solution_judge(A,x,b,0.001);
   if(a<0.001)
  { 
      return x;
  }
 else
{
    printf("100次迭代后Ax-b二范数误差大于0.001");
}
}
int main()
{
    int m,n;
    printf("行数与列数：\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_t b=new_mat_vec(m);
    mat_set_all(A);
    mat_set_all(b);
   mat_t x=S_D(A,b);
    mat_print(x);
    free_mat(A);
    free_mat(b);
    return 0;
}
