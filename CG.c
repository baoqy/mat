#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
mat_t CG(mat_t A,mat_t b)//共轭梯度法(要求A正定)
{
    int m=A.m;
    int i,j,k;
    double a,a1,a2,b1;
    mat_t x=new_mat_vec(m);
    mat_t R=new_mat_vec(m);
    mat_t R1=new_mat_vec(m);
    mat_t r=new_mat_vec(m);
    mat_t rk=new_mat_vec(m);
    mat_t p=new_mat_vec(m);
    mat_t pk=new_mat_vec(m);
    mat_t xk=new_mat_vec(m);

    mat_mul(R,A,x);
    mat_sub(r,b,R);
    mat_copy(p,r);

    for(k=0;k<8;k++)
    {
        a1=dot_product(r,p);
        mat_mul(R,A,p);
        a2=dot_product(R,p);
        a=a1/a2;

        mat_scaler(R,p,a);
        mat_add(xk,x,R);

        mat_mul(R,A,p);
        mat_scaler(R1,R,a);
        mat_sub(rk,r,R1);

        a1=dot_product(rk,R);
        a2=dot_product(p,R);
        b1=-a1/a2;

        mat_scaler(R,p,b1);
        mat_add(pk,rk,R);

        mat_copy(r,rk);
        mat_copy(p,pk);
        mat_copy(x,xk);
        a=solution_judge(A,x,b,0.001);//Ax-b二范数精确度0.001
        if(a<0.001)
        {
            return x;
            break;
        }
    }

}
int main()
{
    int m,n;
    printf("行数与列数:\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_t b=new_mat_vec(m);
    mat_set_all(A);
    mat_set_all(b);
  
   mat_t x=CG(A,b);
  mat_print(x);
   free_mat(A);
   free_mat(b);
  free_mat(x);
   return 0;
}
