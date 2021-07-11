#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
double in_power(mat_t A)//反幂法非奇异方阵A的按模最小特征值
{
    int i,k;
    int n=A.m;
    double mk;
    mat_t y=new_mat_vec(n);
    mat_t yk=new_mat_vec(n);
    mat_t z=new_mat_vec(n);
    mat_t zk=new_mat_vec(n);
    for(i=0;i<n;i++)
    {
        z.mat[i][0]=1;
    }
    mat_t L=L_LU(A);
    mat_t U=U_LU(A);
    for(k=0;k<100;k++)//进行100次迭代计算
    {
        mat_t y=mat_L_solve(L,z);
        mat_t yk=mat_U_solve(U,y);
        mat_copy(y,mat_absolute(yk));
        mk=mat_max(y);
        mat_scaler(zk,yk,1/mk);
        mat_copy(z,zk);
    }
    return mk;
}
int main()
{
    int n;
    printf("方阵阶数:\n");
    scanf("%d",&n);
    mat_t A=new_mat(n,n);
    mat_set_all(A);
    printf("方阵的按模最小特征值为%lf\n",in_power(A));
    free_mat(A);
    return 0;
}
