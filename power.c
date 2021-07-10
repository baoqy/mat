#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
double power(mat_t A,double a)//A为非亏损的n阶实矩阵,乘幂法求A得按模最大特征值,精确度要求为a
{
    int m=A.m;
    int k;
    double max,maxk;
    mat_t y=new_mat_vec(m);
    //取向量元素全为1的向量z；
    mat_t z=new_mat_vec(m);
    for(k=0;k<m;k++)
    {
        z.mat[k][0]=1;
    }
    mat_t zk=new_mat_vec(m);
    for(k=0;k<m;k++)
    {
        z.mat[k][0]=1;
    }
    //计算部分
    for(k=0;k<100;k++)  //设定最大计算次数为100次
    {
        mat_mul(y,A,z);
        mat_t Y=mat_absolute(y);
        max=mat_max(Y);
    //精度判断
        if(k>0)
        {
         if((max-maxk)<a&&(max-maxk)>(-a))
        {
            return max;
            break;
        }
        }

        mat_scaler(zk,y,1/max);
        mat_copy(z,zk);
        maxk=max;
    }
}
int main()
{
    int m,n;
    double c;
    printf("矩阵行数与列数：\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    printf("设定精度:\n");
    scanf("%lf",&c);
    mat_set_all(A);
    printf("该精度要求下此矩阵按模最大特征值为%lf\n",power(A,c));
    free_mat(A);
  
    return 0;
}
