#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
double d_Li(mat_t v,int k)//计算Hermite插值公式中用到的Li（xi）的导数,v的元素凉凉不相等
{
    int i;
    int n=v.m;
    double a=0;
    for(i=0;i<n;i++)
    {
        if(i==k)
        {;}
        else
        {
            a+=1/(v.mat[k][0]-v.mat[i][0]);
        }
    }
    return a;
}
double Li(mat_t v,int k,double c)//计算Hermite（同时也是Lagrange）插值公式的Li(x）(i=k时)
{
    int i;
    int n=v.m;
    double a1=1,a2=1;
    for(i=0;i<n;i++)
    {
        if(i==k)
        {;}
        else
        {
           a1*=(c-v.mat[i][0]);
           a2*=(v.mat[k][0]-v.mat[i][0]);
        }
    }
return a1/a2;
 }

double Hermite(mat_t x,mat_t y,mat_t d_y ,double c)//用Hermite插值公式计算新坐标为c时对应的函数值,d_y存放已知点的导数值
{
    int n=x.m;
    int i,j,k;
    double a1=0,a2=0;
    for(j=0;j<n;j++)
    {
        a1+=y.mat[j][0]*pow(Li(x,j,c),2)*(1-2*d_Li(x,j)*(c-x.mat[j][0]));
        a2+=d_y.mat[j][0]*pow(Li(x,j,c),2)*(c-x.mat[j][0]);
    }
    return a1+a2;
}
int main()
{
    int n;
    printf("方阵的阶数:\n");
    scanf("%d",&n);
    mat_t x=new_mat_vec(n);
    printf("横坐标向量为,");
    mat_set_all(x);
    mat_t y=new_mat_vec(n);
    printf("纵坐标向量,");
    mat_set_all(y);
    mat_t d_y=new_mat_vec(n);
    printf("导数向量,");
    mat_set_all(d_y);
    double c;
    printf("输入新坐标；\n");
    scanf("%lf",&c);
    printf("Hermite计算结果为%lf\n",Hermite(x,y,d_y,c));
            free_mat(x);
            free_mat(y);
            free_mat(d_y);
            return 0;
}
